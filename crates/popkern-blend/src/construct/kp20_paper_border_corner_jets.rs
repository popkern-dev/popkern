//! Setback corner classification: per corner `c_k`, whether the two
//! pieces meeting there may share one KP20 corner block, and the
//! corner-compatible cross jets they take when they may (PK-398,
//! PK-399, PK-401). Nested in `kp20_paper_border_extract` for its
//! segment/calibration substrate; split out for the file budget.

use popkern_core::geometry::gregory::eval_sextic_row;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};

use super::section_path::SectionPath;
use super::{
    ChiPieceCalibration, CornerJet, Kp20SetbackArc, OrientedSegment, kp20_paper_border_corners,
    rim_corner_jet,
};
use crate::fillet_surface::FilletSurfaceResult;

/// A GAPPED corner — the corner point lies off both strips by more
/// than this fraction of the rim leg (1.6e-3 against a 1e-2 leg on
/// the blunt pentagon G2 sweeps before PK-399, whose contact curves
/// "missed" each other by 3.4e-3) — carries no trustworthy corner
/// data: the rim row bends to reach the off-surface point and its
/// tangent there leaves the strips' planes by ~9°. Such corners keep
/// the χ̃-calibrated ends with no shared corner block (the assembly's
/// overlap average and the rim reconcile remain the containment) and
/// are counted in the report. No fixture reaches it since PK-399.
const GAPPED_CORNER_LEG_FRACTION: f64 = 1.0e-2;
/// Minimum sine of the cap boundary's turn at a corner for a convex
/// tensor corner (10°). Cross-sectional rim ends (PK-401) make every
/// turn the face's apex angle, so this guards degenerate faces only.
const REFLEX_CORNER_SIN_FLOOR: f64 = 0.17;

/// The corner classification for one ring.
pub(super) struct CornerClassification {
    /// Per corner `c_k`: `.0` serves arc k's second half (ending at
    /// `c_k`, jet from arc k+1's rim), `.1` serves arc k+1's first
    /// half (jet from arc k's rim); `None` at a gapped or reflex
    /// corner (χ̃-calibrated ends, no shared block).
    pub jets: Vec<Option<(CornerJet, CornerJet)>>,
    /// Per corner: the cap boundary's signed turn (radians, positive
    /// towards the vertex), NaN when degenerate.
    pub turns: Vec<f64>,
    /// Corners classified gapped.
    pub gapped: usize,
    /// Corners classified flat or reflex.
    pub reflex: usize,
}

/// Classify every corner of the ring. At an exact corner the two
/// strips' tangent planes agree; each corner jet is the neighbouring
/// rim's own Hermite block rotated into the BISECTING plane of the
/// two strips' tangent planes (a no-op at exact corners; splits a
/// gapped corner's disagreement evenly — measured: the raw neighbour
/// tangent left the rim 7.1e-2 out of plane against a 5e-2 bar).
#[allow(clippy::too_many_arguments)]
pub(super) fn classify_corners(
    arcs: &[Kp20SetbackArc],
    fillet_surfaces: &[FilletSurfaceResult],
    paths: &[SectionPath],
    shims: &[OrientedSegment],
    calibration: &[ChiPieceCalibration; 2],
    vertex_position: Point3,
    trace: bool,
) -> CornerClassification {
    let n = arcs.len();
    let surface_of = |k: usize| &fillet_surfaces[paths[k].fillet_index].surface;
    let station_normal = |k: usize, row0: &[Point3; 7], t: f64| -> Vec3 {
        let p = eval_sextic_row(row0, t);
        let surface = surface_of(k);
        project_point_to_uv(surface, &p)
            .and_then(|(u, v)| surface.normal_at(u, v).normalized())
            .unwrap_or(Vec3::new(0.0, 0.0, 0.0))
    };
    let foot = |k: usize, p: Point3| -> f64 {
        let surface = surface_of(k);
        project_point_to_uv(surface, &p)
            .map(|(u, v)| surface.point_at(u, v).distance_to(&p))
            .unwrap_or(f64::INFINITY)
    };
    let mut gapped = 0_usize;
    let mut reflex = 0_usize;
    let mut turns = vec![f64::NAN; n];
    let jets = (0..n)
        .map(|k| {
            let next = (k + 1) % n;
            let a_row0 = arcs[k].section_halves[1];
            let b_row0 = arcs[next].section_halves[0];
            let corner = a_row0[6];
            let n_a = station_normal(k, &a_row0, 1.0);
            let n_b = station_normal(next, &b_row0, 0.0);
            let bisector = kp20_paper_border_corners::tangent_plane_bisector(n_a, n_b);
            let jet_a = kp20_paper_border_corners::jet_into_plane(
                rim_corner_jet(&b_row0, true),
                bisector,
            );
            let jet_b = kp20_paper_border_corners::jet_into_plane(
                rim_corner_jet(&a_row0, false),
                bisector,
            );
            let a_march =
                shims[k].calibrated_march(&a_row0, 6, &calibration[0].end, vertex_position);
            let b_march =
                shims[next].calibrated_march(&b_row0, 0, &calibration[1].start, vertex_position);
            turns[k] = kp20_paper_border_corners::cap_corner_turn(
                bisector,
                a_row0[6] - a_row0[5],
                a_march,
                jet_a.d1 * -1.0,
            )
            .unwrap_or(f64::NAN);
            if trace {
                eprintln!("  [corner] c_{k}: turn {:+.2}°", turns[k].to_degrees());
            }
            let leg = (a_row0[5] - a_row0[6]).length().max((b_row0[1] - b_row0[0]).length());
            let gap = foot(k, corner).max(foot(next, corner));
            if gap > GAPPED_CORNER_LEG_FRACTION * leg {
                gapped += 1;
                if trace {
                    eprintln!(
                        "  [corner] c_{k}: gapped ({gap:.2e} off the strips, leg {leg:.2e}) — calibrated ends, no shared block"
                    );
                }
                return None;
            }
            // A tensor-product cap corner needs a CONVEX boundary turn
            // at c (see `cap_corner_is_convex`); flat or reflex turns
            // keep the χ̃-calibrated ends as recorded containment.
            let convex = |rim_out: Vec3, station_march: Vec3, corner_march: Vec3| -> bool {
                kp20_paper_border_corners::cap_corner_is_convex(
                    bisector,
                    rim_out,
                    station_march,
                    corner_march,
                    REFLEX_CORNER_SIN_FLOOR,
                )
            };
            let a_convex = convex(a_row0[6] - a_row0[5], a_march, jet_a.d1 * -1.0);
            let b_convex = convex(b_row0[0] - b_row0[1], b_march, jet_b.d1 * -1.0);
            if !(a_convex && b_convex) {
                reflex += 1;
                if trace {
                    eprintln!(
                        "  [corner] c_{k}: reflex or flat cap corner (a {a_convex}, b {b_convex}) — calibrated ends, no shared block"
                    );
                }
                return None;
            }
            Some((jet_a, jet_b))
        })
        .collect();
    CornerClassification {
        jets,
        turns,
        gapped,
        reflex,
    }
}
