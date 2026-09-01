//! PK-090: corner-compatible two-piece KP20 border from a blend
//! boundary's 2n segments.
//!
//! The apex blend boundary alternates fillet-end arcs and face-edge
//! segments — two genuine surface pieces per KP20 sector, meeting at
//! the sector corner `c`. This extractor builds the paper cap's
//! [`Kp20PaperGuideBorder`] directly from those segments, with depth
//! rows sampled from the OWNING surfaces (exact ring marches for the
//! analytic cylinder/plane cases; first-order tangent-plane march
//! otherwise, recorded in the report). This replaces the invalid
//! smooth-cubic midpoint split, whose halves are not corner-compatible
//! at `c` and kink the assembled rim.
//!
//! Layout: sector `k` = (`FilletEnd[k]`, `FaceEdge[k]`); `c` is the
//! junction where a fillet end meets its trailing face segment, and
//! the spoke points `m` are the face→fillet junctions. Depth rows are
//! the cap-side C² prolongation of the ring (the WP1b convention)
//! expressed as the BB coefficients of the cross-derivative field
//! sampled at true rim points (PK-398): `row1 = row0 − chord`,
//! `row2 − 2 row1 + row0 = 1.5·(q2 − 2 q1 + q0)` with `q_d` ring
//! samples at depth `d` from each rim point.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::{Kp20PaperGuideBorder, Kp20PaperGuideBorderPiece};
use popkern_core::geometry::{Point3, SurfaceGeometry, Vec3};
use popkern_core::topology::BRep;

use crate::blend::{BlendBoundary, BoundarySegment};
use crate::fillet_surface::FilletSurfaceResult;

#[path = "kp20_paper_border_recovery.rs"]
mod recovery;
use recovery::recover_face_edge_surface;
#[path = "kp20_paper_border_section.rs"]
mod section;
#[path = "kp20_paper_border_section_fit.rs"]
mod section_fit;
#[path = "kp20_paper_border_setback.rs"]
mod setback;
use popkern_core::geometry::gregory::elevate_cubic_row_to_sextic;
pub use setback::kp20_paper_border_from_fillet_setback;

/// Extraction evidence for diagnostics and gating.
#[derive(Debug, Clone)]
pub struct Kp20PaperBorderExtractReport {
    /// Valence (number of sectors).
    pub n: usize,
    /// Uniform ring-march depth used for the cross rows.
    pub depth: f64,
    /// Worst gap between consecutive segment endpoints along the loop.
    pub max_loop_gap: f64,
    /// Segments whose owning surface had no exact ring march and fell
    /// back to a first-order tangent-plane march.
    pub first_order_marches: usize,
    /// `true` when the loop order was reversed to make the boundary
    /// wind counterclockwise about its mean normal.
    pub reversed_for_orientation: bool,
    /// Setback path only: the TRUE trim-intersection points per
    /// junction (before any corner pullback) — the natural face
    /// wedge apexes for downstream trimming/rendering.
    pub trim_corners: Vec<Point3>,
    /// Setback path only: per-fillet setback arc data for topology
    /// insertion (the full arc cubic and its shared corner points).
    pub setback_arcs: Vec<Kp20SetbackArc>,
    /// Face-edge segments whose recorded face handle was stale (the
    /// n-gon boundary is computed before the corner-face topology
    /// surgery) and whose owning surface was recovered by curve
    /// incidence instead.
    pub recovered_face_surfaces: usize,
    /// Setback path only: worst sample distance between the fitted
    /// degree-6 section rows and the strips' sections (PK-398).
    pub section_fit_residual: f64,
    /// Setback path only: arcs whose section fit was rejected and
    /// served as the elevated cubic interpolant.
    pub section_fit_fallbacks: usize,
    /// Setback path only: corners whose shared corner point lies off
    /// the strips by more than 1 % of the rim leg (the strips' contact
    /// curves miss each other), served with χ̃-calibrated ends and no
    /// shared corner block — the corner data there is not trustworthy.
    pub gapped_corners: usize,
    /// Setback path only: corners where the cap boundary's turn is
    /// flat or reflex (the neighbour rim leaves `c` on the far side of
    /// this rim, or within 10° of it), where a tensor-product corner
    /// would fold; served with χ̃-calibrated ends and no shared corner
    /// block (PK-399 containment — mixed-radius setback placement).
    pub reflex_corners: usize,
    /// Setback path only: per corner `c_s`, the eq (14) corner-block
    /// disagreement between the two pieces meeting there, as the
    /// assembly's overlap average will see it (PK-399). Machine zero
    /// when the shared corner block reconciled; a residual here is
    /// rim movement the assembly is about to apply.
    pub corner_block_mismatch: Vec<f64>,
}

#[path = "kp20_paper_border_corners.rs"]
mod kp20_paper_border_corners;

/// One fillet's setback arc, exported for BRep topology insertion:
/// the arc is the fillet face's new apex-side boundary, running from
/// `corner_prev` to `corner_next` (shared with the neighbouring
/// fillets' arcs).
#[derive(Debug, Clone)]
pub struct Kp20SetbackArc {
    /// Index into the input `fillet_surfaces` slice.
    pub fillet_index: usize,
    /// Full arc as a cubic Bezier, `corner_prev` → `corner_next`
    /// (the coarse arc used for crossing classification and topology
    /// splitting; the served rim follows the cap rows).
    pub cubic: [Point3; 4],
    /// The exact section spline halves `[corner_prev → m]`,
    /// `[m → corner_next]` (degree 6, C² at `m`) that the border rows
    /// carry (PK-398).
    pub section_halves: [[Point3; 7]; 2],
    /// Shared corner with the previous fillet (trim intersection).
    pub corner_prev: Point3,
    /// Shared corner with the next fillet.
    pub corner_next: Point3,
}

/// Build the paper cap's two-piece border from an apex blend
/// boundary.
pub fn kp20_paper_border_from_blend_boundary(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
) -> Result<(Kp20PaperGuideBorder, Kp20PaperBorderExtractReport), String> {
    let segment_count = boundary.segments.len();
    if segment_count < 6 || !segment_count.is_multiple_of(2) {
        return Err(format!(
            "paper border needs an alternating 2n-segment boundary with n >= 3, got {segment_count} segments"
        ));
    }
    let n = segment_count / 2;

    // Orient every segment along the loop: segment k must start where
    // segment k-1 ends.
    let mut recovered_face_surfaces = 0_usize;
    let mut curves: Vec<OrientedSegment> = Vec::with_capacity(segment_count);
    for segment in &boundary.segments {
        curves.push(OrientedSegment::new(
            segment,
            fillet_surfaces,
            brep,
            boundary,
            &mut recovered_face_surfaces,
        )?);
    }
    let mut max_loop_gap = 0.0_f64;
    for k in 0..segment_count {
        let previous_end = curves[(k + segment_count - 1) % segment_count].end();
        let (gap_fwd, gap_rev) = {
            let this = &curves[k];
            (
                previous_end.distance_to(&this.start()),
                previous_end.distance_to(&this.end()),
            )
        };
        if gap_rev < gap_fwd {
            curves[k].reversed = !curves[k].reversed;
            max_loop_gap = max_loop_gap.max(gap_rev);
        } else {
            max_loop_gap = max_loop_gap.max(gap_fwd);
        }
    }

    // Normalize the phase so segment 0 is a FilletEnd (sector k =
    // FilletEnd[k] then FaceEdge[k], with c at their junction).
    if matches!(boundary.segments[0], BoundarySegment::FaceEdge { .. }) {
        curves.rotate_left(1);
    }

    // Wind counterclockwise about the boundary's mean normal so the
    // sector order matches the cap chart's handedness.
    let centroid = loop_centroid(&curves);
    let mean_normal = loop_mean_normal(&curves, centroid);
    let winding = loop_winding(&curves, centroid, mean_normal);
    let reversed_for_orientation = winding < 0.0;
    if reversed_for_orientation {
        curves.reverse();
        for segment in curves.iter_mut() {
            segment.reversed = !segment.reversed;
        }
        // Keep the FilletEnd-first phase after reversal.
        if !curves[0].is_fillet_end {
            curves.rotate_left(1);
        }
    }

    // Uniform ring-march depth from the mean segment chord.
    let mean_chord = curves
        .iter()
        .map(|segment| segment.start().distance_to(&segment.end()))
        .sum::<f64>()
        / segment_count as f64;
    let depth = 0.8 * mean_chord;

    // Boundary conditioning to the g∘χ̃ shape class via macro-edge
    // fitting. The characteristic layout corners ONLY at the n sector
    // corners c and is C²-smooth at the n spoke points m; a raw
    // fillet boundary corners at all 2n junctions, and a single bi-3
    // span lacks the freedom to be C² at m without disturbing its c
    // corner (measured: second-difference unification flattens c
    // from 119° to 28°). Each macro-edge (c_k → c_{k+1}, through m)
    // is therefore refit as a 2-span cubic spline with:
    // - both c endpoints and the m junction interpolated exactly,
    // - the TRUE segment tangent directions preserved at both c ends
    //   (the real ~119° corner, matching χ̃'s corner class),
    // - C¹/C² at m by construction (the constraints eliminate all
    //   but two scalar unknowns — the c-end tangent magnitudes),
    // - the remaining two scalars least-squares fit to dense samples
    //   of the true boundary segments.
    // Junction j sits at the start of segment j; with the
    // FilletEnd-first phase, ODD junctions are the sector corners c.
    let mut rims: Vec<[Point3; 4]> = curves.iter().map(|segment| segment.rim_row()).collect();
    for k in 0..n {
        let j = 2 * k + 1; // c junction between segment j-1 (arc) and j (bridge)
        let bridge_index = j % segment_count;
        let arc_index = (j + 1) % segment_count;
        let c_start = rims[bridge_index][0];
        let m_point = rims[arc_index][0];
        let c_end = rims[arc_index][3];
        // True corner directions from the raw interpolated rims.
        let t0 = (rims[bridge_index][1] - rims[bridge_index][0])
            .normalized()
            .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
        let t1 = (rims[arc_index][3] - rims[arc_index][2])
            .normalized()
            .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
        // Solve min over (alpha, beta) of the sample distance, with
        // P1 = c_start + alpha t0, Q2 = c_end - beta t1,
        // P2 = (c_start + 4m - c_end + alpha t0 + beta t1)/4,
        // Q1 = 2m - P2 (the C1/C2-at-m eliminations).
        const SAMPLES_PER_SEGMENT: usize = 9;
        let mut ata = [[0.0_f64; 2]; 2];
        let mut atb = [0.0_f64; 2];
        for (span, segment_index) in [(0_usize, bridge_index), (1, arc_index)] {
            for sample in 0..SAMPLES_PER_SEGMENT {
                let t = (sample + 1) as f64 / (SAMPLES_PER_SEGMENT + 1) as f64;
                let target = curves[segment_index].point(t);
                // Cubic Bernstein weights on this span.
                let sfp = 1.0 - t;
                let w = [
                    sfp * sfp * sfp,
                    3.0 * sfp * sfp * t,
                    3.0 * sfp * t * t,
                    t * t * t,
                ];
                // Express the span point as constant + alpha·a + beta·b.
                // Span 0 (P): [c_start, P1, P2, m]; span 1 (Q): [m, Q1, Q2, c_end].
                let combine = |points: [(Point3, f64); 3]| -> [f64; 3] {
                    let mut acc = [0.0_f64; 3];
                    for (point, weight) in points {
                        acc[0] += point.x * weight;
                        acc[1] += point.y * weight;
                        acc[2] += point.z * weight;
                    }
                    acc
                };
                let p2_const = Point3::new(
                    (c_start.x + 4.0 * m_point.x - c_end.x) / 4.0,
                    (c_start.y + 4.0 * m_point.y - c_end.y) / 4.0,
                    (c_start.z + 4.0 * m_point.z - c_end.z) / 4.0,
                );
                let (constant, a_vec, b_vec) = if span == 0 {
                    let constant =
                        combine([(c_start, w[0] + w[1]), (p2_const, w[2]), (m_point, w[3])]);
                    (constant, t0 * (w[1] + w[2] / 4.0), t1 * (w[2] / 4.0))
                } else {
                    let q1_const = Point3::new(
                        2.0 * m_point.x - p2_const.x,
                        2.0 * m_point.y - p2_const.y,
                        2.0 * m_point.z - p2_const.z,
                    );
                    let constant =
                        combine([(m_point, w[0]), (q1_const, w[1]), (c_end, w[2] + w[3])]);
                    (constant, t0 * (-w[1] / 4.0), t1 * (-w[1] / 4.0 - w[2]))
                };
                let r = Vec3::new(
                    target.x - constant[0],
                    target.y - constant[1],
                    target.z - constant[2],
                );
                ata[0][0] += a_vec.dot(&a_vec);
                ata[0][1] += a_vec.dot(&b_vec);
                ata[1][0] += a_vec.dot(&b_vec);
                ata[1][1] += b_vec.dot(&b_vec);
                atb[0] += a_vec.dot(&r);
                atb[1] += b_vec.dot(&r);
            }
        }
        let det = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0];
        if det.abs() < 1.0e-18 {
            continue;
        }
        let alpha = (atb[0] * ata[1][1] - atb[1] * ata[0][1]) / det;
        let beta = (ata[0][0] * atb[1] - ata[1][0] * atb[0]) / det;
        let p1 = c_start + t0 * alpha;
        let q2 = c_end - t1 * beta;
        let p2 = Point3::new(
            (c_start.x + 4.0 * m_point.x - c_end.x + alpha * t0.x + beta * t1.x) / 4.0,
            (c_start.y + 4.0 * m_point.y - c_end.y + alpha * t0.y + beta * t1.y) / 4.0,
            (c_start.z + 4.0 * m_point.z - c_end.z + alpha * t0.z + beta * t1.z) / 4.0,
        );
        let q1 = Point3::new(
            2.0 * m_point.x - p2.x,
            2.0 * m_point.y - p2.y,
            2.0 * m_point.z - p2.z,
        );
        rims[bridge_index] = [c_start, p1, p2, m_point];
        rims[arc_index] = [m_point, q1, q2, c_end];
    }

    // Pass 2: cross rows calibrated to the characteristic border's
    // own parameterization structure. For each of the four piece-end
    // sites (piece A/B × start/end), the χ̃ border for this valence
    // supplies the cross-leg's length ratio to the rim leg and its
    // SIGNED in-plane angle from the rim tangent (measured: ~90° at
    // m with ratio ~1.08; ~94°/86° mirror-skewed at c with ratio
    // 1.0 at n = 5). The real border's cross rows reproduce those
    // ratios and angles in each surface's tangent plane — the
    // ρ-reparameterization freedom applied so the input reaches the
    // assembly in exactly the parameterization the τ/b machinery was
    // built around. Marches are projected onto the owning surfaces.
    let calibration = chi_march_calibration(n)
        .ok_or_else(|| "characteristic march calibration failed".to_string())?;
    let mut first_order_marches = 0_usize;
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let rim_a = elevate_cubic_row_to_sextic(rims[2 * k]);
        let rim_b = elevate_cubic_row_to_sextic(rims[2 * k + 1]);
        // Corner-compatible cross jets at c (shared by the two pieces),
        // then one shared twist step next to c (see `section`).
        let corner_a = rim_corner_jet(&rim_b, true);
        let corner_b = rim_corner_jet(&rim_a, false);
        let samples_a = curves[2 * k].border_piece_calibrated(
            &rim_a,
            &calibration[0],
            boundary.vertex_position,
            &mut first_order_marches,
            [None, Some(corner_a)],
        );
        let samples_b = curves[2 * k + 1].border_piece_calibrated(
            &rim_b,
            &calibration[1],
            boundary.vertex_position,
            &mut first_order_marches,
            [Some(corner_b), None],
        );
        let (pins_a, pins_a2) = corner_pins(corner_a, true, None);
        let (pins_b, pins_b2) = corner_pins(corner_b, false, None);
        let natural_a = section::fit_sextic_field_pinned(&samples_a.jets, pins_a);
        let natural_b = section::fit_sextic_field_pinned(&samples_b.jets, pins_b);
        let station_normal = |segment: &OrientedSegment, row0: &[Point3; 7], t: f64| -> Vec3 {
            segment
                .surface_normal(popkern_core::geometry::gregory::eval_sextic_row(row0, t))
                .and_then(|n| n.normalized())
                .unwrap_or(Vec3::new(0.0, 0.0, 0.0))
        };
        let twist = section::corner_twist_step(
            natural_a[5] - natural_a[6],
            natural_b[1] - natural_b[0],
            corner_leg(&rim_a, true),
            corner_leg(&rim_b, false),
            station_normal(&curves[2 * k], &rim_a, 11.0 / 12.0),
            station_normal(&curves[2 * k + 1], &rim_b, 1.0 / 12.0),
        );
        let (pins_a, _) = corner_pins(corner_a, true, Some(twist));
        let (pins_b, _) = corner_pins(corner_b, false, Some(twist));
        let piece_a = rows_from_cross_fields(
            &rim_a,
            &section::fit_sextic_field_pinned(&samples_a.jets, pins_a),
            &section::fit_sextic_field_pinned(&samples_a.curvatures, pins_a2),
        );
        let piece_b = rows_from_cross_fields(
            &rim_b,
            &section::fit_sextic_field_pinned(&samples_b.jets, pins_b),
            &section::fit_sextic_field_pinned(&samples_b.curvatures, pins_b2),
        );
        // Sector rows run [m^k → c] then [c → m^{k+1}]: the FilletEnd
        // piece starts the sector at m and ends at c.
        sectors.push([
            Kp20PaperGuideBorderPiece { rows: piece_a },
            Kp20PaperGuideBorderPiece { rows: piece_b },
        ]);
    }

    Ok((
        Kp20PaperGuideBorder { sectors },
        Kp20PaperBorderExtractReport {
            n,
            depth,
            max_loop_gap,
            first_order_marches,
            reversed_for_orientation,
            trim_corners: Vec::new(),
            setback_arcs: Vec::new(),
            recovered_face_surfaces,
            section_fit_residual: 0.0,
            section_fit_fallbacks: 0,
            gapped_corners: 0,
            reflex_corners: 0,
            corner_block_mismatch: Vec::new(),
        },
    ))
}

/// Cross-leg calibration at one piece end: length ratio to the rim
/// leg and signed in-plane angle from the rim tangent.
#[derive(Debug, Clone, Copy)]
struct ChiEndCalibration {
    ratio: f64,
    angle: f64,
}

/// Per-piece calibration (start and end).
#[derive(Debug, Clone, Copy)]
struct ChiPieceCalibration {
    start: ChiEndCalibration,
    end: ChiEndCalibration,
}

/// Measure the characteristic border's cross-leg structure for this
/// valence: `[piece A ([m→c]-like near_c), piece B ([c→m]-like
/// near_m)]`, each with start/end calibrations.
fn chi_march_calibration(n: usize) -> Option<[ChiPieceCalibration; 2]> {
    use popkern_core::geometry::gregory::cc_lifted_ring_border_pieces;
    let pieces = cc_lifted_ring_border_pieces(n, |_, _| 0.0).ok()?;
    let [near_c, near_m] = &pieces[0];
    let plane_normal = Vec3::new(0.0, 0.0, 1.0);
    let measure = |rows: &[[Point3; 4]; 3], col: usize| -> Option<ChiEndCalibration> {
        let rim = if col == 0 {
            rows[0][1] - rows[0][0]
        } else {
            rows[0][3] - rows[0][2]
        };
        let cross = rows[1][col] - rows[0][col];
        let rim_len = rim.length();
        let cross_len = cross.length();
        if rim_len < 1.0e-14 || cross_len < 1.0e-14 {
            return None;
        }
        let t_hat = rim * (1.0 / rim_len);
        let side = plane_normal.cross(&t_hat);
        let angle = cross.dot(&side).atan2(cross.dot(&t_hat));
        Some(ChiEndCalibration {
            ratio: cross_len / rim_len,
            angle,
        })
    };
    Some([
        ChiPieceCalibration {
            start: measure(near_c, 0)?,
            end: measure(near_c, 3)?,
        },
        ChiPieceCalibration {
            start: measure(near_m, 0)?,
            end: measure(near_m, 3)?,
        },
    ])
}

/// Cap-side cross jet at a piece end, in cubic-depth Bernstein form:
/// `d1 = row1 − row0`, `d2 = row2 − 2·row1 + row0`.
#[derive(Debug, Clone, Copy)]
pub(super) struct CornerJet {
    pub d1: Vec3,
    pub d2: Vec3,
}

/// The corner-compatible cross jet at a sector corner `c`, read from
/// the NEIGHBOURING rim row (PK-398). At `c` the corner piece's two
/// parameter directions are the two rims, so piece A's cap-side
/// depth row at `c` IS piece B's rim row in cubic-equivalent Bernstein
/// form (`leg ×2`, second difference ×5 relative to the degree-6
/// row). The KP20 corner overlap averages A's rim coefficients next
/// to `c` with B's cross coefficients at `c`; a χ̃-calibrated guess
/// there dragged the served rim off the exact section (measured
/// along-seam curvature 6.83 vs 6.70 at 30 % of the quarter arc from
/// `c`, g2 1.4e-1), while this jet makes the overlap a no-op on the
/// rim by construction. `neighbour_starts_at_c` selects which end of
/// the neighbour's row touches `c`.
pub(super) fn rim_corner_jet(
    neighbour_row0: &[Point3; 7],
    neighbour_starts_at_c: bool,
) -> CornerJet {
    let r = neighbour_row0;
    let (a, b, c) = if neighbour_starts_at_c {
        (r[0], r[1], r[2])
    } else {
        (r[6], r[5], r[4])
    };
    CornerJet {
        d1: b - a,
        d2: ((c - b) - (b - a)) * 1.25,
    }
}

/// Sampled cap-side cross fields along one rim row: `(t, value)`.
pub(super) struct CrossFieldSamples {
    /// Cross-jet samples `row1 − row0`.
    pub jets: Vec<(f64, Vec3)>,
    /// Cross-curvature samples `row2 − 2·row1 + row0`.
    pub curvatures: Vec<(f64, Vec3)>,
}

/// Assemble the three depth rows from the rim row and the fitted
/// cross-jet / cross-curvature rows.
pub(super) fn rows_from_cross_fields(
    row0: &[Point3; 7],
    d1: &[Vec3; 7],
    d2: &[Vec3; 7],
) -> [[Point3; 7]; 3] {
    let mut rows = [*row0, *row0, *row0];
    for i in 0..7 {
        rows[1][i] = row0[i] + d1[i];
        rows[2][i] = rows[1][i] + (rows[1][i] - row0[i]) + d2[i];
    }
    rows
}

/// Corner pins for a piece ending at `c` (`at_end`) or starting there:
/// the jet value and, when given, the shared twist step next to `c`.
pub(super) fn corner_pins(
    corner: CornerJet,
    at_end: bool,
    twist: Option<Vec3>,
) -> ([Option<Vec3>; 7], [Option<Vec3>; 7]) {
    let mut d1 = [None; 7];
    let mut d2 = [None; 7];
    if at_end {
        d1[6] = Some(corner.d1);
        d2[6] = Some(corner.d2);
        if let Some(step) = twist {
            d1[5] = Some(corner.d1 + step);
        }
    } else {
        d1[0] = Some(corner.d1);
        d2[0] = Some(corner.d2);
        if let Some(step) = twist {
            d1[1] = Some(corner.d1 + step);
        }
    }
    (d1, d2)
}

/// The c-quarter rim leg from `c` (the corner block's rim leg).
pub(super) fn corner_leg(row0: &[Point3; 7], at_end: bool) -> Vec3 {
    if at_end {
        (row0[5] - row0[6]) * 0.5
    } else {
        (row0[1] - row0[0]) * 0.5
    }
}

/// One boundary segment with loop orientation and its owning surface.
struct OrientedSegment {
    curve: popkern_core::geometry::NurbsCurve,
    surface: SurfaceGeometry,
    is_fillet_end: bool,
    reversed: bool,
}

impl OrientedSegment {
    fn new(
        segment: &BoundarySegment,
        fillet_surfaces: &[FilletSurfaceResult],
        brep: &BRep,
        boundary: &BlendBoundary,
        recovered_face_surfaces: &mut usize,
    ) -> Result<Self, String> {
        match segment {
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => {
                let surface = fillet_surfaces
                    .get(*fillet_index)
                    .ok_or_else(|| format!("fillet index {fillet_index} out of range"))?
                    .surface
                    .clone();
                Ok(Self {
                    curve: curve.clone(),
                    surface,
                    is_fillet_end: true,
                    reversed: false,
                })
            }
            BoundarySegment::FaceEdge { face, curve } => {
                let surface = match brep.faces.get(*face).and_then(|f| f.surface.as_ref()) {
                    Some(on_face) => on_face.geometry.clone(),
                    // The n-gon boundary is computed before the
                    // corner-face topology surgery, so its face
                    // handles may be stale at dispatch time. Recover
                    // the owning surface by curve incidence, counted
                    // in the report.
                    None => {
                        *recovered_face_surfaces += 1;
                        recover_face_edge_surface(brep, curve, boundary)?
                    }
                };
                Ok(Self {
                    curve: curve.clone(),
                    surface,
                    is_fillet_end: false,
                    reversed: false,
                })
            }
        }
    }

    fn parameter(&self, t: f64) -> f64 {
        let (t0, t1) = self.curve.parameter_range();
        let s = if self.reversed { 1.0 - t } else { t };
        t0 + s * (t1 - t0)
    }

    fn point(&self, t: f64) -> Point3 {
        self.curve.point_at(self.parameter(t))
    }

    fn start(&self) -> Point3 {
        self.point(0.0)
    }

    fn end(&self) -> Point3 {
        self.point(1.0)
    }

    /// Closest surface point to `p` on the owning surface.
    fn surface_foot(&self, p: Point3) -> Option<Point3> {
        use popkern_core::geometry::surface::SurfaceEval;
        use popkern_core::geometry::uv_projection::project_point_to_uv;
        let uv = project_point_to_uv(&self.surface, &p)?;
        Some(self.surface.point_at(uv.0, uv.1))
    }

    /// Surface normal near `p` on the owning surface (used to build
    /// the shared tangent plane at smooth spoke junctions).
    fn surface_normal(&self, p: Point3) -> Option<Vec3> {
        use popkern_core::geometry::surface::SurfaceEval;
        use popkern_core::geometry::uv_projection::project_point_to_uv;
        let uv = project_point_to_uv(&self.surface, &p)?;
        Some(self.surface.normal_at(uv.0, uv.1))
    }

    /// Rim row: cubic through four curve samples at t = 0, 1/3, 2/3,
    /// 1 (endpoint-tangent Hermite is unreliable on the boundary
    /// NURBS, whose end derivatives are not chord-scaled).
    fn rim_row(&self) -> [Point3; 4] {
        let p0 = self.point(0.0);
        let q1 = self.point(1.0 / 3.0);
        let q2 = self.point(2.0 / 3.0);
        let p3 = self.point(1.0);
        let interior = |a: f64, b: f64, c: f64, d: f64| {
            Point3::new(
                (a * p0.x + b * q1.x + c * q2.x + d * p3.x) / 6.0,
                (a * p0.y + b * q1.y + c * q2.y + d * p3.y) / 6.0,
                (a * p0.z + b * q1.z + c * q2.z + d * p3.z) / 6.0,
            )
        };
        [
            p0,
            interior(-5.0, 18.0, -9.0, 2.0),
            interior(2.0, -9.0, 18.0, -5.0),
            p3,
        ]
    }

    /// χ̃-calibrated march vector at one end of the rim row (`col` 0 or
    /// 6): the rim tangent rotated by the characteristic signed angle
    /// in this surface's tangent plane and scaled by the
    /// characteristic cross/rim ratio (measured on cubic χ̃ rows, so
    /// the rim leg is the CUBIC-EQUIVALENT `2·(b₁ − b₀)` of the
    /// degree-6 row).
    fn calibrated_march(
        &self,
        row0: &[Point3; 7],
        col: usize,
        cal: &ChiEndCalibration,
        vertex: Point3,
    ) -> Vec3 {
        let rim = if col == 0 {
            (row0[1] - row0[0]) * 2.0
        } else {
            (row0[6] - row0[5]) * 2.0
        };
        let rim_len = rim.length().max(1.0e-14);
        let t_hat = rim * (1.0 / rim_len);
        let p = row0[col];
        let normal = self.surface_normal(p).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        // Orient the frame so `normal × t_hat` points away from the
        // blend vertex (into the ring), matching the χ̃ frame's
        // outward side.
        let mut side = normal.cross(&t_hat);
        if side.dot(&(p - vertex)) < 0.0 {
            side = -side;
        }
        let side = side.normalized().unwrap_or(side);
        (t_hat * cal.angle.cos() + side * cal.angle.sin()) * (cal.ratio * rim_len)
    }

    /// Sample the cap-side cross-jet and cross-curvature fields along
    /// the rim row (PK-398). The march vector blends `march_ends`
    /// linearly, kept in the tangent plane at each true curve point
    /// `S(t)`. Marching from the row's control points — the earlier
    /// construction — started every interior march off the surface by
    /// the control polygon's sagitta (1.4e-3 at r = 0.1), which the
    /// second difference read as cross-curvature (−1.4e-3/ℓ² ≈ −2.0
    /// against a straight sweep's 0) and served as the cap↔strip
    /// transverse G2 class (1.4 normalized on every rim+verticals seam).
    fn march_cross_field(
        &self,
        row0: &[Point3; 7],
        march_ends: [Vec3; 2],
        first_order_marches: &mut usize,
    ) -> CrossFieldSamples {
        const FIELD_SAMPLES: usize = 25;
        let mut jets = Vec::with_capacity(FIELD_SAMPLES);
        let mut curvatures = Vec::with_capacity(FIELD_SAMPLES);
        for j in 0..FIELD_SAMPLES {
            let t = j as f64 / (FIELD_SAMPLES - 1) as f64;
            let rim_point = popkern_core::geometry::gregory::eval_sextic_row(row0, t);
            // March from the rim point's SURFACE foot: at a gapped
            // corner (trim curves that miss each other — 3.4e-3 on the
            // blunt pentagon G2 sweeps) the rim row interpolates a
            // corner 1.6e-3 off the strip, and a march started there
            // reads the offset as cross-curvature (offset/ℓ² ≈ 16)
            // exactly like the control-point march did. The fields are
            // relative to the row, so the rim reconcile's later landing
            // of the row carries them along.
            let q0 = self.surface_foot(rim_point).unwrap_or(rim_point);
            let n_hat = self
                .surface_normal(q0)
                .and_then(|n| n.normalized())
                .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
            // Smootherstep blend: zero slope AND curvature at both
            // ends, so the corner block's third-order slots agree
            // between the two pieces by construction. The χ̃ scale
            // convention makes the cross leg ~2× longer at the spoke
            // point m than the corner-compatible leg at c, so a linear
            // blend leaves a length KINK at m (both halves rising
            // toward it) and a twist step at c; the field must be C²
            // across m and twist-free at c for the shared junction
            // pins to agree with the samples (measured: linear blend
            // forced the constrained fit 1.1e-2 off its samples).
            let w = t * t * t * (10.0 + t * (-15.0 + 6.0 * t));
            let blend = march_ends[0] * (1.0 - w) + march_ends[1] * w;
            let in_plane = blend - n_hat * blend.dot(&n_hat);
            let v = in_plane
                .normalized()
                .map(|u| u * blend.length())
                .unwrap_or(blend);
            let (q1, q2) = self.ring_march_vector(q0, v, first_order_marches);
            // Cap-side C² prolongation. The WP1b point reflection
            // r1 = 2q0 − q1 mirrors a SURFACE CHORD through q0: the
            // chord leaves the strip's tangent plane by ~κℓ²/2, so
            // the cap's cross-jet tilted out of the surface by
            // ~κℓ/2 radians — measured as a UNIFORM 2.6°–6.9° rim
            // tangency error along every setback arc (PK-212). G1
            // across the rim requires the jet IN the tangent plane:
            // project the chord's normal component out before
            // reflecting.
            let chord = q1 - q0;
            let jet = -(chord - n_hat * chord.dot(&n_hat));
            // Curvature semantics of a Bernstein depth row: derivatives
            // at the rim are 3·(r1 − r0) and 6·(r2 − 2r1 + r0). With
            // the leg equal to the sample spacing ℓ the row represents
            // the surface at depth parameter α = 3ℓ, so its second
            // difference is α²f''/6 = 1.5·ℓ²f'' — 1.5× the SAMPLE
            // second difference (≈ ℓ²f''). Exact for straight marches
            // either way; on curved-spine strips the factor carries the
            // strip's full curvature into the cap.
            let second = (q2 - q1) - (q1 - q0);
            jets.push((t, jet));
            curvatures.push((t, second * 1.5));
        }
        CrossFieldSamples { jets, curvatures }
    }

    /// Cross-field samples with χ̃-calibrated ends, or corner-compatible
    /// jets at `c` where `corner_jets` are given (the cap-side jet is
    /// the negated march vector).
    fn border_piece_calibrated(
        &self,
        row0: &[Point3; 7],
        calibration: &ChiPieceCalibration,
        vertex: Point3,
        first_order_marches: &mut usize,
        corner_jets: [Option<CornerJet>; 2],
    ) -> CrossFieldSamples {
        let start = match corner_jets[0] {
            Some(corner) => -corner.d1,
            None => self.calibrated_march(row0, 0, &calibration.start, vertex),
        };
        let end = match corner_jets[1] {
            Some(corner) => -corner.d1,
            None => self.calibrated_march(row0, 6, &calibration.end, vertex),
        };
        self.march_cross_field(row0, [start, end], first_order_marches)
    }

    /// March into the ring by the leg vector `v` on the owning
    /// surface: returns surface points near `q0 + v` and `q0 + 2v`
    /// (exact in-plane; radially projected onto cylinders;
    /// closest-point-projected onto NURBS; tangent-plane first-order
    /// otherwise).
    fn ring_march_vector(
        &self,
        q0: Point3,
        v: Vec3,
        first_order_marches: &mut usize,
    ) -> (Point3, Point3) {
        let raw1 = q0 + v;
        let raw2 = q0 + v * 2.0;
        match &self.surface {
            SurfaceGeometry::Plane { origin, normal, .. } => {
                let snap = |p: Point3| {
                    let d = (p - *origin).dot(normal);
                    p - *normal * d
                };
                (snap(raw1), snap(raw2))
            }
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
                ..
            } => {
                let snap = |p: Point3| {
                    let along = (p - *origin).dot(axis);
                    let foot = *origin + *axis * along;
                    let radial = p - foot;
                    let unit = radial.normalized().unwrap_or(radial);
                    foot + unit * *radius
                };
                (snap(raw1), snap(raw2))
            }
            SurfaceGeometry::Nurbs(nurbs) => {
                // Snap the march onto the NURBS sweep via closest-point
                // projection (second-order accurate ring sampling).
                use popkern_core::geometry::surface::SurfaceEval;
                use popkern_core::geometry::uv_projection::project_point_to_uv;
                let snap = |p: Point3| -> Option<Point3> {
                    let (u, v) = project_point_to_uv(&self.surface, &p)?;
                    if popkern_core::trace::trace_enabled("kp-march") {
                        let out = nurbs.point_at(u, v);
                        eprintln!(
                            "      nurbs snap: raw=({:.5},{:.5},{:.5}) uv=({:.4},{:.4}) out=({:.5},{:.5},{:.5})",
                            p.x, p.y, p.z, u, v, out.x, out.y, out.z,
                        );
                    }
                    Some(nurbs.point_at(u, v))
                };
                match (snap(raw1), snap(raw2)) {
                    (Some(q1), Some(q2)) => (q1, q2),
                    _ => {
                        *first_order_marches += 1;
                        (raw1, raw2)
                    }
                }
            }
            _ => {
                *first_order_marches += 1;
                (raw1, raw2)
            }
        }
    }
}

fn loop_centroid(curves: &[OrientedSegment]) -> Point3 {
    let mut acc = Vec3::new(0.0, 0.0, 0.0);
    let mut count = 0.0;
    for segment in curves {
        let p = segment.start();
        acc = acc + Vec3::new(p.x, p.y, p.z);
        count += 1.0;
    }
    Point3::new(acc.x / count, acc.y / count, acc.z / count)
}

fn loop_mean_normal(curves: &[OrientedSegment], centroid: Point3) -> Vec3 {
    let mut acc = Vec3::new(0.0, 0.0, 0.0);
    for segment in curves {
        let a = segment.start() - centroid;
        let b = segment.end() - centroid;
        acc = acc + a.cross(&b);
    }
    acc.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0))
}

fn loop_winding(curves: &[OrientedSegment], centroid: Point3, normal: Vec3) -> f64 {
    let mut acc = 0.0;
    for segment in curves {
        let a = segment.start() - centroid;
        let b = segment.end() - centroid;
        acc += a.cross(&b).dot(&normal);
    }
    acc
}

/// De Casteljau halves of one cubic row.
fn split_cubic_row(row: &[Point3; 4]) -> ([Point3; 4], [Point3; 4]) {
    let mid =
        |a: Point3, b: Point3| Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
    let p01 = mid(row[0], row[1]);
    let p12 = mid(row[1], row[2]);
    let p23 = mid(row[2], row[3]);
    let p012 = mid(p01, p12);
    let p123 = mid(p12, p23);
    let center = mid(p012, p123);
    ([row[0], p01, p012, center], [center, p123, p23, row[3]])
}

#[cfg(test)]
#[path = "kp20_paper_border_extract_tests.rs"]
mod tests;
