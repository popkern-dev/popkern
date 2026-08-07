use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::gregory::GregoryPatch;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

#[test]
fn heal_holes_closes_the_open_box_with_the_single_patch_g1_cap() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let report = heal_holes(&mut brep, Tolerance::default().linear());

    assert_eq!(report.holes_detected, 1);
    assert_eq!(report.holes_committed, 1);
    assert_eq!(report.attempts.len(), 1);
    let attempt = &report.attempts[0];
    // PK-224: the 4-sided hole serves the single-patch G1 cap — no
    // internal seams, ribbon-tangent to every wall. The polar-fan/
    // PnS/Coons/ruled rungs are the recorded fallbacks and never
    // fire here.
    assert_eq!(attempt.committed_by, Some(HoleFillerKind::SinglePatchG1));
    assert_eq!(attempt.num_sides, 4);
    assert!(attempt.single_patch_decline.is_none());
    assert!(attempt.polar_fan_decline.is_none());
    assert!(attempt.pns_decline.is_none());
    assert!(attempt.coons_decline.is_none());
    assert!(attempt.ruled_decline.is_none());

    assert!(brep.boundary_edges.is_empty());
    brep.validate_euler_poincare().expect("strict Euler holds");
    brep.validate().expect("full fast validation holds");
}

#[test]
fn heal_holes_on_a_closed_body_is_a_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("unit cube");
    let before = serde_json::to_string(&brep).expect("serializes");
    let report = heal_holes(&mut brep, Tolerance::default().linear());
    assert_eq!(report.holes_detected, 0);
    assert_eq!(report.holes_committed, 0);
    assert!(report.attempts.is_empty());
    assert_eq!(serde_json::to_string(&brep).expect("serializes"), before);
}

/// PK-248 acceptance: N-sided open-prism holes heal with the
/// polar-wedge G1 fan, and the SERVED body measures G1 at the
/// commercial 5e-2 bar across every internal spoke seam and every
/// cap↔wall rim seam (serving-boundary measurement — the PK-220
/// lesson: commit-time numbers are not the serving truth).
fn assert_open_prism_serves_polar_fan_g1(n_sides: usize) {
    let mut brep = BRep::open_prism(n_sides, 1.0, 1.0, Tolerance::default()).expect("open prism");
    let report = heal_holes(&mut brep, Tolerance::default().linear());
    assert_eq!(report.holes_detected, 1);
    assert_eq!(report.holes_committed, 1);
    let attempt = &report.attempts[0];
    assert_eq!(
        attempt.committed_by,
        Some(HoleFillerKind::PolarFanG1),
        "n={n_sides} must serve the polar fan; declines: single={:?} polar={:?}",
        attempt.single_patch_decline,
        attempt.polar_fan_decline,
    );
    assert!(
        attempt.single_patch_decline.is_some(),
        "the 4-sided-only rung must record its side-count decline"
    );
    assert!(brep.boundary_edges.is_empty());
    brep.validate().expect("healed prism validates");

    let (worst_spoke, worst_wall) = worst_served_cap_seam_g1(&brep);
    eprintln!(
        "[pk248] n={n_sides} served worst G1 |sin|: spoke={worst_spoke:.3e} wall={worst_wall:.3e}"
    );
    assert!(
        worst_spoke <= 5.0e-2 && worst_wall <= 5.0e-2,
        "served cap seams must clear the commercial G1 bar, got spoke={worst_spoke:.3e} wall={worst_wall:.3e}"
    );
}

#[test]
fn open_prism3_serves_the_polar_fan_g1_cap() {
    assert_open_prism_serves_polar_fan_g1(3);
}

#[test]
fn open_prism5_serves_the_polar_fan_g1_cap() {
    assert_open_prism_serves_polar_fan_g1(5);
}

#[test]
fn open_prism6_serves_the_polar_fan_g1_cap() {
    assert_open_prism_serves_polar_fan_g1(6);
}

/// Worst mid-band `|sin|` over every edge touching a committed
/// polar-wedge face, as `(worst_spoke_seam, worst_wall_seam)`:
/// spoke seams compare the two wedges' normals, rim seams compare
/// the wedge against its wall. Corner bands are excluded
/// (structural corner argument; fractions 0.25–0.75).
fn worst_served_cap_seam_g1(brep: &BRep) -> (f64, f64) {
    let is_wedge = |face: crate::topology::FaceHandle| -> bool {
        matches!(
            brep.faces.get(face).and_then(|f| f.surface.as_ref()),
            Some(on_face)
                if matches!(
                    on_face.geometry,
                    SurfaceGeometry::Gregory(GregoryPatch::RectangularDeg5(_))
                )
        )
    };
    let face_of = |half_edge: crate::topology::HalfEdgeHandle| {
        brep.half_edges
            .get(half_edge)
            .and_then(|he| brep.loops.get(he.loop_))
            .map(|loop_data| loop_data.face)
    };
    let mut worst_spoke = 0.0_f64;
    let mut worst_wall = 0.0_f64;
    for (_, edge) in &brep.edges {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let (Some(face_a), Some(face_b)) = (face_of(edge.half_edge), face_of(he.twin)) else {
            continue;
        };
        if face_a == face_b || (!is_wedge(face_a) && !is_wedge(face_b)) {
            continue;
        }
        let spoke_seam = is_wedge(face_a) && is_wedge(face_b);
        let surface = |face: crate::topology::FaceHandle| {
            brep.faces
                .get(face)
                .and_then(|f| f.surface.as_ref())
                .map(|on_face| &on_face.geometry)
        };
        let (Some(sa), Some(sb)) = (surface(face_a), surface(face_b)) else {
            continue;
        };
        for &f in &[0.25, 0.375, 0.5, 0.625, 0.75] {
            let t = curve.t_start + f * (curve.t_end - curve.t_start);
            let point = (&curve.geometry as &dyn CurveEval).point_at(t);
            let normal = |s: &SurfaceGeometry| {
                s.closest_point(&point)
                    .and_then(|(u, v)| s.normal_at(u, v).normalized())
            };
            let (Some(na), Some(nb)) = (normal(sa), normal(sb)) else {
                continue;
            };
            let sin = na.cross(&nb).length();
            if spoke_seam {
                worst_spoke = worst_spoke.max(sin);
            } else {
                worst_wall = worst_wall.max(sin);
            }
        }
    }
    (worst_spoke, worst_wall)
}
