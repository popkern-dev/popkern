//! Gates for the corner-patched cap: a spherical hexagon with one short
//! side, its neighbour the sphere itself.

use super::super::tests::{spherical_ngon, spherical_polygon};
use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::uv_projection::project_point_to_uv;

/// Side 1 spans 0.15 rad of azimuth between neighbours of 1.0 and 1.05.
const SHORT_SIDE: [f64; 6] = [0.0, 1.0, 1.15, 2.2, 3.5, 4.8];

fn targets() -> VariationalCapTargets {
    VariationalCapTargets::g1(1e-6, 1e-6)
}

/// The blend's corner-cell tangent options (PK-615).
fn options() -> VariationalCapOptions {
    VariationalCapOptions {
        max_refinement_levels: 0,
        energy: super::super::CapEnergy::Bending,
        reference_passes: 3,
        corner_cells: true,
        degree: 9,
        layout_refit_curvature: Some(1e9),
        ..Default::default()
    }
}

#[test]
fn the_short_side_is_the_thin_side() {
    let sides = spherical_polygon(&SHORT_SIDE, 0.6);
    assert_eq!(thin_sides(&sides, 0.6).expect("thin"), vec![1]);
    match thin_sides(&spherical_ngon(6, 0.6), 0.6) {
        Err(VariationalCapError::NoThinSide {
            shortest_ratio,
            threshold,
        }) => assert!((shortest_ratio - 1.0).abs() < 1e-9 && threshold == 0.6),
        other => panic!("expected NoThinSide: {other:?}"),
    }
}

#[test]
fn the_patch_and_the_main_cap_meet_on_the_seam() {
    let sides = spherical_polygon(&SHORT_SIDE, 0.6);
    let cap = fit_corner_patched_cap(
        &sides,
        &targets(),
        &options(),
        &CornerPatchOptions::default(),
        1e-7,
    )
    .expect("corner-patched cap");
    assert_eq!(cap.corners.len(), 1);
    let patch = &cap.corners[0];
    assert_eq!((patch.side, patch.piece.sides.len()), (1, 4));
    assert_eq!(cap.main.sides.len(), 6);
    // The seam closes both loops: it runs from the patch's head end to
    // its tail start, and the main cap carries it reversed.
    let (a, b) = patch.seam.parameter_range();
    let (k0, k1) = (patch.seam.point_at(a), patch.seam.point_at(b));
    let head = &patch.piece.sides[2].curve;
    let tail = &patch.piece.sides[0].curve;
    assert!(k0.distance_to(&head.point_at(head.parameter_range().1)) < 1e-12);
    assert!(k1.distance_to(&tail.point_at(tail.parameter_range().0)) < 1e-12);
    let main_k = &cap.main.sides[1].curve;
    assert!(main_k.point_at(main_k.parameter_range().0).distance_to(&k1) < 1e-12);
    // The ribbon contains K, and its normal along K is the thin side's
    // normal carried across: on the sphere, within a few degrees of the
    // radius at K (the side's own points lie 0.05 rad away).
    let ribbon = &patch.piece.sides[3].neighbour;
    for k in 0..=8 {
        let p = patch.seam.point_at(a + (b - a) * k as f64 / 8.0);
        let (u, v) = project_point_to_uv(ribbon, &p).expect("on the ribbon");
        assert!(ribbon.point_at(u, v).distance_to(&p) < 1e-9);
        let radial = (p - crate::geometry::Point3::new(0.0, 0.0, 0.0))
            .normalized()
            .expect("off centre");
        let normal = ribbon.normal_at(u, v).normalized().expect("regular ribbon");
        let dot = normal.dot(&radial);
        assert!(dot > 0.99, "sample {k}: {dot}");
    }
    // Both pieces meet the ribbon to the G1 target, so each other to
    // twice it; the sphere's caps stay fair.
    for piece in std::iter::once(&cap.main).chain(cap.corners.iter().map(|c| &c.piece)) {
        let g1 = piece
            .cap
            .report
            .rim
            .iter()
            .map(|r| r.g1)
            .fold(0.0, f64::max);
        assert!(g1 <= 1e-6, "rim g1 {g1:.3e}");
        assert!(
            piece.interior_curvature < 3.0,
            "{}",
            piece.interior_curvature
        );
    }
    assert_eq!(
        cap.interior_curvature(),
        cap.main
            .interior_curvature
            .max(cap.corners[0].piece.interior_curvature)
    );
}
