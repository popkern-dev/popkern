//! Gates for the first-class polyhedral-net spline patch.
//!
//! The central claim this variant rests on is that it is *not* a
//! geometric improvement — the NURBS image is exact — so the first
//! test proves exactly that. If direct evaluation ever disagreed with
//! the converted form, the variant's whole justification (identity,
//! not accuracy) would be wrong and the doc would need revisiting.

use super::*;
use crate::geometry::nurbs_surface::NurbsSurface;

fn sample_patch() -> PolyhedralSplinePatch {
    // A gently curved net; nothing symmetric, so a transposed or
    // mis-indexed evaluation would show up immediately.
    let mut net = [[Point3::new(0.0, 0.0, 0.0); 5]; 5];
    for (j, row) in net.iter_mut().enumerate() {
        for (i, cp) in row.iter_mut().enumerate() {
            let (x, y) = (i as f64 * 0.25, j as f64 * 0.25);
            *cp = Point3::new(x, y, 0.4 * x * y + 0.1 * x * x - 0.05 * y);
        }
    }
    PolyhedralSplinePatch::new(
        net,
        PolyhedralSplineProvenance {
            valence: 5,
            sector: 2,
        },
    )
}

#[test]
fn direct_evaluation_matches_the_nurbs_image_exactly() {
    let patch = sample_patch();
    let nurbs: NurbsSurface = patch.to_nurbs();
    for iu in 0..=8 {
        for iv in 0..=8 {
            let (u, v) = (iu as f64 / 8.0, iv as f64 / 8.0);
            let direct = patch.point_at(u, v);
            let converted = nurbs.point_at(u, v);
            assert!(
                direct.distance_to(&converted) < 1e-12,
                "at ({u}, {v}) direct {direct:?} != nurbs {converted:?}; \
                 the variant's justification depends on these agreeing"
            );
        }
    }
}

/// Derivatives agree to finite-difference accuracy.
///
/// The bound here is looser than the positional one on purpose: this
/// patch differentiates its Bernstein basis in closed form, while the
/// generic NURBS evaluator differences its own evaluation, so the
/// residual measured (~1e-9) is the NURBS side's truncation error,
/// not a disagreement about the surface. That is a small bonus of
/// keeping the patch in its own form, not a justification for the
/// variant — see the module docs.
#[test]
fn derivatives_match_the_nurbs_image() {
    let patch = sample_patch();
    let nurbs = patch.to_nurbs();
    let fd_bound = 1e-6;
    for iu in 1..8 {
        for iv in 1..8 {
            let (u, v) = (iu as f64 / 8.0, iv as f64 / 8.0);
            let (du_a, du_b) = (patch.partial_u(u, v), nurbs.partial_u(u, v));
            let (dv_a, dv_b) = (patch.partial_v(u, v), nurbs.partial_v(u, v));
            assert!(
                (du_a - du_b).length() < fd_bound,
                "d/du disagrees at ({u}, {v}) beyond finite-difference error: {du_a:?} vs {du_b:?}"
            );
            assert!(
                (dv_a - dv_b).length() < fd_bound,
                "d/dv disagrees at ({u}, {v}) beyond finite-difference error: {dv_a:?} vs {dv_b:?}"
            );
        }
    }
}

#[test]
fn corners_interpolate_the_control_net() {
    let patch = sample_patch();
    let corners = [
        ((0.0, 0.0), patch.net[0][0]),
        ((1.0, 0.0), patch.net[0][4]),
        ((0.0, 1.0), patch.net[4][0]),
        ((1.0, 1.0), patch.net[4][4]),
    ];
    for ((u, v), expected) in corners {
        assert!(
            patch.point_at(u, v).distance_to(&expected) < 1e-12,
            "corner ({u}, {v}) must interpolate its control point"
        );
    }
}

#[test]
fn pole_is_the_zero_corner() {
    let patch = sample_patch();
    assert!(patch.pole().distance_to(&patch.net[0][0]) < 1e-15);
    assert!(patch.pole().distance_to(&patch.point_at(0.0, 0.0)) < 1e-12);
}

#[test]
fn provenance_survives_a_round_trip_through_serde() {
    // Provenance is the entire point of the variant: it has to
    // survive persistence, or a reloaded body loses the identity this
    // was built to keep.
    let patch = sample_patch();
    let json = serde_json::to_string(&patch).expect("patch serialises");
    let back: PolyhedralSplinePatch = serde_json::from_str(&json).expect("patch deserialises");
    assert_eq!(back.provenance.valence, 5);
    assert_eq!(back.provenance.sector, 2);
    assert_eq!(back.net[2][3], patch.net[2][3]);
}

#[test]
fn parameter_range_is_the_unit_square() {
    let patch = sample_patch();
    assert_eq!(patch.parameter_range(), ((0.0, 1.0), (0.0, 1.0)));
}

#[test]
fn normal_is_perpendicular_to_both_partials() {
    let patch = sample_patch();
    let (u, v) = (0.37, 0.61);
    let n = patch.normal_at(u, v).normalized().expect("non-degenerate");
    assert!(n.dot(&patch.partial_u(u, v)).abs() < 1e-9);
    assert!(n.dot(&patch.partial_v(u, v)).abs() < 1e-9);
}
