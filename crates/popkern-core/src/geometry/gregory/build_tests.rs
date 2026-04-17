//! Tests for the `RectangularGregory::from_boundary_beziers*` and
//! `TriangularGregory::from_boundary_beziers` constructors. Kept in a
//! separate file from the evaluator tests so neither grows past the
//! 500-line structural cap.

use super::*;
use crate::geometry::Point3;

/// §5.3 C1 — a rectangular patch built via `from_boundary_beziers_no_twist`
/// on a flat bicubic evaluates identically to directly-constructed CPs
/// and produces a flat surface. Confirms the builder's corner-indexing
/// and grid-assembly logic.
#[test]
fn rect_from_boundary_beziers_no_twist_produces_flat_patch() {
    let top_edge = [
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0 / 3.0, 1.0, 0.0),
        Point3::new(2.0 / 3.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ];
    let left_edge = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 1.0 / 3.0, 0.0),
        Point3::new(0.0, 2.0 / 3.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let right_edge_reversed = [
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(1.0, 2.0 / 3.0, 0.0),
        Point3::new(1.0, 1.0 / 3.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    ];
    let bottom_edge_reversed = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0 / 3.0, 0.0, 0.0),
        Point3::new(1.0 / 3.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
    ];
    let boundaries = [
        left_edge,
        top_edge,
        right_edge_reversed,
        bottom_edge_reversed,
    ];
    let patch = RectangularGregory::from_boundary_beziers_no_twist(boundaries, 1e-9)
        .expect("flat bicubic should build");
    for u in [0.0, 0.3, 0.5, 0.7, 1.0] {
        for v in [0.0, 0.3, 0.5, 0.7, 1.0] {
            let p = patch.point_at(u, v);
            assert!(p.z.abs() < 1e-10, "z should be 0, got {} at ({u},{v})", p.z);
            assert!(
                (p.x - u).abs() < 1e-10 && (p.y - v).abs() < 1e-10,
                "bilinear flat patch: point_at({u},{v}) = ({}, {}, {})",
                p.x,
                p.y,
                p.z,
            );
        }
    }
}

/// §5.3 C1 — corner mismatch detection: if the two boundaries through
/// a corner disagree on that corner's position by more than `tol`,
/// the builder rejects with a `CornerMismatch`.
#[test]
fn rect_from_boundary_beziers_rejects_inconsistent_corners() {
    let mut boundaries = [
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0 / 3.0, 0.0),
            Point3::new(0.0, 2.0 / 3.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        [
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0 / 3.0, 1.0, 0.0),
            Point3::new(2.0 / 3.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        [
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 2.0 / 3.0, 0.0),
            Point3::new(1.0, 1.0 / 3.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        [
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0 / 3.0, 0.0, 0.0),
            Point3::new(1.0 / 3.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
        ],
    ];
    boundaries[1][3] = Point3::new(1.0, 1.0, 0.1);
    let inner_cps = [(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)); 4];
    let err = RectangularGregory::from_boundary_beziers(boundaries, inner_cps, 1e-9).unwrap_err();
    match err {
        GregoryBuildError::CornerMismatch { corner, distance } => {
            assert_eq!(
                corner, 2,
                "mismatch between boundaries 1 and 2 should be corner 2"
            );
            assert!(distance > 0.05, "distance = {distance}");
        }
    }
}

/// §5.3 C1 — a flat triangular patch built via `from_boundary_beziers`
/// reproduces the flat-simplex result. Verifies boundary CP ordering
/// in the Longhi storage scheme.
#[test]
fn tri_from_boundary_beziers_produces_flat_simplex() {
    let corner0 = Point3::new(1.0, 0.0, 0.0);
    let corner1 = Point3::new(0.0, 1.0, 0.0);
    let corner2 = Point3::new(0.0, 0.0, 0.0);
    let edge_0 = [
        corner1,
        Point3::new(0.0, 2.0 / 3.0, 0.0),
        Point3::new(0.0, 1.0 / 3.0, 0.0),
        corner2,
    ];
    let edge_1 = [
        corner2,
        Point3::new(1.0 / 3.0, 0.0, 0.0),
        Point3::new(2.0 / 3.0, 0.0, 0.0),
        corner0,
    ];
    let edge_2 = [
        corner0,
        Point3::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
        Point3::new(1.0 / 3.0, 2.0 / 3.0, 0.0),
        corner1,
    ];
    let interior = Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.0);
    let patch = TriangularGregory::from_boundary_beziers(
        [edge_0, edge_1, edge_2],
        [interior, interior, interior],
        1e-9,
    )
    .expect("flat tri should build");
    for (u, v) in [(0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.25, 0.25), (0.4, 0.4)] {
        let p = patch.point_at_uv(u, v);
        assert!(
            p.z.abs() < 1e-10 && (p.x - u).abs() < 1e-10 && (p.y - v).abs() < 1e-10,
            "flat triangular patch point_at_uv({u},{v}) = ({}, {}, {})",
            p.x,
            p.y,
            p.z,
        );
    }
}
