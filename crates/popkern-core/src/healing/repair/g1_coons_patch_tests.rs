use super::*;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};
use crate::healing::repair::coons_patch::CoonsBoundary;

/// Build a rectangular boundary from four line segments.
fn rectangular_boundary(x0: f64, y0: f64, x1: f64, y1: f64) -> CoonsBoundary {
    let bottom = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(x0, y0, 0.0), Point3::new(x1, y0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let top = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(x0, y1, 0.0), Point3::new(x1, y1, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let left = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(x0, y0, 0.0), Point3::new(x0, y1, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let right = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(x1, y0, 0.0), Point3::new(x1, y1, 0.0)],
        weights: vec![1.0, 1.0],
    };
    CoonsBoundary {
        bottom,
        top,
        left,
        right,
    }
}

/// Build uniform normal cross-derivatives (pointing in +Z direction).
fn uniform_normal_cross_derivs(n_samples: usize) -> CrossDerivatives {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    CrossDerivatives {
        bottom: vec![normal; n_samples],
        top: vec![normal; n_samples],
        left: vec![normal; n_samples],
        right: vec![normal; n_samples],
    }
}

#[test]
fn g1_coons_patch_from_rectangle_with_normal_cross_derivs() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    let cross_derivs = uniform_normal_cross_derivs(5);

    let result = construct_g1_coons_patch(&boundary, &cross_derivs, 8, 8, 1e-6);
    assert!(result.is_some(), "G1 construction should succeed");

    let r = result.as_ref();
    // Boundary deviation should remain small for a planar case.
    assert!(
        r.map(|r| r.max_boundary_deviation).unwrap_or(f64::MAX) < 1.0,
        "boundary deviation {} is too large",
        r.map(|r| r.max_boundary_deviation).unwrap_or(f64::MAX),
    );
}

#[test]
fn g1_coons_patch_corners_match_boundary() {
    let boundary = rectangular_boundary(0.0, 0.0, 2.0, 3.0);
    let cross_derivs = uniform_normal_cross_derivs(5);

    let result = construct_g1_coons_patch(&boundary, &cross_derivs, 8, 8, 1e-6);
    assert!(result.is_some());

    let r = result.as_ref();
    let corners = r.map(|r| r.corners);
    assert!(corners.is_some());
    let c = corners.as_ref();
    assert!(
        c.map(|c| c[0].distance_to(&Point3::new(0.0, 0.0, 0.0)))
            .unwrap_or(1.0)
            < 1e-10
    );
    assert!(
        c.map(|c| c[1].distance_to(&Point3::new(2.0, 0.0, 0.0)))
            .unwrap_or(1.0)
            < 1e-10
    );
    assert!(
        c.map(|c| c[2].distance_to(&Point3::new(0.0, 3.0, 0.0)))
            .unwrap_or(1.0)
            < 1e-10
    );
    assert!(
        c.map(|c| c[3].distance_to(&Point3::new(2.0, 3.0, 0.0)))
            .unwrap_or(1.0)
            < 1e-10
    );
}

#[test]
fn g1_coons_patch_rejects_empty_cross_derivs() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    let empty = CrossDerivatives {
        bottom: vec![],
        top: vec![Vec3::new(0.0, 0.0, 1.0)],
        left: vec![Vec3::new(0.0, 0.0, 1.0)],
        right: vec![Vec3::new(0.0, 0.0, 1.0)],
    };
    assert!(construct_g1_coons_patch(&boundary, &empty, 8, 8, 1e-6).is_none());
}

#[test]
fn g1_coons_patch_rejects_small_grid() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    let cross_derivs = uniform_normal_cross_derivs(5);
    assert!(construct_g1_coons_patch(&boundary, &cross_derivs, 3, 8, 1e-6).is_none());
    assert!(construct_g1_coons_patch(&boundary, &cross_derivs, 8, 3, 1e-6).is_none());
}

#[test]
fn g1_coons_patch_non_closing_boundary_returns_none() {
    let bottom = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let top = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let left = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(5.0, 5.0, 5.0), // wrong start
            Point3::new(0.0, 1.0, 0.0),
        ],
        weights: vec![1.0, 1.0],
    };
    let right = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        weights: vec![1.0, 1.0],
    };

    let boundary = CoonsBoundary {
        bottom,
        top,
        left,
        right,
    };
    let cross_derivs = uniform_normal_cross_derivs(5);
    assert!(construct_g1_coons_patch(&boundary, &cross_derivs, 8, 8, 1e-6).is_none());
}

#[test]
fn g1_patch_produces_surface_with_correct_grid_size() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    let cross_derivs = uniform_normal_cross_derivs(5);

    let result = construct_g1_coons_patch(&boundary, &cross_derivs, 10, 10, 1e-6);
    assert!(result.is_some());
    let r = result.as_ref();
    assert!(r.map(|r| r.surface.num_u).unwrap_or(0) == 10);
    assert!(r.map(|r| r.surface.num_v).unwrap_or(0) == 10);
}
