//! §11.1 analytic cylinder / skip-path tests for `build_fillet_surface`.
//!
//! Extracted from `fillet_surface_analytic_tests.rs` to keep that
//! file under the 500-line structural cap. Covers:
//!
//! - parallel-axis cylinder-cylinder → analytic Cylinder fillet
//! - non-parallel cylinder-cylinder → falls through to NURBS sweep
//! - oblique plane-cylinder → falls through to NURBS sweep

use super::*;
use crate::center_curve::CenterCurveResult;
use crate::types::ContinuityLevel;
use popkern_core::geometry::{Point3, Vec3};

/// Stand-in `CenterCurveResult` for the negative (skip) tests —
/// the detector rejects the surface pair before sampling.
fn linear_center_result() -> CenterCurveResult {
    let fp1 = (0..5)
        .map(|i| Point3::new(i as f64, 0.0, 0.0))
        .collect::<Vec<_>>();
    let fp2 = (0..5)
        .map(|i| Point3::new(i as f64, 0.0, 1.0))
        .collect::<Vec<_>>();
    let center_pts: Vec<Point3> = (0..5).map(|i| Point3::new(i as f64, 0.0, 0.5)).collect();
    let center_curve = crate::center_curve::fit_nurbs_through_points(&center_pts).unwrap();
    CenterCurveResult {
        center_curve,
        footprint_1: fp1,
        footprint_2: fp2,
        contact_solutions: None,
    }
}

#[test]
fn fillet_surface_parallel_cylinder_cylinder_produces_analytic_cylinder() {
    let cyl_a = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let cyl_b = SurfaceGeometry::Cylinder {
        origin: Point3::new(3.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let r = 0.3;
    // Rolling-ball center position in the cross-section plane.
    let y_c = ((2.0 + r) * (2.0 + r) - 1.5 * 1.5_f64).sqrt();
    let n = 8usize;
    // Center curve: straight line at (1.5, y_c, z) for z ∈ [0, 2].
    let center_pts: Vec<Point3> = (0..n)
        .map(|i| Point3::new(1.5, y_c, 2.0 * (i as f64) / ((n - 1) as f64)))
        .collect();
    let center_curve = crate::center_curve::fit_nurbs_through_points(&center_pts).unwrap();
    // Footprints: contact points on each cylinder, straight lines
    // parallel to the z axis (doesn't matter for the surface-type
    // assertion, but has to look plausible).
    let fp1: Vec<Point3> = (0..n)
        .map(|i| {
            let z = 2.0 * (i as f64) / ((n - 1) as f64);
            let norm_a = (1.5_f64 * 1.5 + y_c * y_c).sqrt();
            Point3::new(1.5 * 2.0 / norm_a, y_c * 2.0 / norm_a, z)
        })
        .collect();
    let fp2: Vec<Point3> = (0..n)
        .map(|i| {
            let z = 2.0 * (i as f64) / ((n - 1) as f64);
            let dx = 1.5 - 3.0;
            let norm_b = (dx * dx + y_c * y_c).sqrt();
            Point3::new(3.0 + dx * 2.0 / norm_b, y_c * 2.0 / norm_b, z)
        })
        .collect();
    let cr = CenterCurveResult {
        center_curve,
        footprint_1: fp1,
        footprint_2: fp2,
        contact_solutions: None,
    };

    for (s1, s2) in [(&cyl_a, &cyl_b), (&cyl_b, &cyl_a)] {
        let result =
            build_fillet_surface(s1, s2, &cr, r, 5, true, ContinuityLevel::G1, false, false)
                .expect("parallel cylinder-cylinder fillet should build");
        match result.surface {
            SurfaceGeometry::Cylinder {
                axis,
                radius: cyl_radius,
                origin,
            } => {
                assert!(
                    axis.dot(&Vec3::new(0.0, 0.0, 1.0)).abs() > 0.999,
                    "cylinder axis should be ±ẑ, got {axis:?}"
                );
                assert!(
                    (cyl_radius - r).abs() < 1e-12,
                    "cylinder radius should equal fillet radius {r}, got {cyl_radius}"
                );
                // Origin should be on the center curve line at
                // (1.5, y_c, *).
                assert!(
                    (origin.x - 1.5).abs() < 0.05 && (origin.y - y_c).abs() < 0.05,
                    "cylinder origin should be near (1.5, {y_c}, *), got {origin:?}"
                );
            }
            other => panic!(
                "parallel cylinder-cylinder expected analytic Cylinder, got {:?}",
                std::mem::discriminant(&other)
            ),
        }
    }
}

/// Two cylinders with NON-parallel axes must not take the analytic
/// cylinder path. The axis-parallel gate rejects the pair before
/// sampling, and the result falls through to the NURBS sweep.
#[test]
fn fillet_surface_non_parallel_cylinder_cylinder_skips_analytic_path() {
    let cyl_a = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    // Second cylinder's axis is +x instead of +z.
    let cyl_b = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 1.0,
    };
    let cr = linear_center_result();
    let result = build_fillet_surface(
        &cyl_a,
        &cyl_b,
        &cr,
        0.3,
        5,
        true,
        ContinuityLevel::G1,
        false,
        false,
    )
    .expect("non-parallel cylinder-cylinder fillet should still build (NURBS sweep)");
    assert!(
        matches!(result.surface, SurfaceGeometry::Nurbs(_)),
        "non-parallel cylinder-cylinder should fall through to the NURBS sweep, \
         got {:?}",
        std::mem::discriminant(&result.surface)
    );
}

/// A plane-cylinder pair where the plane is NOT perpendicular to the
/// cylinder axis must NOT take the analytic torus path — the
/// intersection is an ellipse and a torus would be geometrically
/// wrong. The detector rejects the pair on the normal-vs-axis gate
/// before sampling the center curve.
#[test]
fn fillet_surface_plane_cylinder_oblique_skips_analytic_torus() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let c = std::f64::consts::FRAC_1_SQRT_2;
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, c, c),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, c, -c),
    };
    let cr = linear_center_result();
    let result = build_fillet_surface(
        &cyl,
        &plane,
        &cr,
        0.5,
        5,
        true,
        ContinuityLevel::G1,
        false,
        false,
    )
    .expect("oblique plane-cylinder fillet should still build (NURBS sweep)");
    assert!(
        matches!(result.surface, SurfaceGeometry::Nurbs(_)),
        "oblique plane-cylinder should fall through to the NURBS sweep, \
         got {:?}",
        std::mem::discriminant(&result.surface)
    );
}
