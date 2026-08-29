//! §11.1 analytic fast-path tests for `build_fillet_surface`.
//!
//! Split from `fillet_surface_tests.rs` so neither file exceeds
//! the 500-line structural cap. Covers the surface pairs whose
//! rolling-ball fillet is an exact torus or cylinder:
//!
//! - plane-sphere (always circular center curve)
//! - plane-cylinder with axis perpendicular to the plane
//! - plane-cone with axis perpendicular to the plane
//! - cylinder-cylinder with parallel axes
//!
//! Plus companion negative tests that non-perpendicular / non-
//! parallel configurations fall through to the generic NURBS sweep.

use super::*;
use crate::center_curve::CenterCurveResult;
use crate::types::ContinuityLevel;
use popkern_core::geometry::{Point3, Vec3};

/// Small stand-in `CenterCurveResult` for the oblique negative
/// test — the detector rejects the surface pair before sampling,
/// so the center curve's contents don't matter.
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
        sweep_center_curve: None,
        sweep_footprint_1: None,
        sweep_footprint_2: None,
        center_curve,
        footprint_1: fp1,
        footprint_2: fp2,
        contact_solutions: None,
    }
}

/// §11.1 plane-sphere analytic fast path: a constant-radius fillet
/// between a plane and a sphere must produce an exact Torus, not a
/// generic NURBS sweep.
///
/// Configuration: sphere of radius 5 centred at origin, plane z = 3
/// so the intersection is a circle. For a fillet radius r = 0.5 the
/// offset plane is at z = 3 − 0.5 = 2.5 and the offset sphere has
/// radius 5.5, so the center curve is a circle at z = 2.5 with
/// radius √(5.5² − 2.5²) = √24 ≈ 4.899.
///
/// We don't run the real center-curve producer here (that needs the
/// offset-surface intersection solver); we hand-build a
/// `CenterCurveResult` with the analytically-known circle and
/// verify `build_fillet_surface` picks the analytic path.
#[test]
fn fillet_surface_plane_sphere_produces_analytic_torus() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 5.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 3.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // Center curve: circle in z = 2.5 plane, radius √24, centred at origin.
    // Use 60 input points so the cubic NURBS fit-+-resample error
    // falls to ~0.1%, below the test tolerance.
    let major_expected = (24.0_f64).sqrt();
    let n = 60usize;
    let center_pts: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            Point3::new(
                major_expected * theta.cos(),
                major_expected * theta.sin(),
                2.5,
            )
        })
        .collect();
    let center_curve = crate::center_curve::fit_nurbs_through_points(&center_pts).unwrap();
    let fp1: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            Point3::new(
                major_expected * theta.cos(),
                major_expected * theta.sin(),
                3.0,
            )
        })
        .collect();
    let fp2: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            let r_in = major_expected - 0.1;
            Point3::new(r_in * theta.cos(), r_in * theta.sin(), 3.0)
        })
        .collect();
    let cr = CenterCurveResult {
        sweep_center_curve: None,
        sweep_footprint_1: None,
        sweep_footprint_2: None,
        center_curve,
        footprint_1: fp1,
        footprint_2: fp2,
        contact_solutions: None,
    };

    let r = 0.5;
    for (s1, s2) in [(&sphere, &plane), (&plane, &sphere)] {
        let result = build_fillet_surface(
            s1,
            s2,
            &cr,
            r,
            5,
            true,
            ContinuityLevel::G1,
            false,
            false,
            crate::types::BlendProfile::Circular,
        )
        .expect("plane-sphere fillet should build");
        match result.surface {
            SurfaceGeometry::Torus {
                axis,
                major_radius,
                minor_radius,
                center,
            } => {
                assert!(
                    axis.dot(&Vec3::new(0.0, 0.0, 1.0)).abs() > 0.999,
                    "torus axis should be ±ẑ, got {axis:?}"
                );
                // Minor radius equals the fillet radius exactly —
                // this value is passed through, not sampled.
                assert!(
                    (minor_radius - r).abs() < 1e-12,
                    "minor_radius should equal fillet radius {r}, got {minor_radius}"
                );
                // Major radius and centre are read from sampling
                // the NURBS-fitted center curve, so they inherit
                // hand-fit interpolation error. Loose 2 % tolerance
                // — the real pipeline has tighter fidelity.
                assert!(
                    (major_radius - major_expected).abs() / major_expected < 2e-2,
                    "major_radius expected {major_expected}, got {major_radius}"
                );
                assert!(
                    (center.z - 2.5).abs() < 0.1,
                    "centre z expected 2.5, got {}",
                    center.z
                );
            }
            other => panic!(
                "plane-sphere expected analytic Torus, got {:?}",
                std::mem::discriminant(&other)
            ),
        }
    }
}

/// §11.1 plane-cylinder (perpendicular axis) analytic fast path.
/// A plane perpendicular to a cylinder's axis has a circular edge
/// and the rolling-ball fillet along that edge is an exact torus
/// with the cylinder's axis.
#[test]
fn fillet_surface_plane_cylinder_perpendicular_produces_analytic_torus() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // Center curve: rolling-ball centers lie on a circle at
    // z = 1 − r, with radial distance (cyl_radius − r) from
    // the axis. For r = 0.3: z ≈ 0.7, R ≈ 1.7.
    let r = 0.3;
    let major_expected = 2.0 - r;
    let n = 60usize;
    let center_pts: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            Point3::new(
                major_expected * theta.cos(),
                major_expected * theta.sin(),
                1.0 - r,
            )
        })
        .collect();
    let center_curve = crate::center_curve::fit_nurbs_through_points(&center_pts).unwrap();
    let fp1: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            Point3::new(2.0 * theta.cos(), 2.0 * theta.sin(), 1.0 - r)
        })
        .collect();
    let fp2: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            Point3::new(
                major_expected * theta.cos(),
                major_expected * theta.sin(),
                1.0,
            )
        })
        .collect();
    let cr = CenterCurveResult {
        sweep_center_curve: None,
        sweep_footprint_1: None,
        sweep_footprint_2: None,
        center_curve,
        footprint_1: fp1,
        footprint_2: fp2,
        contact_solutions: None,
    };

    for (s1, s2) in [(&cyl, &plane), (&plane, &cyl)] {
        let result = build_fillet_surface(
            s1,
            s2,
            &cr,
            r,
            5,
            true,
            ContinuityLevel::G1,
            false,
            false,
            crate::types::BlendProfile::Circular,
        )
        .expect("plane-cylinder (perp) fillet should build");
        match result.surface {
            SurfaceGeometry::Torus {
                axis, minor_radius, ..
            } => {
                assert!(
                    axis.dot(&Vec3::new(0.0, 0.0, 1.0)).abs() > 0.999,
                    "torus axis should be ±ẑ, got {axis:?}"
                );
                assert!(
                    (minor_radius - r).abs() < 1e-12,
                    "minor_radius should equal {r}, got {minor_radius}"
                );
            }
            other => panic!(
                "perpendicular plane-cylinder expected analytic Torus, got {:?}",
                std::mem::discriminant(&other)
            ),
        }
    }
}

/// §11.1 plane-cone (perpendicular axis) analytic fast path. A
/// plane perpendicular to a cone's axis cuts the cone in a circle;
/// by rotational symmetry the rolling-ball center curve is also a
/// 3D circle, and the fillet surface is an exact torus.
///
/// Configuration: cone with apex at (0, 0, 2) opening downward with
/// half-angle arctan(0.5) (base radius 1 at z = 0), plane z = 0.
/// For a fillet radius r = 0.1, the rolling-ball 2D tangency in the
/// (radial, axial) half-plane gives the ball center at some (ρ, z)
/// with ρ slightly less than 1 and z slightly above 0 — exact value
/// depends on the cone half-angle, but by rotational symmetry the
/// 3D locus is a horizontal circle at that (ρ, z).
///
/// Instead of deriving ρ and z analytically, hand-build a plausible
/// circle and verify `build_fillet_surface` takes the analytic path.
/// The end-to-end test in `contact_cone_demo_tests.rs` drives the
/// real walker → analytic-torus pipeline on `BRep::cone`.
#[test]
fn fillet_surface_plane_cone_perpendicular_produces_analytic_torus() {
    let half_angle = 0.5_f64.atan();
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 2.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // EXACT rolling-ball geometry for the base fillet of a solid
    // cone standing on the plane (concave groove, ball in the
    // void): ball centers sit at z = r and at radial distance
    // ρ_cone(z=r) + r/cos α from the axis, where ρ_cone(z) =
    // (2 − z)·tan α is the cone's section radius. The exact-torus
    // detector derives this circle from the surface pair and
    // verifies the spine/footprints actually lie on it, so
    // hand-plausible-but-inconsistent values now (correctly)
    // decline to the NURBS sweep.
    let r = 0.1;
    let z = r;
    let cos_a = half_angle.cos();
    let sin_a = half_angle.sin();
    let rho_cone_at_z = (2.0 - z) * half_angle.tan();
    let rho = rho_cone_at_z + r / cos_a;
    let n = 60usize;
    let center_pts: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            Point3::new(rho * theta.cos(), rho * theta.sin(), z)
        })
        .collect();
    let center_curve = crate::center_curve::fit_nurbs_through_points(&center_pts).unwrap();
    // Exact contacts: on the cone, r along the inward surface
    // normal (cos α radial, sin α axial) from the ball center; on
    // the plane, directly below the center.
    let fp1: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            let rc = rho - r * cos_a;
            Point3::new(rc * theta.cos(), rc * theta.sin(), z - r * sin_a)
        })
        .collect();
    let fp2: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            Point3::new(rho * theta.cos(), rho * theta.sin(), 0.0)
        })
        .collect();
    let cr = CenterCurveResult {
        sweep_center_curve: None,
        sweep_footprint_1: None,
        sweep_footprint_2: None,
        center_curve,
        footprint_1: fp1,
        footprint_2: fp2,
        contact_solutions: None,
    };

    for (s1, s2) in [(&cone, &plane), (&plane, &cone)] {
        let result = build_fillet_surface(
            s1,
            s2,
            &cr,
            r,
            5,
            true,
            ContinuityLevel::G1,
            false,
            false,
            crate::types::BlendProfile::Circular,
        )
        .expect("plane-cone (perp) fillet should build");
        match result.surface {
            SurfaceGeometry::Torus {
                axis, minor_radius, ..
            } => {
                assert!(
                    axis.dot(&Vec3::new(0.0, 0.0, 1.0)).abs() > 0.999,
                    "torus axis should be ±ẑ (parallel to cone axis), got {axis:?}"
                );
                assert!(
                    (minor_radius - r).abs() < 1e-12,
                    "minor_radius should equal fillet radius {r}, got {minor_radius}"
                );
            }
            other => panic!(
                "perpendicular plane-cone expected analytic Torus, got {:?}",
                std::mem::discriminant(&other)
            ),
        }
    }
}

/// A plane-cone pair where the plane is NOT perpendicular to the
/// cone axis must NOT take the analytic torus path — the
/// intersection is an ellipse / parabola / hyperbola and the
/// rolling-ball fillet is a canal surface with no closed-form
/// popkern primitive. Verify the detector bails to the NURBS
/// sweep.
#[test]
fn fillet_surface_plane_cone_oblique_skips_analytic_torus() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 2.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: 0.5_f64.atan(),
    };
    // Plane normal tilted 30° from the cone axis.
    let c30 = 30f64.to_radians().cos();
    let s30 = 30f64.to_radians().sin();
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, s30, c30),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, c30, -s30),
    };
    let cr = linear_center_result();
    let result = build_fillet_surface(
        &cone,
        &plane,
        &cr,
        0.1,
        5,
        true,
        ContinuityLevel::G1,
        false,
        false,
        crate::types::BlendProfile::Circular,
    )
    .expect("oblique plane-cone fillet should still build (NURBS sweep)");
    assert!(
        matches!(result.surface, SurfaceGeometry::Nurbs(_)),
        "oblique plane-cone should fall through to the NURBS sweep, got {:?}",
        std::mem::discriminant(&result.surface)
    );
}
