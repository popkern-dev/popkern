//! Plane×torus analytic-arm tests: the circle classes (perpendicular,
//! axis-containing, Villarceau), tangency-silent bands, dispatch
//! integration, and operand-order equivalence (PK-279).
use super::*;
use crate::geometry::intersection::surface_surface::intersect_surfaces;
use crate::geometry::surface::SurfaceEval;

fn torus() -> SurfaceGeometry {
    SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    }
}

fn plane_z(z: f64) -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, z),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn full() -> (f64, f64) {
    (0.0, TAU)
}

fn torus_residual(p: &Point3, big_r: f64, r: f64) -> f64 {
    // Both profile sheets: a spindle torus's inner lobe corresponds to
    // negative signed radius (the antipodal chart preimage).
    let rho = (p.x * p.x + p.y * p.y).sqrt();
    let outer = (((rho - big_r) * (rho - big_r) + p.z * p.z).sqrt() - r).abs();
    let inner = (((rho + big_r) * (rho + big_r) + p.z * p.z).sqrt() - r).abs();
    outer.min(inner)
}

fn curves_of(
    result: Vec<SurfaceSurfaceIntersection>,
) -> Vec<crate::geometry::curve_on_surface::CurveOnSurface> {
    result
        .into_iter()
        .filter_map(|ix| match ix {
            SurfaceSurfaceIntersection::Curve(c) => Some(c),
            SurfaceSurfaceIntersection::CoincidentRegion { .. } => None,
        })
        .collect()
}

/// Verify every curve lies on plane and torus and its pcurves are
/// co-parameterized with the 3D curve on both operands.
fn verify_curves(
    curves: &[crate::geometry::curve_on_surface::CurveOnSurface],
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    plane: &SurfaceGeometry,
    surface_gate: f64,
) {
    let SurfaceGeometry::Plane { origin, normal, .. } = plane else {
        panic!("plane expected");
    };
    for c in curves {
        for i in 0..=256 {
            let t = i as f64 / 256.0;
            let p = c.curve_3d.point_at(t);
            assert!(
                (p - *origin).dot(normal).abs() < surface_gate,
                "point off plane by {}",
                (p - *origin).dot(normal).abs()
            );
            assert!(
                torus_residual(&p, 2.0, 0.5) < surface_gate,
                "point off torus by {}",
                torus_residual(&p, 2.0, 0.5)
            );
            let uv1 = c.curve_uv1.point_at(t);
            let uv2 = c.curve_uv2.point_at(t);
            let sync1 = s1.point_at(uv1.x, uv1.y).distance_to(&p);
            let sync2 = s2.point_at(uv2.x, uv2.y).distance_to(&p);
            assert!(sync1 < 5e-6, "uv1 sync {sync1}");
            assert!(sync2 < 5e-6, "uv2 sync {sync2}");
        }
    }
}

/// The README-style usage example from the design review: a plane cut
/// through a fillet-like torus, straight from `intersect_surfaces`.
#[test]
fn perpendicular_plane_cuts_two_exact_circles() {
    let tol = Tolerance::default();
    let torus = torus();
    let plane = plane_z(0.0);
    let result = intersect_surfaces(
        &plane,
        (-4.0, 4.0),
        (-4.0, 4.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("intersect_surfaces");
    let curves = curves_of(result.intersections);
    assert_eq!(curves.len(), 2, "equatorial plane: outer + inner circle");
    for c in &curves {
        assert!(c.is_closed);
        assert_eq!(c.curve_3d.degree, 2, "exact rational circle expected");
    }
    verify_curves(&curves, &plane, &torus, &plane, 1e-9);
    let radii: Vec<f64> = curves
        .iter()
        .map(|c| {
            let p = c.curve_3d.point_at(0.0);
            (p.x * p.x + p.y * p.y).sqrt()
        })
        .collect();
    assert!(radii.iter().any(|r| (r - 2.5).abs() < 1e-9));
    assert!(radii.iter().any(|r| (r - 1.5).abs() < 1e-9));
}

#[test]
fn perpendicular_offset_plane_cuts_two_circles() {
    let tol = Tolerance::default();
    let torus = torus();
    let plane = plane_z(0.3);
    let out = try_plane_torus(
        &plane,
        (-4.0, 4.0),
        (-4.0, 4.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("plane×torus must be claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2);
    verify_curves(&curves, &plane, &torus, &plane, 1e-9);
}

#[test]
fn tangent_and_missing_perpendicular_planes_are_silent() {
    let tol = Tolerance::default();
    let torus = torus();
    for z in [0.5, -0.5, 0.7] {
        let plane = plane_z(z);
        let out = try_plane_torus(
            &plane,
            (-4.0, 4.0),
            (-4.0, 4.0),
            &torus,
            full(),
            full(),
            &tol,
        )
        .expect("plane×torus must be claimed");
        assert!(curves_of(out).is_empty(), "z={z} must be silent");
    }
}

#[test]
fn axis_containing_plane_cuts_two_tube_circles() {
    let tol = Tolerance::default();
    let torus = torus();
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let out = try_plane_torus(
        &plane,
        (-4.0, 4.0),
        (-4.0, 4.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2, "axis plane: two tube circles");
    verify_curves(&curves, &plane, &torus, &plane, 1e-9);
    for c in &curves {
        assert_eq!(c.curve_3d.degree, 2, "exact rational circle expected");
        let p = c.curve_3d.point_at(0.0);
        // Tube circles centred at (0, ±2, 0) with radius 0.5.
        let d_pos = ((p.y - 2.0) * (p.y - 2.0) + p.z * p.z + p.x * p.x).sqrt();
        let d_neg = ((p.y + 2.0) * (p.y + 2.0) + p.z * p.z + p.x * p.x).sqrt();
        assert!(
            (d_pos - 0.5).abs() < 1e-9 || (d_neg - 0.5).abs() < 1e-9,
            "point not on a tube circle"
        );
    }
}

#[test]
fn villarceau_plane_cuts_two_exact_circles_of_major_radius() {
    let tol = Tolerance::default();
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 3.0,
    };
    // Bitangent plane 3x = 4z through the centre: |n·a| = 4/5 =
    // √(R²−r²)/R.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.6, 0.0, -0.8),
        u_axis: Vec3::new(0.8, 0.0, 0.6),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let out = try_plane_torus(
        &plane,
        (-9.0, 9.0),
        (-9.0, 9.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2, "Villarceau plane: two circles");
    for c in &curves {
        assert_eq!(
            c.curve_3d.degree, 2,
            "the exact Villarceau arm must fire (not the sampled spiric arm)"
        );
        // Radius R = 5 about a centre at (0, ±3, 0).
        let p0 = c.curve_3d.point_at(0.0);
        let d_pos = (p0 - Point3::new(0.0, 3.0, 0.0)).length();
        let d_neg = (p0 - Point3::new(0.0, -3.0, 0.0)).length();
        assert!(
            (d_pos - 5.0).abs() < 1e-9 || (d_neg - 5.0).abs() < 1e-9,
            "not a Villarceau circle: d_pos={d_pos} d_neg={d_neg}"
        );
        for i in 0..=128 {
            let t = i as f64 / 128.0;
            let p = c.curve_3d.point_at(t);
            let rho = (p.x * p.x + p.y * p.y).sqrt();
            let res = (((rho - 5.0) * (rho - 5.0) + p.z * p.z).sqrt() - 3.0).abs();
            assert!(res < 1e-9, "Villarceau sample off torus by {res}");
        }
    }
}

#[test]
fn operand_order_is_equivalent_for_the_circle_classes() {
    let tol = Tolerance::default();
    let torus = torus();
    let plane = plane_z(0.3);
    let ab = curves_of(
        try_plane_torus(
            &plane,
            (-4.0, 4.0),
            (-4.0, 4.0),
            &torus,
            full(),
            full(),
            &tol,
        )
        .unwrap(),
    );
    let ba = curves_of(
        try_plane_torus(
            &torus,
            full(),
            full(),
            &plane,
            (-4.0, 4.0),
            (-4.0, 4.0),
            &tol,
        )
        .unwrap(),
    );
    assert_eq!(ab.len(), ba.len());
    verify_curves(&ba, &torus, &plane, &plane, 1e-9);
    for a in &ab {
        let pa = a.curve_3d.point_at(0.41);
        let best = ba
            .iter()
            .map(|b| {
                (0..=4096)
                    .map(|i| b.curve_3d.point_at(i as f64 / 4096.0).distance_to(&pa))
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(f64::INFINITY, f64::min);
        // Gate at the sampling spacing of the densely-sampled circle.
        assert!(best < 5e-3, "curve sets differ between operand orders");
    }
}

#[test]
fn spindle_torus_perpendicular_plane_keeps_only_real_circles() {
    let tol = Tolerance::default();
    // Spindle torus: minor > major. The equatorial plane still cuts
    // the outer circle R+r and the inner lobe |R−r| (antipodal side).
    let spindle = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 0.5,
        minor_radius: 0.8,
    };
    let plane = plane_z(0.0);
    let out = try_plane_torus(
        &plane,
        (-4.0, 4.0),
        (-4.0, 4.0),
        &spindle,
        full(),
        full(),
        &tol,
    )
    .expect("claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2);
    for c in &curves {
        let p = c.curve_3d.point_at(0.125);
        let rho = (p.x * p.x + p.y * p.y).sqrt();
        let res = (((rho - 0.5) * (rho - 0.5) + p.z * p.z).sqrt() - 0.8)
            .abs()
            .min((((rho + 0.5) * (rho + 0.5) + p.z * p.z).sqrt() - 0.8).abs());
        assert!(res < 1e-9, "spindle circle off torus by {res}");
    }
}
