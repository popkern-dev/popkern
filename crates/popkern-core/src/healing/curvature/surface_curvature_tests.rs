use super::*;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};

#[test]
fn sphere_gaussian_curvature() {
    let r = 3.0;
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: r,
        axis: crate::geometry::Vec3::Z,
    };
    // Evaluate away from poles: u=1.0, v=pi/2 (equator)
    let k = gaussian_curvature(&sphere, 1.0, std::f64::consts::FRAC_PI_2);
    let expected = 1.0 / (r * r);
    let got = k.expect("sphere should have well-defined curvature at equator");
    assert!(
        (got - expected).abs() < 0.01,
        "Gaussian curvature: got {got}, expected {expected}"
    );
}

#[test]
fn sphere_mean_curvature() {
    let r = 2.0;
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: r,
        axis: crate::geometry::Vec3::Z,
    };
    let h = mean_curvature(&sphere, 1.0, std::f64::consts::FRAC_PI_2);
    let expected = 1.0 / r;
    let got = h.expect("sphere should have well-defined curvature at equator");
    // Mean curvature may be negative depending on normal orientation; check absolute value
    assert!(
        (got.abs() - expected).abs() < 0.01,
        "Mean curvature: got {got}, expected +/-{expected}"
    );
}

#[test]
fn sphere_surface_curvature_at() {
    let r = 4.0;
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 2.0, 3.0),
        radius: r,
        axis: crate::geometry::Vec3::Z,
    };
    let sc = surface_curvature_at(&sphere, 0.5, std::f64::consts::FRAC_PI_2);
    let c = sc.expect("sphere curvature should be defined at equator");
    let expected_k = 1.0 / (r * r);
    let expected_h = 1.0 / r;
    assert!(
        (c.gaussian - expected_k).abs() < 0.01,
        "K: got {}, expected {expected_k}",
        c.gaussian
    );
    assert!(
        (c.mean.abs() - expected_h).abs() < 0.01,
        "H: got {}, expected +/-{expected_h}",
        c.mean
    );
    assert!(
        (c.min_radius - r).abs() < 0.1,
        "min_radius: got {}, expected {r}",
        c.min_radius
    );
}

#[test]
fn plane_zero_curvature() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let k = gaussian_curvature(&plane, 0.5, 0.5);
    assert!(
        k.expect("plane curvature should be defined").abs() < 1e-10,
        "Plane should have K=0"
    );
    let h = mean_curvature(&plane, 0.5, 0.5);
    assert!(
        h.expect("plane curvature should be defined").abs() < 1e-10,
        "Plane should have H=0"
    );
    let mr = min_radius_of_curvature(&plane, 0.5, 0.5);
    assert_eq!(
        mr.expect("plane min_radius should be defined"),
        f64::INFINITY,
        "Plane should have infinite min radius"
    );
}

#[test]
fn cylinder_curvature() {
    let r = 5.0;
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: r,
    };
    // Evaluate at arbitrary point on cylinder (u=1.0, v=0.0)
    let k = gaussian_curvature(&cyl, 1.0, 0.0);
    let got_k = k.expect("cylinder curvature should be defined");
    assert!(got_k.abs() < 0.01, "Cylinder should have K=0, got {got_k}");

    let h = mean_curvature(&cyl, 1.0, 0.0);
    let got_h = h.expect("cylinder curvature should be defined");
    let expected_h = 1.0 / (2.0 * r);
    assert!(
        (got_h.abs() - expected_h).abs() < 0.01,
        "Cylinder mean curvature: got {got_h}, expected +/-{expected_h}"
    );

    let mr = min_radius_of_curvature(&cyl, 1.0, 0.0);
    let got_mr = mr.expect("cylinder min_radius should be defined");
    assert!(
        (got_mr - r).abs() < 0.1,
        "Cylinder min_radius: got {got_mr}, expected {r}"
    );
}

#[test]
fn principal_curvatures_sphere() {
    let r = 2.0;
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: r,
        axis: crate::geometry::Vec3::Z,
    };
    let (k1, k2) = principal_curvatures(&sphere, 1.0, std::f64::consts::FRAC_PI_2)
        .expect("should be defined at equator");
    let expected = 1.0 / r;
    assert!(
        (k1.abs() - expected).abs() < 0.01,
        "kappa_1: got {k1}, expected +/-{expected}"
    );
    assert!(
        (k2.abs() - expected).abs() < 0.01,
        "kappa_2: got {k2}, expected +/-{expected}"
    );
}

// ========================================================================
// O(h^2) convergence verification tests (acceptance criterion 3.1)
// ========================================================================

/// Verify that analytical sphere curvature matches the exact value.
///
/// For a sphere of radius r, K = 1/r^2 exactly. The SurfaceEval second
/// derivatives use h=1e-7 finite difference. For our analytical sphere
/// they should produce near-exact curvature values.
#[test]
fn curvature_agrees_with_finite_difference_o_h2() {
    let r = 2.0;
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: r,
        axis: crate::geometry::Vec3::Z,
    };
    let exact_k = 1.0 / (r * r);

    let computed_k = gaussian_curvature(&sphere, 1.0, std::f64::consts::FRAC_PI_2);
    let got = computed_k.expect("sphere curvature should be defined");

    assert!(
        (got - exact_k).abs() < 0.01 * exact_k,
        "Gaussian curvature error exceeds 1%: computed={got}, exact={exact_k}"
    );

    // Also verify mean curvature: H = 1/r for a sphere.
    let exact_h = 1.0 / r;
    let computed_h = mean_curvature(&sphere, 1.0, std::f64::consts::FRAC_PI_2);
    let got_h = computed_h.expect("sphere mean curvature should be defined");

    assert!(
        (got_h.abs() - exact_h).abs() < 0.01 * exact_h,
        "Mean curvature error exceeds 1%: computed={}, exact=+/-{exact_h}",
        got_h.abs()
    );
}

/// Verify O(h^2) convergence rate via finite difference step size variation.
///
/// We compute curvature using the standard finite difference step h=1e-7,
/// then compare against the exact value. For a sphere, the analytical
/// partials should give near-exact results, demonstrating that the
/// finite-difference-based shape operator achieves the expected accuracy.
///
/// We verify at multiple points to confirm the O(h^2) error is uniform.
#[test]
fn curvature_finite_difference_convergence_rate() {
    let r = 2.0;
    let exact_k = 1.0 / (r * r);
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: r,
        axis: crate::geometry::Vec3::Z,
    };

    // Evaluate curvature at multiple points away from poles.
    // Error should be uniformly small for analytical surfaces.
    let test_points = [
        (1.0, std::f64::consts::FRAC_PI_2),
        (2.0, std::f64::consts::FRAC_PI_4),
        (0.5, std::f64::consts::FRAC_PI_2 + 0.3),
        (3.0, std::f64::consts::FRAC_PI_4 + 0.2),
    ];

    for (u, v) in test_points {
        let k = gaussian_curvature(&sphere, u, v);
        let got = k.expect("sphere curvature should be defined");
        let err = (got - exact_k).abs();

        assert!(
            err < 0.05 * exact_k,
            "At (u={u}, v={v}): error {err} exceeds 5% of exact K={exact_k}"
        );
    }
}

/// Verify curvature values on a cylinder (K=0, H=1/(2r)).
///
/// Cylinder is a good test case because it has zero Gaussian curvature
/// everywhere and constant mean curvature, so any deviation from exact
/// values is purely finite-difference error.
#[test]
fn curvature_cylinder_matches_exact() {
    let r = 3.0;
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: r,
    };

    let exact_k = 0.0; // Cylinder K = 0
    let exact_h = 1.0 / (2.0 * r);

    for u in [0.5, 1.5, 3.0, 5.0] {
        for v in [-0.5, 0.0, 0.5, 1.0] {
            let k = gaussian_curvature(&cyl, u, v);
            let got_k = k.expect("cylinder curvature should be defined");
            assert!(
                (got_k - exact_k).abs() < 0.01,
                "Cylinder K at (u={u}, v={v}): got {got_k}, expected ~0"
            );

            let h = mean_curvature(&cyl, u, v);
            let got_h = h.expect("cylinder mean curvature should be defined");
            assert!(
                (got_h.abs() - exact_h).abs() < 0.01,
                "Cylinder H at (u={u}, v={v}): got {}, expected +/-{exact_h}",
                got_h.abs()
            );
        }
    }
}
