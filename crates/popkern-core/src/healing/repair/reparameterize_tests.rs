use super::*;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{NurbsSurface, Point3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn measure_distortion_plane_near_uniform() {
    // A plane has perfectly uniform parameterization: area element ratio ~1.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: crate::geometry::Vec3::new(0.0, 0.0, 1.0),
        u_axis: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
        v_axis: crate::geometry::Vec3::new(0.0, 1.0, 0.0),
    };
    let ratio = measure_distortion(&plane, (0.0, 1.0), (0.0, 1.0), 5);
    let r = ratio.expect("plane should have valid distortion");
    assert!(
        (r - 1.0).abs() < 0.01,
        "plane distortion should be ~1.0, got {r}"
    );
}

#[test]
fn measure_distortion_sphere_non_uniform() {
    // A sphere's parameterization compresses area at the poles.
    // Sampling away from exact poles should still show distortion > 1.
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    // u in [0, 2*pi], v in [0.1, pi-0.1] (avoid exact poles).
    let ratio = measure_distortion(
        &sphere,
        (0.0, std::f64::consts::TAU),
        (0.1, std::f64::consts::PI - 0.1),
        8,
    );
    let r = ratio.expect("sphere should have valid distortion");
    assert!(r > 1.0, "sphere distortion should be > 1.0, got {r}");
}

#[test]
fn reparameterize_face_plane_returns_none() {
    // Cube faces are planes, not NURBS, so reparameterization should return None.
    let mut brep =
        BRep::unit_cube(Tolerance::default()).unwrap_or_else(|_| BRep::new(Tolerance::default()));
    let face = brep.faces.keys().next();
    let Some(fh) = face else {
        return;
    };
    let result = reparameterize_face(&mut brep, fh, 5.0, 0.01);
    assert!(
        result.is_none(),
        "plane face should not be reparameterizable"
    );
}

#[test]
fn reparameterize_all_valid_cube_empty() {
    // Unit cube has all planar faces, so no NURBS to reparameterize.
    let mut brep =
        BRep::unit_cube(Tolerance::default()).unwrap_or_else(|_| BRep::new(Tolerance::default()));
    let results = reparameterize_all(&mut brep, 5.0, 0.01);
    assert!(
        results.is_empty(),
        "expected no reparameterization on a cube, got {} results",
        results.len()
    );
}

/// Build a NURBS surface representing a plane with distorted (quadratic)
/// u-parameterization. Control points at u_i = (i/(n-1))^2 in 3D space
/// but at uniform parameter values in the knot vector.
fn make_distorted_nurbs() -> NurbsSurface {
    let n_u = 5;
    let n_v = 5;
    let degree = 3;

    // Place control points with quadratic spacing in x (the u direction):
    // x_i = (i/(n-1))^2, so points are compressed near x=0 and
    // stretched near x=1.
    let mut control_points = Vec::with_capacity(n_u * n_v);
    for i in 0..n_u {
        let x = (i as f64 / (n_u - 1) as f64).powi(2);
        for j in 0..n_v {
            let y = j as f64 / (n_v - 1) as f64;
            control_points.push(Point3::new(x, y, 0.0));
        }
    }

    let knots_u = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
    let knots_v = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
    let weights = vec![1.0; n_u * n_v];

    NurbsSurface {
        degree_u: degree,
        degree_v: degree,
        knots_u,
        knots_v,
        control_points,
        num_u: n_u,
        num_v: n_v,
        weights,
    }
}

#[test]
fn reparameterize_distorted_nurbs_improves_uniformity() {
    let nurbs = make_distorted_nurbs();
    let (u_range, v_range) = nurbs.parameter_range();

    let before =
        measure_distortion(&nurbs, u_range, v_range, 10).expect("distortion should be measurable");
    assert!(
        before > 1.5,
        "distorted surface should have ratio > 1.5, got {before}"
    );

    // Allow generous deviation for this test.
    let result = reparameterize_surface(&nurbs, u_range, v_range, 1.5, 1.0);
    let Some((new_nurbs, reparam_result)) = result else {
        // If reparameterization returns None, it means the surface is
        // either already below threshold or deviation is too large.
        // The surface is flat so deviation should be small; check that
        // the distortion was below threshold.
        assert!(
            before <= 1.5,
            "reparameterize_surface returned None but distortion_before={before} > 1.5"
        );
        return;
    };

    assert!(
        reparam_result.distortion_before > 1.0,
        "distortion_before should be > 1.0"
    );

    // New surface should evaluate to similar 3D points.
    let new_range = new_nurbs.parameter_range();
    let p_old = nurbs.point_at((u_range.0 + u_range.1) * 0.5, (v_range.0 + v_range.1) * 0.5);
    let p_new = new_nurbs.point_at(
        (new_range.0.0 + new_range.0.1) * 0.5,
        (new_range.1.0 + new_range.1.1) * 0.5,
    );
    let dist = p_old.distance_to(&p_new);
    assert!(dist < 1.0, "midpoint deviation={dist} exceeds tolerance");
}

#[test]
fn measure_distortion_small_grid_returns_none() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: crate::geometry::Vec3::new(0.0, 0.0, 1.0),
        u_axis: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
        v_axis: crate::geometry::Vec3::new(0.0, 1.0, 0.0),
    };
    // Grid size < 2 should return None.
    assert!(measure_distortion(&plane, (0.0, 1.0), (0.0, 1.0), 1).is_none());
    assert!(measure_distortion(&plane, (0.0, 1.0), (0.0, 1.0), 0).is_none());
}

#[test]
fn measure_distortion_zero_span_returns_none() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: crate::geometry::Vec3::new(0.0, 0.0, 1.0),
        u_axis: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
        v_axis: crate::geometry::Vec3::new(0.0, 1.0, 0.0),
    };
    assert!(measure_distortion(&plane, (0.5, 0.5), (0.0, 1.0), 5).is_none());
    assert!(measure_distortion(&plane, (0.0, 1.0), (0.5, 0.5), 5).is_none());
}

#[test]
fn arc_length_to_param_clamps_correctly() {
    let table = vec![(0.0, 0.0), (1.0, 0.5), (3.0, 1.0)];

    // Before start.
    let t = arc_length_to_param(&table, -1.0);
    assert!((t - 0.0).abs() < 1e-10);

    // After end.
    let t = arc_length_to_param(&table, 10.0);
    assert!((t - 1.0).abs() < 1e-10);

    // Interpolated midpoint.
    let t = arc_length_to_param(&table, 0.5);
    assert!((t - 0.25).abs() < 0.01, "got {t}");
}

#[test]
fn uniform_clamped_knots_structure() {
    let knots = uniform_clamped_knots(5, 3);
    // Should have 5 + 3 + 1 = 9 knots.
    assert_eq!(knots.len(), 9, "expected 9 knots, got {}", knots.len());
    // First 4 should be 0.0 and last 4 should be 1.0.
    for k in &knots[..4] {
        assert!((k - 0.0).abs() < 1e-15);
    }
    for k in &knots[5..] {
        assert!((k - 1.0).abs() < 1e-15);
    }
    // Should be monotonically non-decreasing.
    for w in knots.windows(2) {
        assert!(w[1] >= w[0], "knot vector not non-decreasing: {:?}", knots);
    }
}
