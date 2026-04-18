use super::*;
use crate::geometry::surface::SurfaceEval;

#[test]
fn plane_is_not_periodic() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert_eq!(plane.is_periodic(), (false, false));
}

#[test]
fn cylinder_is_periodic_in_u_only() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert_eq!(cyl.is_periodic(), (true, false));
}

#[test]
fn cone_is_periodic_in_u_only() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 1.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    };
    assert_eq!(cone.is_periodic(), (true, false));
}

#[test]
fn sphere_is_periodic_in_u_only() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    assert_eq!(sphere.is_periodic(), (true, false));
}

#[test]
fn torus_is_periodic_in_both() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    assert_eq!(torus.is_periodic(), (true, true));
}

#[test]
fn periodicity_returns_correct_periods() {
    use std::f64::consts::TAU;
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert_eq!(plane.periodicity(), (None, None));

    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert_eq!(cyl.periodicity(), (Some(TAU), None));

    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    assert_eq!(torus.periodicity(), (Some(TAU), Some(TAU)));
}

#[test]
fn nurbs_is_not_periodic() {
    let nurbs = crate::geometry::nurbs_surface::NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };
    let geom = SurfaceGeometry::Nurbs(nurbs);
    assert_eq!(geom.is_periodic(), (false, false));
}

#[test]
fn plane_to_nurbs_preserves_shape() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let nurbs = plane.to_nurbs_in_range((-1.0, 1.0), (-1.0, 1.0));
    for u in [-1.0_f64, -0.5, 0.0, 0.5, 0.99] {
        for v in [-1.0_f64, -0.5, 0.0, 0.5, 0.99] {
            let p1 = plane.point_at(u, v);
            let p2 = nurbs.point_at(u, v);
            assert!(
                p1.distance_to(&p2) < 1e-8,
                "Plane-to-NURBS mismatch at u={u}, v={v}: p1={p1:?}, p2={p2:?}"
            );
        }
    }
}

#[test]
fn nurbs_to_nurbs_is_clone() {
    let nurbs = NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };
    let geom = SurfaceGeometry::Nurbs(nurbs.clone());
    let promoted = geom.to_nurbs_in_range((0.0, 1.0), (0.0, 1.0));
    assert_eq!(promoted.degree_u, nurbs.degree_u);
}

#[test]
fn sphere_to_nurbs_approximate() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let nurbs = sphere.to_nurbs_in_range((0.5, 2.5), (0.5, 2.5));
    for u in [0.6_f64, 1.0, 1.5, 2.0, 2.4] {
        for v in [0.6_f64, 1.0, 1.5, 2.0, 2.4] {
            let p = nurbs.point_at(u, v);
            let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
            assert!(
                (r - 1.0).abs() < 0.1,
                "Point not on unit sphere: r={r} at u={u}, v={v}"
            );
        }
    }
}

#[test]
fn sphere_to_nurbs_default_range() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let nurbs = sphere.to_nurbs();
    let p = nurbs.point_at(std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
    let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
    assert!((r - 1.0).abs() < 0.15, "Point not on unit sphere: r={r}");
}

#[test]
fn cylinder_to_nurbs_in_range() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let nurbs = cyl.to_nurbs_in_range((0.0, std::f64::consts::TAU), (-1.0, 1.0));
    let p = nurbs.point_at(std::f64::consts::FRAC_PI_2, 0.0);
    let r = (p.x * p.x + p.y * p.y).sqrt();
    assert!((r - 1.0).abs() < 0.15, "Point not on cylinder: r={r}");
}

/// §5.3 item #9 — Gregory sub-range `to_nurbs_in_range` must
/// preserve the exact rational form (cropped from the full-range
/// exact NURBS via `split_at_u` / `split_at_v`). Previously the
/// sub-range path fell through to `sample_to_nurbs`, a lossy
/// degree-3 interpolation. The cropped NURBS's point evaluations
/// must match the Gregory point evaluations to near-f64 precision
/// within the cropped domain.
#[test]
fn gregory_to_nurbs_sub_range_preserves_exact_form() {
    use crate::geometry::gregory::{GregoryPatch, RectangularGregory};

    // Build a twisted Gregory patch (non-planar, non-trivial inner CPs).
    let lin = |u: f64, v: f64| Point3::new(u, v, 0.0);
    let mut boundary_cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (i, row) in boundary_cps.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            *cp = lin(i as f64 / 3.0, j as f64 / 3.0);
        }
    }
    // Lift the top-right corner to create a twist.
    boundary_cps[3][3] = Point3::new(1.0, 1.0, 0.3);
    boundary_cps[2][3] = Point3::new(2.0 / 3.0, 1.0, 0.2);
    boundary_cps[3][2] = Point3::new(1.0, 2.0 / 3.0, 0.2);
    let inner_cps = [
        (
            Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.05),
            Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.02),
        ),
        (
            Point3::new(1.0 / 3.0, 2.0 / 3.0, 0.08),
            Point3::new(1.0 / 3.0, 2.0 / 3.0, 0.10),
        ),
        (
            Point3::new(2.0 / 3.0, 2.0 / 3.0, 0.15),
            Point3::new(2.0 / 3.0, 2.0 / 3.0, 0.12),
        ),
        (
            Point3::new(2.0 / 3.0, 1.0 / 3.0, 0.07),
            Point3::new(2.0 / 3.0, 1.0 / 3.0, 0.04),
        ),
    ];
    let rg = RectangularGregory {
        boundary_cps,
        inner_cps,
    };
    let geom = SurfaceGeometry::Gregory(GregoryPatch::Rectangular(Box::new(rg.clone())));

    // Sub-range request: interior rectangle [0.2, 0.8] × [0.3, 0.7].
    let (u_lo, u_hi) = (0.2, 0.8);
    let (v_lo, v_hi) = (0.3, 0.7);
    let cropped = geom.to_nurbs_in_range((u_lo, u_hi), (v_lo, v_hi));

    // The cropped NURBS's parameter domain is [u_lo, u_hi] × [v_lo, v_hi].
    // Evaluate at a 5×5 grid inside that domain and compare against the
    // Gregory patch's native evaluation at the same (u, v).
    let n = 5;
    let mut max_dev: f64 = 0.0;
    for i in 0..=n {
        for j in 0..=n {
            let u = u_lo + (u_hi - u_lo) * (i as f64) / (n as f64);
            let v = v_lo + (v_hi - v_lo) * (j as f64) / (n as f64);
            let p_gregory = rg.point_at(u, v);
            let p_nurbs = cropped.point_at(u, v);
            let d = (p_gregory.x - p_nurbs.x)
                .abs()
                .max((p_gregory.y - p_nurbs.y).abs())
                .max((p_gregory.z - p_nurbs.z).abs());
            max_dev = max_dev.max(d);
        }
    }
    // Exact NURBS + split is lossless modulo float noise; sampling
    // was ~1e-3 accurate at best. 1e-6 cleanly separates the two.
    assert!(
        max_dev < 1e-6,
        "§5.3 item #9 sub-range exact crop: max deviation {max_dev} between Gregory \
         eval and cropped NURBS eval must be < 1e-6 (old sample-fit path was ~1e-3)"
    );
}
