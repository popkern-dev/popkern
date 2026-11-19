//! Trig sphere pairs on a turned sphere (PK-858): a configuration solved
//! upright and the same configuration under a rigid rotation (sphere
//! pole off world +Z) both serve the analytic dispatch with the
//! identity parameter map and the factory's guarantees, and the turned
//! triple's 3-D curve is the rotated upright one.

use super::tests_support::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::edge_triple::{
    EdgeAnchors, EdgeTripleOptions, EdgeTripleProvenance, ParameterMap, compute_edge_geometry,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::transform::Transform3;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

fn turn() -> Transform3 {
    Transform3::rotation_around_axis(Point3::new(0.3, -0.2, 0.5), Vec3::new(0.6, 0.0, 0.8), 0.9)
        .expect("rotation")
}

fn sphere(center: Point3, radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center,
        radius,
        axis: Vec3::Z,
    }
}

/// Solve upright and turned; check both and that they agree.
fn check(a: &SurfaceGeometry, b: &SurfaceGeometry, start: Point3, end: Point3, what: &str) {
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let upright = compute_edge_geometry(a, b, EdgeAnchors::Open { start, end }, None, &tol, &opts)
        .unwrap_or_else(|e| panic!("{what} upright: {e:?}"));
    let t = turn();
    let (ta, tb) = (
        crate::topology::brep::transform_surface(a, &t, 1.0),
        crate::topology::brep::transform_surface(b, &t, 1.0),
    );
    let turned_pole = [&ta, &tb]
        .iter()
        .any(|s| matches!(s, SurfaceGeometry::Sphere { axis, .. } if (axis.z - 1.0).abs() > 1e-6));
    assert!(turned_pole, "{what}: the turn moves the sphere's pole");
    let (ts, te) = (t.apply_point(&start), t.apply_point(&end));
    let turned = compute_edge_geometry(
        &ta,
        &tb,
        EdgeAnchors::Open { start: ts, end: te },
        None,
        &tol,
        &opts,
    )
    .unwrap_or_else(|e| panic!("{what} turned: {e:?}"));
    for triple in [&upright, &turned] {
        assert!(
            matches!(
                triple.provenance,
                EdgeTripleProvenance::AnalyticIntersection { .. }
            ),
            "{what}: analytic dispatch, got {:?}",
            triple.provenance
        );
        assert_eq!(triple.parameter_map, ParameterMap::Identity, "{what}");
    }
    assert_guarantees(&upright, a, b, start, end, &tol);
    assert_guarantees(&turned, &ta, &tb, ts, te, &tol);
    // The turned curve is the rotated upright one (same shape, same
    // sense): compare points at matching fractions of each range.
    for k in 0..=8 {
        let f = f64::from(k) / 8.0;
        let at = |tr: &crate::geometry::edge_triple::EdgeTriple| {
            tr.curve_3d
                .point_at(tr.t_start + (tr.t_end - tr.t_start) * f)
        };
        let d = t.apply_point(&at(&upright)).distance_to(&at(&turned));
        assert!(
            d < 1e-9,
            "{what}: turned curve off the rotated upright by {d:e} at {f}"
        );
    }
}

#[test]
fn turned_sphere_pairs_serve_the_analytic_dispatch() {
    let (r3, s15) = (3.0_f64.sqrt(), 3.75_f64.sqrt());
    // Plane z = 1 through a sphere of radius 2: the circle r = √3.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::Z,
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = sphere(Point3::new(0.0, 0.0, 0.0), 2.0);
    check(
        &plane,
        &s2,
        Point3::new(r3, 0.0, 1.0),
        Point3::new(0.0, r3, 1.0),
        "plane x sphere",
    );
    // A unit cylinder on the sphere's axis: circles at z = ±√3.
    check(
        &z_axis_cyl(1.0),
        &s2,
        Point3::new(1.0, 0.0, r3),
        Point3::new(0.0, 1.0, r3),
        "cylinder x sphere",
    );
    // Two spheres on one pole: the circle z = 0.5, r = √3.75.
    check(
        &s2,
        &sphere(Point3::new(0.0, 0.0, 1.0), 2.0),
        Point3::new(s15, 0.0, 0.5),
        Point3::new(0.0, s15, 0.5),
        "sphere x sphere",
    );
    // A sphere off a cone's axis: the quartic of the quadric-pair suite.
    let v_at = |u: f64| (3.0_f64.sqrt() / 4.0) * (u.cos() + (u.cos() * u.cos() + 4.0).sqrt());
    let pt_at = |u: f64| {
        let v = v_at(u);
        let r = v / 3.0_f64.sqrt();
        Point3::new(r * u.cos(), r * u.sin(), v)
    };
    check(
        &sphere(Point3::new(1.0, 0.0, 0.0), 2.0_f64.sqrt()),
        &z_axis_cone(std::f64::consts::FRAC_PI_6),
        pt_at(0.0),
        pt_at(std::f64::consts::FRAC_PI_3),
        "sphere x cone",
    );
}
