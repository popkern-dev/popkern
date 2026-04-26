//! Tests for the EdgeTriple factory.
//!
//! Stage 1 baseline checks (cancellation, options, batch, provenance
//! serde) plus stage 2a analytic dispatch tests for plane-pair
//! specialisations (plane×plane, plane×cylinder, plane×sphere,
//! plane×cone).
//!
//! Each pair-class test asserts the hard guarantees G1–G6 at
//! `tol.linear()`.

use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

// ─── Fixtures ─────────────────────────────────────────────────────

fn xy_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn xz_plane_at(z0: f64) -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, z0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn z_axis_cyl(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

fn unit_sphere_at(c: Point3) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: c,
        radius: 1.0,
    }
}

fn z_axis_cone(half_angle: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle,
    }
}

fn nurbs_torus_dummy() -> SurfaceGeometry {
    SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    }
}

// ─── Hard-guarantee assertions ────────────────────────────────────

fn assert_guarantees(
    triple: &EdgeTriple,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors_start: Point3,
    anchors_end: Point3,
    tol: &Tolerance,
) {
    // G3 / G5: t_start ≤ t_end and curve.point_at(t_start) is the
    // start anchor.
    assert!(
        triple.t_start <= triple.t_end,
        "G3/G5: t_start={} t_end={}",
        triple.t_start,
        triple.t_end
    );
    let p_start = triple.curve_3d.point_at(triple.t_start);
    let p_end = triple.curve_3d.point_at(triple.t_end);
    assert!(
        p_start.distance_to(&anchors_start) <= tol.linear(),
        "G5 start: {:?} vs {:?}",
        p_start,
        anchors_start
    );
    assert!(
        p_end.distance_to(&anchors_end) <= tol.linear(),
        "G5 end: {:?} vs {:?}",
        p_end,
        anchors_end
    );

    // G1: endpoint CP equality. For NURBS variants, check first/last
    // CPs explicitly. For analytic primitives the endpoint comes from
    // the surface evaluation and is exact-by-construction; the
    // distance check above already covers it.
    if let crate::geometry::curve::CurveGeometry::Nurbs(n) = &triple.curve_3d {
        let cps = &n.control_points;
        assert_eq!(*cps.first().unwrap(), anchors_start, "G1 first CP");
        assert_eq!(*cps.last().unwrap(), anchors_end, "G1 last CP");
    }

    // G2: per-sample positional invariant. Verify already runs 32
    // samples — call again to assert the achieved residual fits.
    let q = triple.verify(surf_a, surf_b, tol, 32).expect("G2 verify");
    assert!(
        q.max_positional_residual <= tol.linear(),
        "G2 residual {:?}",
        q
    );

    // G4: UV-box containment for bounded surfaces is checked
    // inside verify(); if we got here without an error, we're fine.

    // G6: determinism — second invocation of verify is bit-identical.
    let q2 = triple.verify(surf_a, surf_b, tol, 32).expect("G6 verify");
    assert_eq!(
        q.max_positional_residual.to_bits(),
        q2.max_positional_residual.to_bits()
    );
}

// ─── Stage 1: still-passing baseline tests ────────────────────────

#[test]
fn unimplemented_for_unsupported_pair() {
    let a = nurbs_torus_dummy();
    let b = nurbs_torus_dummy();
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let anchors = EdgeAnchors::Open {
        start: Point3::new(0.0, 0.0, 0.0),
        end: Point3::new(1.0, 0.0, 0.0),
    };
    let result = compute_edge_geometry(&a, &b, anchors, None, &tol, &opts);
    assert!(matches!(result, Err(EdgeTripleError::Unimplemented)));
}

#[test]
fn batch_passes_through_dispatch() {
    let a = nurbs_torus_dummy();
    let b = nurbs_torus_dummy();
    let tol = Tolerance::default();
    let opts = EdgeTripleBatchOptions::default();
    let req = EdgeRequest {
        surf_a: &a,
        surf_b: &b,
        anchors: EdgeAnchors::Open {
            start: Point3::new(0.0, 0.0, 0.0),
            end: Point3::new(1.0, 0.0, 0.0),
        },
        hint: None,
    };
    let results = compute_edge_geometry_batch(&[req.clone(), req], &tol, &opts);
    assert_eq!(results.len(), 2);
    for r in results {
        assert!(matches!(r, Err(EdgeTripleError::Unimplemented)));
    }
}

#[test]
fn cancellation_token_round_trip() {
    let token = CancellationToken::new();
    assert!(!token.is_cancelled());
    token.cancel();
    assert!(token.is_cancelled());
    let cloned = token.clone();
    assert!(cloned.is_cancelled());
}

#[test]
fn options_default_values() {
    let opts = EdgeTripleOptions::default();
    assert_eq!(opts.max_marching_iterations, 250);
    assert_eq!(opts.prefer_branches, BranchPolicy::ClosestToAnchors);
    assert!(!opts.verify_in_release);
    assert!(opts.cancel.is_none());
    assert!(opts.progress.is_none());
}

#[test]
fn provenance_serialises_round_trip() {
    use serde_json;
    let p = EdgeTripleProvenance::AnalyticIntersection {
        kind: AnalyticPair::CylinderTorus,
        branch: BranchSelector::Indexed(2),
    };
    let s = serde_json::to_string(&p).unwrap();
    let p2: EdgeTripleProvenance = serde_json::from_str(&s).unwrap();
    assert_eq!(p, p2);
}

// ─── Stage 2a: per-pair analytic dispatch ─────────────────────────

#[test]
fn plane_plane_intersection_line() {
    let a = xy_plane();
    let b = xz_plane_at(0.0); // y = 0 plane through origin
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let start = Point3::new(-2.0, 0.0, 0.0);
    let end = Point3::new(2.0, 0.0, 0.0);
    let anchors = EdgeAnchors::Open { start, end };

    let triple = compute_edge_geometry(&a, &b, anchors, None, &tol, &opts)
        .expect("plane×plane analytic dispatch");

    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlanePlane,
            ..
        }
    ));
    assert_guarantees(&triple, &a, &b, start, end, &tol);
}

#[test]
fn plane_cylinder_perpendicular_analytic_circle() {
    // Stage 2b.3: plane perpendicular to cyl axis returns an
    // analytic Circle 3D + Linear cyl pcurve + Elliptical plane
    // pcurve, with `ParameterMap::Identity`.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.5),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let cyl = z_axis_cyl(1.0);
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let start = Point3::new(1.0, 0.0, 0.5);
    let end = Point3::new(-1.0, 0.0, 0.5);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&plane, &cyl, anchors, None, &tol, &opts)
        .expect("plane×cyl analytic dispatch");

    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneCylinder,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &plane, &cyl, start, end, &tol);
}

#[test]
fn plane_cylinder_oblique_analytic_ellipse() {
    // Stage 2b.11: tilt the cutting plane 30° off horizontal. The
    // intersection is a 3-D ellipse; the dispatcher must return
    // `CurveGeometry::Ellipse` with the cyl side as a Sinusoidal
    // pcurve and the plane side as an Elliptical pcurve.
    let cyl = z_axis_cyl(1.0);
    // Plane normal tilted 30° from +Z toward +X: passes through the
    // origin, normal = (sin 30°, 0, cos 30°) = (0.5, 0, √3/2).
    let n = Vec3::new(0.5, 0.0, 3.0_f64.sqrt() / 2.0);
    // Choose plane u in the cutting plane, perpendicular to the
    // global Y axis projection. We pick u = (cos 30°, 0, -sin 30°)
    // and v = (0, 1, 0). Both are perpendicular to n.
    let pu = Vec3::new(3.0_f64.sqrt() / 2.0, 0.0, -0.5);
    let pv = Vec3::new(0.0, 1.0, 0.0);
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: n,
        u_axis: pu,
        v_axis: pv,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();

    // Anchor at θ = 0 on the cylinder: x = 1, y = 0; the cutting
    // plane through origin with normal n places z so that
    // (x, y, z)·n = 0  ⇒  z = -x · (sin30°/cos30°) = -1/√3.
    let z_for = |x: f64, y: f64| -(x * n.x + y * n.y) / n.z;
    let start = Point3::new(1.0, 0.0, z_for(1.0, 0.0));
    let end = Point3::new(-1.0, 0.0, z_for(-1.0, 0.0));
    let anchors = EdgeAnchors::Open { start, end };

    let triple = compute_edge_geometry(&plane, &cyl, anchors, None, &tol, &opts)
        .expect("plane×cyl oblique analytic dispatch");

    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneCylinder,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Ellipse { .. }
    ));
    // The cyl-side pcurve must be Sinusoidal (analytic, not NURBS).
    assert!(
        triple.pcurve_b.is_analytic(),
        "cyl pcurve should be analytic"
    );
    assert!(matches!(
        triple.pcurve_b.geometry,
        crate::geometry::pcurve::PcurveGeometry::Sinusoidal { .. }
    ));
    assert_guarantees(&triple, &plane, &cyl, start, end, &tol);
}

#[test]
fn plane_cone_oblique_analytic_ellipse() {
    // Stage 2b.12: cone half-angle 30°, plane tilted 30° from
    // perpendicular. With angle α = 30° from axis and cone θ = 30°,
    // |cos α| = √3/2 > tan(30°)·sin(30°) = 1/(2√3) — ellipse case.
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_6); // 30°
    // Plane passes through (0, 0, 1) with normal tilted 30° from
    // +Z toward +X: (sin 30°, 0, cos 30°) = (0.5, 0, √3/2).
    let n = Vec3::new(0.5, 0.0, 3.0_f64.sqrt() / 2.0);
    let pu = Vec3::new(3.0_f64.sqrt() / 2.0, 0.0, -0.5);
    let pv = Vec3::new(0.0, 1.0, 0.0);
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: n,
        u_axis: pu,
        v_axis: pv,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();

    // Anchors at cone u=0 and u=π. Solve cone v from plane equation:
    //   point = (±v/√3, 0, v); plane: 0.5·x + √3/2·(z-1) = 0
    //   ⇒ v(0)   = 3/4         (u=0):  (√3/4, 0, 3/4)
    //   ⇒ v(π)   = 3/2         (u=π):  (−√3/2, 0, 3/2)
    let start = Point3::new(3.0_f64.sqrt() / 4.0, 0.0, 0.75);
    let end = Point3::new(-3.0_f64.sqrt() / 2.0, 0.0, 1.5);
    let anchors = EdgeAnchors::Open { start, end };

    let triple = compute_edge_geometry(&plane, &cone, anchors, None, &tol, &opts)
        .expect("plane×cone oblique analytic dispatch");

    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Ellipse { .. }
    ));
    assert_guarantees(&triple, &plane, &cone, start, end, &tol);
}

#[test]
fn plane_cone_oblique_analytic_parabola() {
    // Stage 2b.13: cone half-angle 30°, plane normal at exactly
    // 60° from the axis — the case-boundary between ellipse and
    // hyperbola. Result: an analytic parabola.
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_6);
    // n = (sin 60°, 0, cos 60°) = (√3/2, 0, 1/2). pu = (cos 60°, 0,
    // -sin 60°) = (1/2, 0, -√3/2). pv = (0, 1, 0). Plane through
    // (1, 0, 0): cuts the cone along a parabola whose vertex is at
    // (1/2, 0, √3/2) — the cone point at u = 0.
    let n = Vec3::new(3.0_f64.sqrt() / 2.0, 0.0, 0.5);
    let pu = Vec3::new(0.5, 0.0, -3.0_f64.sqrt() / 2.0);
    let pv = Vec3::new(0.0, 1.0, 0.0);
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 0.0, 0.0),
        normal: n,
        u_axis: pu,
        v_axis: pv,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();

    // Anchors at cone u = ±π/3, derived analytically from the plane
    // equation v(u) = √3 / (1 + cos u). At u = π/3, v = 2√3/3, so
    // the cone point is (1/3, √3/3, 2√3/3).
    let v_at_third = 2.0 * 3.0_f64.sqrt() / 3.0;
    let start = Point3::new(1.0 / 3.0, -3.0_f64.sqrt() / 3.0, v_at_third);
    let end = Point3::new(1.0 / 3.0, 3.0_f64.sqrt() / 3.0, v_at_third);
    let anchors = EdgeAnchors::Open { start, end };

    let triple = compute_edge_geometry(&plane, &cone, anchors, None, &tol, &opts)
        .expect("plane×cone parabola analytic dispatch");

    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Parabola { .. }
    ));
    assert_guarantees(&triple, &plane, &cone, start, end, &tol);
}

#[test]
fn plane_cone_oblique_analytic_hyperbola() {
    // Stage 2b.14: cone half-angle 30°, plane normal at 75° from
    // axis (steeper than the case-boundary at 60°). Result:
    // analytic hyperbola.
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_6);
    let n_angle = 75.0_f64.to_radians();
    let n = Vec3::new(n_angle.sin(), 0.0, n_angle.cos());
    let pu = Vec3::new(n_angle.cos(), 0.0, -n_angle.sin());
    let pv = Vec3::new(0.0, 1.0, 0.0);
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 0.0, 0.0),
        normal: n,
        u_axis: pu,
        v_axis: pv,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();

    // Anchors at cone u = ±π/3. Plane equation: 0.9659·x +
    // 0.2588·z = 0.9659 ⇒ z = (0.9659 − 0.9659·x)/0.2588 (when
    // y = 0). For points on cone (u, v): x = v·cos(u)/√3, y =
    // v·sin(u)/√3, z = v. Combine to get v as a function of u.
    let cn = n.dot(&Vec3::new(0.0, 0.0, 1.0));
    let p1 = Vec3::new(1.0, 0.0, 0.0).dot(&n);
    let p2 = Vec3::new(0.0, 1.0, 0.0).dot(&n);
    let t = std::f64::consts::FRAC_PI_6.tan();
    let d_signed = (Point3::new(1.0, 0.0, 0.0) - Point3::new(0.0, 0.0, 0.0)).dot(&n);

    let cone_pt = |u: f64| -> Point3 {
        let v = d_signed / (cn + t * (p1 * u.cos() + p2 * u.sin()));
        Point3::new(
            v * u.cos() / 3.0_f64.sqrt(),
            v * u.sin() / 3.0_f64.sqrt(),
            v,
        )
    };
    let start = cone_pt(-std::f64::consts::FRAC_PI_3);
    let end = cone_pt(std::f64::consts::FRAC_PI_3);
    let anchors = EdgeAnchors::Open { start, end };

    let triple = compute_edge_geometry(&plane, &cone, anchors, None, &tol, &opts)
        .expect("plane×cone hyperbola analytic dispatch");

    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Hyperbola { .. }
    ));
    assert_guarantees(&triple, &plane, &cone, start, end, &tol);
}

#[test]
fn plane_sphere_oblique_analytic_circle() {
    // Stage 2b.15: plane normal NOT aligned with the sphere's polar
    // axis. The polar dispatcher returns None and we fall through
    // to the general dispatcher, which still produces an exact
    // analytic Arc (every plane-sphere intersection is a circle).
    let sphere = unit_sphere_at(Point3::new(0.0, 0.0, 0.0));
    // Plane through origin with normal (1, 0, 1)/√2 — a non-axis-
    // aligned great-circle slice that doesn't pass through either
    // sphere pole, so the (longitude, colatitude) NURBS fit on the
    // sphere stays well-conditioned.
    let inv_root2 = 1.0 / 2.0_f64.sqrt();
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(inv_root2, 0.0, inv_root2),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(-inv_root2, 0.0, inv_root2),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Two points on the great circle (orthogonal to plane normal,
    // unit distance from sphere centre); both stay in the equatorial
    // band on the sphere, away from poles.
    let start = Point3::new(0.0, 1.0, 0.0);
    let end = Point3::new(-inv_root2, 0.0, inv_root2);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&plane, &sphere, anchors, None, &tol, &opts)
        .expect("plane×sphere general analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneSphere,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &plane, &sphere, start, end, &tol);
}

#[test]
fn sphere_sphere_oblique_analytic_circle() {
    // Stage 2b.20: two unit spheres with centres NOT aligned on the
    // polar axis. Polar dispatcher returns None, general dispatcher
    // produces an exact analytic Arc with NURBS pcurves on both
    // spheres. Periodicity-aware check_uv_box (added to the surface
    // API in this stage) absorbs minor LSQ-fit overshoots in the
    // longitude coord.
    let sphere_a = unit_sphere_at(Point3::new(0.0, 0.0, 0.0));
    let sphere_b = unit_sphere_at(Point3::new(1.0, 1.0, 0.0));
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let inv_root2 = 1.0 / 2.0_f64.sqrt();
    let start = Point3::new(0.5, 0.5, inv_root2);
    let end = Point3::new(0.5, 0.5, -inv_root2);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&sphere_a, &sphere_b, anchors, None, &tol, &opts)
        .expect("sphere×sphere general analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereSphere,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &sphere_a, &sphere_b, start, end, &tol);
}

#[test]
fn plane_sphere_polar_analytic_circle() {
    // Stage 2b.4: equator slice of a unit sphere by the XY plane.
    // Expected: analytic Arc 3D + Linear sphere pcurve at v=π/2 +
    // Elliptical plane pcurve, ParameterMap::Identity.
    let plane = xy_plane();
    let sphere = unit_sphere_at(Point3::new(0.0, 0.0, 0.0));
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let start = Point3::new(1.0, 0.0, 0.0);
    let end = Point3::new(-1.0, 0.0, 0.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&plane, &sphere, anchors, None, &tol, &opts)
        .expect("plane×sphere analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneSphere,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert_guarantees(&triple, &plane, &sphere, start, end, &tol);
}

#[test]
fn plane_torus_perpendicular_analytic_circle() {
    // Stage 2b.16: cutting plane perpendicular to the torus axis at
    // height h = 0.5·minor_radius. Two intersection circles appear;
    // the dispatcher selects whichever branch the anchors lie on.
    let torus = nurbs_torus_dummy(); // R = 5, r = 1
    // Plane at z = 0.5 (h = 0.5 · r), perpendicular to axis (0,0,1).
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.5),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // sin v = 0.5 ⇒ v_outer = π/6 (cos v = √3/2). Outer circle radius
    // = R + r·cos v = 5 + √3/2 ≈ 5.866.
    let r_circ = 5.0 + 3.0_f64.sqrt() / 2.0;
    // Anchors at u = 0 and u = π/2 on the outer circle.
    let start = Point3::new(r_circ, 0.0, 0.5);
    let end = Point3::new(0.0, r_circ, 0.5);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&plane, &torus, anchors, None, &tol, &opts)
        .expect("plane×torus axis-perpendicular analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneTorus,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &plane, &torus, start, end, &tol);
}

#[test]
fn sphere_cone_axial_analytic_circle() {
    // Stage 2b.18: unit sphere centred at apex of cone (axis (0,0,1),
    // half-angle 30°). Solving v²·sec²θ − 2dv + (d² − 1) = 0 with
    // d = 0 gives v = cos(30°) = √3/2; intersection circle radius
    // = v·tan(30°) = 1/2.
    let sphere = unit_sphere_at(Point3::new(0.0, 0.0, 0.0));
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_6);
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let v = 3.0_f64.sqrt() / 2.0;
    let r_circ = 0.5;
    let start = Point3::new(r_circ, 0.0, v);
    let end = Point3::new(0.0, r_circ, v);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&sphere, &cone, anchors, None, &tol, &opts)
        .expect("sphere×cone axial analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &sphere, &cone, start, end, &tol);
}

#[test]
fn cone_cone_coaxial_analytic_circle() {
    // Stage 2b.19: cone A apex at origin, half-angle 30°; cone B
    // apex at (0, 0, 1), half-angle 60°. Solving v_a·tan30 =
    // (v_a−1)·tan60 gives v_a · (1/√3 − √3) = -√3, so
    // v_a = √3/(√3 − 1/√3) = √3/(2/√3) = 3/2. circle_radius =
    // (3/2)·tan(30°) = √3/2. Centre at (0, 0, 1.5).
    let cone_a = z_axis_cone(std::f64::consts::FRAC_PI_6); // apex (0,0,0), 30°
    let cone_b = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 1.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_3, // 60°
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let r_circ = 3.0_f64.sqrt() / 2.0;
    let start = Point3::new(r_circ, 0.0, 1.5);
    let end = Point3::new(0.0, r_circ, 1.5);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cone_a, &cone_b, anchors, None, &tol, &opts)
        .expect("cone×cone coaxial analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::ConeCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &cone_a, &cone_b, start, end, &tol);
}

#[test]
fn cyl_cone_coaxial_analytic_circle() {
    // Stage 2b.17: cylinder of radius 1 along z-axis, cone with apex
    // at origin, axis (0,0,1), half-angle 30°. The intersection is a
    // single circle at height v_cone = 1 / tan(30°) = √3, radius 1.
    let cyl = z_axis_cyl(1.0);
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_6);
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let h = 3.0_f64.sqrt();
    let start = Point3::new(1.0, 0.0, h);
    let end = Point3::new(0.0, 1.0, h);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl, &cone, anchors, None, &tol, &opts)
        .expect("cyl×cone coaxial analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    assert_guarantees(&triple, &cyl, &cone, start, end, &tol);
}

#[test]
fn plane_torus_villarceau_analytic_circle() {
    // Stage 2b.22: Villarceau plane bitangent to a torus, factoring
    // the quartic into 2 circles of major radius. Torus R=5, r=1.
    // Villarceau angle α = arcsin(1/5) ≈ 0.2014 rad.
    // Plane normal: (sin α, 0, cos α) = (1/5, 0, √24/5).
    let torus = nurbs_torus_dummy(); // R=5, r=1
    let s_a = 0.2_f64; // sin α = 1/5
    let c_a = (1.0_f64 - s_a * s_a).sqrt(); // cos α = √24/5
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(s_a, 0.0, c_a),
        u_axis: Vec3::new(c_a, 0.0, -s_a),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Anchors on the +y Villarceau circle (centre (0, 1, 0), R=5).
    // At θ=0 in arc-basis frame: (R, 0, 0) projected — but the arc
    // basis depends on plane normal. Just pick two known points on
    // the +y circle: outer-equator point (0, 6, 0) (= (0, 1, 0) +
    // 5·(0, 1, 0))? Wait, circle centre is (0, 1, 0) and radius 5,
    // lies in the Villarceau plane. So circle points satisfy
    // x·sin α + z·cos α = 0. At (0, 6, 0): 0+0=0 ✓ and distance
    // from centre = (0, 5, 0) magnitude 5 ✓.
    let start = Point3::new(0.0, 6.0, 0.0);
    // Other known point: where x·s_a + z·c_a = 0, distance 5 from
    // (0, 1, 0). Try x = 5·c_a, z = -5·s_a, y = 1: distance from
    // centre = (5·c_a, 0, -5·s_a), magnitude = 5·sqrt(c_a²+s_a²) = 5 ✓.
    let end = Point3::new(5.0 * c_a, 1.0, -5.0 * s_a);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&plane, &torus, anchors, None, &tol, &opts)
        .expect("plane×torus Villarceau analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneTorus,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Arc { .. }
    ));
    // Major radius circle.
    if let crate::geometry::curve::CurveGeometry::Arc { radius, .. } = triple.curve_3d {
        assert!(
            (radius - 5.0).abs() < 1e-9,
            "Villarceau circle should have radius = R = 5, got {radius}"
        );
    }
    assert_guarantees(&triple, &plane, &torus, start, end, &tol);
}

#[test]
fn sphere_cone_non_axial_quartic_curve() {
    // Stage 2b.24: sphere centred OFF the cone axis. The intersection
    // is a genuine 4th-degree curve. The dispatcher samples via the
    // closed-form quadratic in v at each cone-u angle and emits a
    // NURBS curve_3d.
    //
    // Cone: apex (0,0,0), axis (0,0,1), half-angle 30°.
    // Sphere: centre (1, 0, 0), radius √2. Solving:
    //   |(v/√3·cos u − 1, v/√3·sin u, v) − (0,0,0)|² = 2
    //   ⇒ 4v²/3 − (2v/√3)·cos u − 1 = 0
    //   ⇒ v = (√3/4)·(cos u + √(cos²u + 4))     (positive branch)
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_6);
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 0.0, 0.0),
        radius: 2.0_f64.sqrt(),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let v_at = |u: f64| -> f64 {
        let cu = u.cos();
        (3.0_f64.sqrt() / 4.0) * (cu + (cu * cu + 4.0).sqrt())
    };
    let pt_at = |u: f64| -> Point3 {
        let v = v_at(u);
        let r = v / 3.0_f64.sqrt();
        Point3::new(r * u.cos(), r * u.sin(), v)
    };
    let start = pt_at(0.0);
    let end = pt_at(std::f64::consts::FRAC_PI_3);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&sphere, &cone, anchors, None, &tol, &opts)
        .expect("sphere×cone non-axial dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    // Genuinely quartic ⇒ NURBS curve_3d.
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Nurbs(_)
    ));
    assert_guarantees(&triple, &sphere, &cone, start, end, &tol);
}

#[test]
fn cyl_cyl_perpendicular_viviani_curve() {
    // Stage 2b.23: two unit cylinders with perpendicular intersecting
    // axes — the classical Viviani's curve. The intersection is
    // genuinely 4th-order and has no analytic Line/Arc/Ellipse/Para/
    // Hyper representation; the dispatcher emits a NURBS curve_3d
    // sampled along the implicit-pair solution.
    let cyl_a = z_axis_cyl(1.0);
    let cyl_b = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 1.0,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Anchors on the upper-z branch: at u_cyl_a = 0 and u_cyl_a = π/4.
    // x = cos u, y = sin u, z = +sqrt(1 − sin²u) = +cos u (Viviani:
    // z = ±x for equal radii). At u=0: (1, 0, 1). At u=π/4:
    // (√2/2, √2/2, +sqrt(1 − 1/2)) = (√2/2, √2/2, √2/2).
    let r2 = 2.0_f64.sqrt() / 2.0;
    let start = Point3::new(1.0, 0.0, 1.0);
    let end = Point3::new(r2, r2, r2);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl_a, &cyl_b, anchors, None, &tol, &opts)
        .expect("cyl×cyl perpendicular dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderCylinder,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    // 4th-order quartic ⇒ NURBS curve_3d (no analytic primitive).
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Nurbs(_)
    ));
    assert_guarantees(&triple, &cyl_a, &cyl_b, start, end, &tol);
}

#[test]
fn cyl_cyl_parallel_analytic_line() {
    // Two unit cylinders along z-axis, centres at (0,0,0) and (1,0,0)
    // (axis offset 1.0 < r1+r2 = 2.0). Two intersection lines at
    // x = 0.5, y = ±√(1 - 0.25) = ±√0.75.
    let cyl_a = z_axis_cyl(1.0);
    let cyl_b = SurfaceGeometry::Cylinder {
        origin: Point3::new(1.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let h = (1.0_f64 - 0.25).sqrt();
    let start = Point3::new(0.5, h, 0.0);
    let end = Point3::new(0.5, h, 2.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl_a, &cyl_b, anchors, None, &tol, &opts)
        .expect("cyl×cyl parallel analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderCylinder,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert!(matches!(
        triple.curve_3d,
        crate::geometry::curve::CurveGeometry::Line { .. }
    ));
    assert_guarantees(&triple, &cyl_a, &cyl_b, start, end, &tol);
}

#[test]
fn cyl_sphere_axial_analytic_circle() {
    // Z-axis cylinder of radius 1; sphere at origin, radius √2.
    // Intersection: two circles at z = ±1, each of radius 1.
    let cyl = z_axis_cyl(1.0);
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0_f64.sqrt(),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    // Anchors on the z=+1 circle.
    let start = Point3::new(1.0, 0.0, 1.0);
    let end = Point3::new(-1.0, 0.0, 1.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&cyl, &sphere, anchors, None, &tol, &opts)
        .expect("cyl×sphere axial analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::CylinderSphere,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert_guarantees(&triple, &cyl, &sphere, start, end, &tol);
}

#[test]
fn sphere_sphere_polar_analytic_circle() {
    // Two unit spheres at z=0 and z=1, centres along the polar axis.
    // Intersection: a circle at z=0.5 with radius √(1 - 0.25) ≈ 0.866.
    let a = unit_sphere_at(Point3::new(0.0, 0.0, 0.0));
    let b = unit_sphere_at(Point3::new(0.0, 0.0, 1.0));
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let r = (1.0_f64 - 0.25).sqrt();
    let start = Point3::new(r, 0.0, 0.5);
    let end = Point3::new(-r, 0.0, 0.5);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&a, &b, anchors, None, &tol, &opts)
        .expect("sphere×sphere analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereSphere,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert_guarantees(&triple, &a, &b, start, end, &tol);
}

#[test]
fn plane_cone_perpendicular_analytic_circle() {
    // Stage 2b.5: 45° cone, plane perpendicular to axis at z=1.
    // Circle radius = h * tan(half_angle) = 1.0 * 1.0 = 1.0.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_4);
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let start = Point3::new(1.0, 0.0, 1.0);
    let end = Point3::new(-1.0, 0.0, 1.0);
    let anchors = EdgeAnchors::Open { start, end };
    let triple = compute_edge_geometry(&plane, &cone, anchors, None, &tol, &opts)
        .expect("plane×cone analytic dispatch");
    assert!(matches!(
        triple.provenance,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneCone,
            ..
        }
    ));
    assert_eq!(triple.parameter_map, ParameterMap::Identity);
    assert_guarantees(&triple, &plane, &cone, start, end, &tol);
}

#[test]
fn analytic_dispatch_returns_no_intersection_for_parallel_planes() {
    let a = xy_plane();
    let b = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let anchors = EdgeAnchors::Open {
        start: Point3::new(0.0, 0.0, 0.0),
        end: Point3::new(1.0, 0.0, 0.0),
    };
    let result = compute_edge_geometry(&a, &b, anchors, None, &tol, &opts);
    assert!(matches!(result, Err(EdgeTripleError::NoIntersection)));
}

#[test]
fn determinism_identical_inputs_identical_output() {
    // G6: same inputs → bit-identical output. Exercised on
    // plane×plane (the only stage-2a working pair).
    let a = xy_plane();
    let b = xz_plane_at(0.0);
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let anchors = EdgeAnchors::Open {
        start: Point3::new(-2.0, 0.0, 0.0),
        end: Point3::new(2.0, 0.0, 0.0),
    };

    let t1 = compute_edge_geometry(&a, &b, anchors.clone(), None, &tol, &opts).unwrap();
    let t2 = compute_edge_geometry(&a, &b, anchors, None, &tol, &opts).unwrap();

    let s1 = serde_json::to_string(&t1).unwrap();
    let s2 = serde_json::to_string(&t2).unwrap();
    assert_eq!(
        s1, s2,
        "G6 determinism: identical inputs → identical output"
    );
}
