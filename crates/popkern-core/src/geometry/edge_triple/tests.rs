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
fn plane_cone_returns_unimplemented_until_stage_2b() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let cone = z_axis_cone(std::f64::consts::FRAC_PI_4);
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, 1.0),
        end: Point3::new(-1.0, 0.0, 1.0),
    };
    let result = compute_edge_geometry(&plane, &cone, anchors, None, &tol, &opts);
    assert!(matches!(result, Err(EdgeTripleError::Unimplemented)));
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
