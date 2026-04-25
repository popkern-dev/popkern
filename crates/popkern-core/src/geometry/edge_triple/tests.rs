//! Stage 1 tests: API skeleton compiles, every dispatch path
//! returns `Unimplemented`. Stages 2+ replace these with the real
//! per-pair-class invariant tests.

use super::*;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

fn dummy_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

#[test]
fn compute_returns_unimplemented_for_every_pair() {
    let surf = dummy_plane();
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let anchors = EdgeAnchors::Open {
        start: Point3::new(0.0, 0.0, 0.0),
        end: Point3::new(1.0, 0.0, 0.0),
    };
    let result = compute_edge_geometry(&surf, &surf, anchors, None, &tol, &opts);
    assert!(matches!(result, Err(EdgeTripleError::Unimplemented)));
}

#[test]
fn batch_returns_unimplemented_for_every_request() {
    let surf = dummy_plane();
    let tol = Tolerance::default();
    let opts = EdgeTripleBatchOptions::default();
    let req = EdgeRequest {
        surf_a: &surf,
        surf_b: &surf,
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
fn curvature_comb_stub_returns_empty() {
    let triple_json = r#"{
        "curve_3d":{"control_points":[],"weights":[],"knots":[],"degree":3},
        "pcurve_a":{"control_points":[],"weights":[],"knots":[],"degree":3},
        "pcurve_b":{"control_points":[],"weights":[],"knots":[],"degree":3},
        "t_start":0.0,"t_end":1.0,
        "topology":"Open",
        "provenance":{"AnalyticIntersection":{"kind":"PlanePlane","branch":"OnlyBranch"}},
        "quality":{"max_positional_residual":0.0,"max_endpoint_error":0.0,"samples":0,"tolerance_used":0.0}
    }"#;
    // Just check serialisation round-trips without exercising the
    // (empty) curvature comb implementation. The triple_json is
    // illustrative only; the curvature comb is a stub at stage 1.
    let _ = triple_json;
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
