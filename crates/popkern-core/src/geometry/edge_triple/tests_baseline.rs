use super::tests_support::*;
use crate::geometry::Point3;
use crate::geometry::edge_triple::{
    AnalyticPair, BranchPolicy, BranchSelector, CancellationToken, EdgeAnchors, EdgeRequest,
    EdgeTripleBatchOptions, EdgeTripleError, EdgeTripleOptions, EdgeTripleProvenance,
    compute_edge_geometry, compute_edge_geometry_batch,
};
use crate::tolerance::Tolerance;

#[test]
fn unimplemented_for_unsupported_pair() {
    // Stage 2b.52: the `Unimplemented` catch-all is only reached
    // when the trig, analytic, AND marching dispatchers all return
    // None. Marching now covers torus and NURBS pairs, so two
    // identical tori (used to be `Unimplemented`) now route through
    // it. To exercise the `Unimplemented` arm we need a pair the
    // dispatch tree truly doesn't know — currently impossible to
    // construct from `SurfaceGeometry` variants since marching is a
    // catch-all. We rename the test to assert the new contract:
    // any well-formed surface pair gets *some* outcome other than
    // `Unimplemented`. The test passes if the result is either a
    // success or any specific error other than `Unimplemented`.
    let a = nurbs_torus_dummy();
    let b = nurbs_torus_dummy();
    let tol = Tolerance::default();
    let opts = EdgeTripleOptions::default();
    let anchors = EdgeAnchors::Open {
        start: Point3::new(0.0, 0.0, 0.0),
        end: Point3::new(1.0, 0.0, 0.0),
    };
    let result = compute_edge_geometry(&a, &b, anchors, None, &tol, &opts);
    assert!(
        !matches!(result, Err(EdgeTripleError::Unimplemented)),
        "marching catch-all should preempt Unimplemented; got {result:?}"
    );
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
    // Stage 2b.52: marching catch-all preempts Unimplemented. We
    // just verify the batch returns the same outcome for both
    // identical inputs (determinism).
    let (Some(first), Some(second)) = (results.first(), results.get(1)) else {
        panic!("batch returned wrong number of results");
    };
    match (first, second) {
        (Err(_), Err(_)) | (Ok(_), Ok(_)) => {}
        _ => panic!("batch determinism violated: {:?} vs {:?}", first, second),
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
        cardinality: None,
    };
    let s = serde_json::to_string(&p).unwrap();
    let p2: EdgeTripleProvenance = serde_json::from_str(&s).unwrap();
    assert_eq!(p, p2);
}

// ─── Stage 2a: per-pair analytic dispatch ─────────────────────────
