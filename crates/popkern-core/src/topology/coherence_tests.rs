use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;

#[test]
fn new_brep_has_empty_pending_coherence() {
    let brep = BRep::new(Tolerance::default());
    assert!(brep.pending_coherence.is_empty());
}

#[test]
fn record_pending_authoring_stores_marker() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
        edge: mev.edge,
        owner_operator: "test_op",
    });

    assert_eq!(brep.pending_coherence.len(), 1);
}

#[test]
fn record_pending_authoring_is_idempotent() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let marker = PendingCoherenceMarker::EdgeCurvePending {
        edge: mev.edge,
        owner_operator: "test_op",
    };

    brep.record_pending_authoring(marker);
    brep.record_pending_authoring(marker);
    brep.record_pending_authoring(marker);

    assert_eq!(brep.pending_coherence.len(), 1);
}

#[test]
fn record_pending_authoring_distinguishes_owner_operator() {
    // Same edge, different owner_operator → two separate markers
    // (the operator label is part of the marker identity).
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
        edge: mev.edge,
        owner_operator: "op_a",
    });
    brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
        edge: mev.edge,
        owner_operator: "op_b",
    });

    assert_eq!(brep.pending_coherence.len(), 2);
}

#[test]
fn clear_pending_authoring_for_edge_removes_only_matching_markers() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();

    brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
        edge: mev1.edge,
        owner_operator: "op_a",
    });
    brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
        edge: mev2.edge,
        owner_operator: "op_b",
    });
    assert_eq!(brep.pending_coherence.len(), 2);

    brep.clear_pending_authoring_for_edge(mev1.edge);

    assert_eq!(brep.pending_coherence.len(), 1);
    assert!(matches!(
        brep.pending_coherence[0],
        PendingCoherenceMarker::EdgeCurvePending { edge, .. } if edge == mev2.edge
    ));
}

#[test]
fn coherence_report_is_empty_by_default() {
    let report = CoherenceReport::default();
    assert!(report.is_empty());
    assert!(report.l1_missing_edge_curves.is_empty());
    assert!(report.l3_missing_pcurves.is_empty());
    assert!(report.l2_face_boundary_disagreements.is_empty());
}
