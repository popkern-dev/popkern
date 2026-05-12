use super::*;
use crate::euler::DiagnosticHandle;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

#[test]
fn live_selection_intent_is_separate_from_topology() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);

    let explicit = cube.diagnose_corner_apex(request(vertex, &edges, "explicit-edge-set"));
    let propagated = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        "inferred-tangent-chain-propagation",
    ));

    assert_eq!(
        explicit.topology.neighborhood.incident_edges,
        propagated.topology.neighborhood.incident_edges
    );
    assert_eq!(
        explicit.selection.source,
        CornerSelectionSource::ExplicitEdgeSet
    );
    assert_eq!(
        explicit.selection.propagation,
        CornerSelectionPropagation::NotRequested
    );
    assert_eq!(
        propagated.selection.source,
        CornerSelectionSource::InferredPropagation
    );
    assert_eq!(
        propagated.selection.propagation,
        CornerSelectionPropagation::TangentChain
    );
}

#[test]
fn live_selection_intent_records_requested_options_and_law_count() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    let mut request = request(vertex, &edges[..2], "explicit-setback-cap");
    request.selected_edges[0].blend_kind = BlendKind::FilletVariableRadius;

    let report = cube.diagnose_corner_apex(request);

    assert_eq!(
        report.selection.requested_setback,
        CornerSelectionOptionState::Requested
    );
    assert_eq!(
        report.selection.requested_cap,
        CornerSelectionOptionState::Requested
    );
    assert_eq!(report.selection.per_edge_law_assignments, 1);
}

fn request(
    vertex: VertexHandle,
    edges: &[EdgeHandle],
    selection_intent: &'static str,
) -> CornerApexDiagnosticRequest {
    CornerApexDiagnosticRequest {
        scope: LocalCornerScope {
            vertex: DiagnosticHandle::Vertex(vertex),
        },
        selected_edges: edges
            .iter()
            .map(|edge| LocalBlendEdgeSelection {
                edge: DiagnosticHandle::Edge(*edge),
                blend_kind: BlendKind::FilletConstantRadius,
                radius_nanometers: Some(10_000_000),
            })
            .collect(),
        selection_intent,
        operation_context: "live-corner-apex-diagnostics",
    }
}

fn vertex_at(brep: &BRep, x: f64, y: f64, z: f64) -> VertexHandle {
    brep.vertices
        .iter()
        .find(|(_, vertex)| {
            vertex.position.x == x && vertex.position.y == y && vertex.position.z == z
        })
        .map(|(handle, _)| handle)
        .expect("vertex")
}

fn incident_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    let mut edges = std::collections::BTreeSet::new();
    for (_, half_edge) in &brep.half_edges {
        if half_edge.origin == vertex
            || brep
                .half_edges
                .get(half_edge.twin)
                .is_some_and(|twin| twin.origin == vertex)
        {
            edges.insert(half_edge.edge);
        }
    }
    edges.into_iter().collect()
}
