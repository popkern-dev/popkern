use super::types::*;
use crate::euler::DiagnosticHandle;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, TopologyRole, VertexHandle};

#[test]
fn live_boundary_vertex_has_boundary_validity_fact() {
    let open_box = BRep::open_box(Tolerance::default()).expect("open box");
    let boundary_edge = *open_box
        .boundary_edges
        .iter()
        .next()
        .expect("boundary edge");
    let vertex = open_box.half_edges[open_box.edges[boundary_edge].half_edge].origin;

    let report = open_box.diagnose_corner_apex(request(
        vertex,
        &[boundary_edge],
        BlendKind::FilletConstantRadius,
    ));

    assert_eq!(
        report.topology.validity.boundary,
        CornerBoundaryMembership::Boundary
    );
    assert_eq!(
        report.topology.validity.radial,
        CornerRadialTopology::Manifold
    );
}

#[test]
fn live_nonmanifold_vertex_has_radial_validity_fact() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    let half_edge = cube.edges[edges[0]].half_edge;
    cube.nonmanifold_radial.insert(half_edge, half_edge);

    let report =
        cube.diagnose_corner_apex(request(vertex, &edges, BlendKind::FilletConstantRadius));

    assert_eq!(
        report.topology.validity.radial,
        CornerRadialTopology::NonManifoldRadial
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::Invalid("non-manifold-corner")
    );
}

#[test]
fn live_construction_residual_face_has_generated_context_fact() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    cube.faces[face].role = TopologyRole::ConstructionResidual;
    let edges = incident_edges(&cube, vertex);

    let report =
        cube.diagnose_corner_apex(request(vertex, &edges, BlendKind::FilletConstantRadius));

    assert_eq!(
        report.topology.validity.generated_context,
        CornerGeneratedTopologyContext::ConstructionResidual
    );
    assert_eq!(
        report.topology.surface_compatibility,
        SurfaceCompatibility::Incompatible("construction-residual-support")
    );
}

fn request(
    vertex: VertexHandle,
    edges: &[EdgeHandle],
    kind: BlendKind,
) -> CornerApexDiagnosticRequest {
    CornerApexDiagnosticRequest {
        scope: LocalCornerScope {
            vertex: DiagnosticHandle::Vertex(vertex),
        },
        selected_edges: edges
            .iter()
            .map(|edge| LocalBlendEdgeSelection {
                edge: DiagnosticHandle::Edge(*edge),
                blend_kind: kind,
                radius_nanometers: Some(10_000_000),
            })
            .collect(),
        selection_intent: "unit-test",
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

fn incident_face(brep: &BRep, vertex: VertexHandle) -> crate::topology::FaceHandle {
    brep.faces
        .iter()
        .find(|(_, face)| loop_vertices(brep, face.outer_loop).contains(&vertex))
        .map(|(handle, _)| handle)
        .expect("incident face")
}

fn loop_vertices(brep: &BRep, loop_: crate::topology::LoopHandle) -> Vec<VertexHandle> {
    let mut vertices = Vec::new();
    let start = brep.loops[loop_].half_edge;
    let mut cursor = start;
    loop {
        let half_edge = &brep.half_edges[cursor];
        vertices.push(half_edge.origin);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    vertices
}
