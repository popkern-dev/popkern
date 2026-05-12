use super::*;
use crate::euler::DiagnosticHandle;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

#[test]
fn live_three_face_tangency_rejects_variable_radius_request() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(variable_radius_request(vertex, &edges));

    assert_not_ready(
        report.strategy.readiness,
        "variable-radius-three-face-tangency-solver",
    );
}

#[test]
fn live_three_face_tangency_rejects_singular_analytic_normals() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let faces = incident_faces(&cube, vertex);
    cube.faces[faces[0]]
        .surface
        .as_mut()
        .expect("surface")
        .geometry = SurfaceGeometry::Sphere {
        center: cube.vertices[vertex].position,
        radius: 1.0,
    };

    let report = cube.diagnose_corner_apex(request(vertex, &incident_edges(&cube, vertex)));

    assert_not_ready(
        report.strategy.readiness,
        "three-face-tangency-singular-normal",
    );
}

#[test]
fn live_three_face_tangency_rejects_impossible_offsets() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let faces = incident_faces(&cube, vertex);
    cube.faces[faces[0]]
        .surface
        .as_mut()
        .expect("surface")
        .geometry = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 0.0, 0.0),
        radius: 1.0,
    };
    cube.faces[faces[1]]
        .surface
        .as_mut()
        .expect("surface")
        .geometry = SurfaceGeometry::Sphere {
        center: Point3::new(3.0, 0.0, 0.0),
        radius: 1.0,
    };
    cube.faces[faces[2]]
        .surface
        .as_mut()
        .expect("surface")
        .geometry = z_plane();

    let report = cube.diagnose_corner_apex(request(vertex, &incident_edges(&cube, vertex)));

    assert_not_ready(report.strategy.readiness, "three-face-tangency-no-solution");
}

fn assert_not_ready(readiness: CornerApexStrategyReadiness, expected_reason: &'static str) {
    assert_eq!(
        readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive(expected_reason)
    );
    assert_ne!(readiness, CornerApexStrategyReadiness::Ready);
}

fn request(vertex: VertexHandle, edges: &[EdgeHandle]) -> CornerApexDiagnosticRequest {
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
        selection_intent: "unit-test",
        operation_context: "live-corner-apex-diagnostics",
    }
}

fn variable_radius_request(
    vertex: VertexHandle,
    edges: &[EdgeHandle],
) -> CornerApexDiagnosticRequest {
    let mut request = request(vertex, edges);
    if let Some(edge) = request.selected_edges.get_mut(0) {
        edge.radius_nanometers = Some(20_000_000);
    }
    request
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

fn incident_faces(brep: &BRep, vertex: VertexHandle) -> Vec<crate::topology::FaceHandle> {
    brep.faces
        .iter()
        .filter(|(_, face)| loop_vertices(brep, face.outer_loop).contains(&vertex))
        .map(|(handle, _)| handle)
        .collect()
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

fn z_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}
