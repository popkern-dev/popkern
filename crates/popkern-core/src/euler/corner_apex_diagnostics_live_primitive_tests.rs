use super::*;
use crate::euler::DiagnosticHandle;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

#[test]
fn live_mutual_trim_names_nonlinear_spine_gap() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    cube.edges[edges[0]].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: std::f64::consts::FRAC_PI_2,
    });

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges[..2],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(report.strategy.family, CornerApexStrategyFamily::MutualTrim);
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("nonlinear-spine-mutual-trim-probe")
    );
}

#[test]
fn live_mutual_trim_names_variable_radius_gap() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    let mut request = request(
        vertex,
        &edges[..2],
        BlendKind::FilletConstantRadius,
        10_000_000,
    );
    request.selected_edges[1].radius_nanometers = Some(20_000_000);

    let report = cube.diagnose_corner_apex(request);

    assert_eq!(report.strategy.family, CornerApexStrategyFamily::MutualTrim);
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("variable-radius-mutual-trim-probe")
    );
}

#[test]
fn live_mutual_trim_keeps_line_spine_planar_support_ready() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges[..2],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(report.strategy.family, CornerApexStrategyFamily::MutualTrim);
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::Ready
    );
}

#[test]
fn live_mutual_trim_names_nonplanar_support_gap() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    cube.faces[face].surface.as_mut().expect("surface").geometry = SurfaceGeometry::Cylinder {
        origin: Point3::new(1.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges[..2],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive(
            "nonplanar-support-mutual-trim-probe"
        )
    );
}

#[test]
fn live_mutual_trim_names_near_tangent_branch() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    let direction = edge_line_direction(&cube, edges[0]);
    cube.edges[edges[1]].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: cube.vertices[vertex].position,
            direction,
        },
        t_start: 0.0,
        t_end: 1.0,
    });

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges[..2],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("near-tangent-mutual-trim")
    );
}

#[test]
fn live_mutual_trim_names_missing_intersection_refinement() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    cube.edges[edges[0]].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: cube.vertices[vertex].position,
            direction: Vec3::zero(),
        },
        t_start: 0.0,
        t_end: 1.0,
    });

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges[..2],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive(
            "missing-mutual-trim-intersection-refinement"
        )
    );
}

fn request(
    vertex: VertexHandle,
    edges: &[EdgeHandle],
    kind: BlendKind,
    radius_nanometers: i64,
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
                radius_nanometers: Some(radius_nanometers),
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

fn edge_line_direction(brep: &BRep, edge: EdgeHandle) -> Vec3 {
    match brep.edges[edge].curve.as_ref().expect("curve").geometry {
        CurveGeometry::Line { direction, .. } => direction,
        _ => panic!("fixture edge is linear"),
    }
}
