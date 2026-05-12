use super::*;
use crate::euler::DiagnosticHandle;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
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
