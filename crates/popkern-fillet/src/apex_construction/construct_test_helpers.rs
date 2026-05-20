use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::VertexHandle;

use super::cap_kind::CapKind;
use super::plan::ApexPlan;
use super::target_topology::{TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology};
use super::telemetry::{ConstructionRung, PlanTelemetry};

/// Build a synthetic ApexPlan whose target_topology has the supplied apex
/// vertex positions and a Sphere cap kind. Used by the DirectApex.2.b find-only tests
/// where the real sphere-fit pipeline is not relevant to the assertion.
pub(super) fn synthetic_sphere_plan(
    apex_vertex: VertexHandle,
    apex_vertices: Vec<Point3>,
) -> ApexPlan {
    let n = apex_vertices.len();
    let cap_kind = CapKind::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.1,
        tangent_arc_endpoints: apex_vertices.clone(),
    };
    let cap_face = TargetCapFace::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.1,
        outer_loop_edge_indices: (0..n).collect(),
    };
    let target_topology = TargetTopology {
        apex_vertices,
        apex_edges: Vec::new(),
        cap_face,
        face_boundary_updates: Vec::new(),
    };
    let mut telemetry = PlanTelemetry::new();
    telemetry.rung = ConstructionRung::Direct;
    ApexPlan {
        apex_vertex,
        cap_kind,
        target_topology,
        telemetry,
    }
}

/// Build a synthetic ApexPlan with explicit apex_edges connecting supplied
/// apex_vertices via degree-1 NURBS curves. Used by DirectApex.2.c find/create tests.
pub(super) fn synthetic_sphere_plan_with_edges(
    apex_vertex: VertexHandle,
    apex_vertices: Vec<Point3>,
    edge_pairs: Vec<(usize, usize)>,
) -> ApexPlan {
    let mut apex_edges: Vec<TargetEdge> = Vec::new();
    for &(origin_idx, dest_idx) in &edge_pairs {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![apex_vertices[origin_idx], apex_vertices[dest_idx]],
            weights: vec![1.0, 1.0],
        };
        apex_edges.push(TargetEdge {
            origin_vertex_idx: origin_idx,
            dest_vertex_idx: dest_idx,
            curve,
            role: TargetEdgeRole::CapVsCyl { fillet_idx: 0 },
        });
    }
    let cap_kind = CapKind::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.1,
        tangent_arc_endpoints: apex_vertices.clone(),
    };
    let cap_face = TargetCapFace::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.1,
        outer_loop_edge_indices: (0..apex_edges.len()).collect(),
    };
    let target_topology = TargetTopology {
        apex_vertices,
        apex_edges,
        cap_face,
        face_boundary_updates: Vec::new(),
    };
    let mut telemetry = PlanTelemetry::new();
    telemetry.rung = ConstructionRung::Direct;
    ApexPlan {
        apex_vertex,
        cap_kind,
        target_topology,
        telemetry,
    }
}
