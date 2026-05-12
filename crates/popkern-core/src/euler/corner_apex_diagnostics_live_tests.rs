use super::*;
use crate::euler::DiagnosticHandle;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, TopologyRole, VertexHandle};

#[test]
fn live_cube_corner_extracts_topology_geometry_sizing_and_strategy() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    assert_eq!(edges.len(), 3);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(report.topology.state, CornerTopologyState::Ready);
    assert_eq!(report.topology.incident_blend_edges, 3);
    assert_eq!(
        report.topology.geometry.tangent,
        CornerGeometryAvailability::Available
    );
    assert_eq!(
        report.topology.geometry.normal,
        CornerGeometryAvailability::Available
    );
    assert_eq!(report.classification.valence, BlendValence::StandardCorner);
    assert_eq!(
        report.classification.convexity,
        VertexConvexityClass::AllConvex
    );
    assert_eq!(report.sizing.state, ApexSizingState::Ready);
    assert_eq!(
        report.support_adequacy.state,
        CornerSupportAdequacyState::Fits
    );
    assert_eq!(
        report.strategy.family,
        CornerApexStrategyFamily::ThreeFaceBlend
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("three-face-tangency-solver")
    );
}

#[test]
fn live_cube_through_vertex_reports_mutual_trim_trial_readiness() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges[..2],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(report.classification.valence, BlendValence::ThroughVertex);
    assert_eq!(
        report.trial_interaction.state,
        TrialInteractionState::MissingIntersectionPrimitive
    );
    assert_eq!(report.strategy.family, CornerApexStrategyFamily::MutualTrim);
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("blend-pipe-intersection-probe")
    );
}

#[test]
fn live_pyramid_apex_reports_high_valence_missing_patch_primitive() {
    let pyramid = BRep::square_pyramid(2.0, 1.5, Tolerance::default()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, vertex)| vertex.position.z > 0.0)
        .map(|(handle, _)| handle)
        .expect("apex");
    let edges = incident_edges(&pyramid, apex);
    assert_eq!(edges.len(), 4);

    let report = pyramid.diagnose_corner_apex(request(
        apex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.classification.valence,
        BlendValence::HighValenceCorner
    );
    assert_eq!(
        report.strategy.family,
        CornerApexStrategyFamily::HighValencePatch
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("n-sided-corner-patch")
    );
}

#[test]
fn live_stale_vertex_handle_is_reported_without_panicking() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let stale = vertex_at(&cube, 0.0, 0.0, 0.0);
    cube.vertices.remove(stale);

    let report = cube.diagnose_corner_apex(request(
        stale,
        &[],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(report.topology.state, CornerTopologyState::StaleHandle);
    assert_eq!(
        report.topology.surface_compatibility,
        SurfaceCompatibility::Unknown("missing-corner-topology")
    );
    assert_eq!(
        report.lifecycle.state,
        LifecycleState::TopologyFingerprintMismatch
    );
}

#[test]
fn live_free_end_selects_cap_runout() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edge = incident_edges(&cube, vertex)[0];

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &[edge],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(report.classification.valence, BlendValence::FreeEnd);
    assert_eq!(report.strategy.family, CornerApexStrategyFamily::CapRunout);
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::Ready
    );
}

#[test]
fn live_oversized_radius_reports_setback_adjustment_dependency() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        2_000_000_000,
    ));

    assert_eq!(report.sizing.state, ApexSizingState::RadiusTooLarge);
    assert_eq!(report.strategy.family, CornerApexStrategyFamily::Setback);
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("radius-or-setback-adjustment")
    );
}

#[test]
fn live_incompatible_blend_families_block_strategy_readiness() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    let request = CornerApexDiagnosticRequest {
        scope: LocalCornerScope {
            vertex: DiagnosticHandle::Vertex(vertex),
        },
        selected_edges: vec![
            LocalBlendEdgeSelection {
                edge: DiagnosticHandle::Edge(edges[0]),
                blend_kind: BlendKind::G2Blend,
                radius_nanometers: Some(10_000_000),
            },
            LocalBlendEdgeSelection {
                edge: DiagnosticHandle::Edge(edges[1]),
                blend_kind: BlendKind::Chamfer,
                radius_nanometers: Some(10_000_000),
            },
        ],
        selection_intent: "unit-test",
        operation_context: "live-corner-apex-diagnostics",
    };

    let report = cube.diagnose_corner_apex(request);

    assert_eq!(
        report.classification.blend_types,
        BlendTypeCompatibility::Incompatible("g2-chamfer-corner")
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::Invalid("g2-chamfer-corner")
    );
}

#[test]
fn live_missing_support_surface_names_geometry_dependency() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = cube
        .faces
        .iter()
        .find(|(_, face)| {
            loop_vertices(&cube, face.outer_loop)
                .into_iter()
                .any(|loop_vertex| loop_vertex == vertex)
        })
        .map(|(handle, _)| handle)
        .expect("incident face");
    cube.faces[face].surface = None;

    let edges = incident_edges(&cube, vertex);
    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.topology.geometry.normal,
        CornerGeometryAvailability::BlockedByDependency(GeometryProbeDependency::SurfaceNormal)
    );
    assert_eq!(
        report.support_adequacy.dependency,
        Some(GeometryProbeDependency::SupportFootprint)
    );
}

#[test]
fn live_open_box_boundary_vertex_reports_cliff_termination() {
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
        10_000_000,
    ));

    assert_eq!(report.topology.state, CornerTopologyState::Boundary);
    assert_eq!(
        report.support_adequacy.state,
        CornerSupportAdequacyState::Cliff
    );
    assert_eq!(
        report.strategy.family,
        CornerApexStrategyFamily::CliffTermination
    );
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::ReadyWithManagedRisk
    );
}

#[test]
fn live_stale_selected_edge_marks_topology_stale() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let stale_edge = incident_edges(&cube, vertex)[0];
    cube.edges.remove(stale_edge);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &[stale_edge],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(report.topology.state, CornerTopologyState::StaleHandle);
    assert_eq!(report.topology.incident_blend_edges, 0);
}

#[test]
fn live_nonmanifold_radial_marker_blocks_strategy() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edges = incident_edges(&cube, vertex);
    let half_edge = cube.edges[edges[0]].half_edge;
    cube.nonmanifold_radial.insert(half_edge, half_edge);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(report.topology.state, CornerTopologyState::NonManifold);
    assert_eq!(
        report.strategy.readiness,
        CornerApexStrategyReadiness::Invalid("non-manifold-corner")
    );
}

#[test]
fn live_inner_loop_support_collision_is_reported() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    let outer_loop = cube.faces[face].outer_loop;
    cube.faces[face].inner_loops.push(outer_loop);
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.support_adequacy.state,
        CornerSupportAdequacyState::HoleOrInnerLoopCollision
    );
}

#[test]
fn live_construction_residual_support_is_incompatible() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let face = incident_face(&cube, vertex);
    cube.faces[face].role = TopologyRole::ConstructionResidual;
    let edges = incident_edges(&cube, vertex);

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &edges,
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.topology.surface_compatibility,
        SurfaceCompatibility::Incompatible("construction-residual-support")
    );
}

#[test]
fn live_near_tangent_face_normals_are_classified() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vertex = vertex_at(&cube, 0.0, 0.0, 0.0);
    let edge = incident_edges(&cube, vertex)[0];
    let faces = edge_faces_for_test(&cube, edge);
    let normal = match cube.faces[faces[0]]
        .surface
        .as_ref()
        .expect("surface")
        .geometry
    {
        SurfaceGeometry::Plane { normal, .. } => normal,
        _ => panic!("expected plane"),
    };
    if let Some(surface) = &mut cube.faces[faces[1]].surface {
        if let SurfaceGeometry::Plane { normal: n, .. } = &mut surface.geometry {
            *n = normal;
        }
    }

    let report = cube.diagnose_corner_apex(request(
        vertex,
        &[edge],
        BlendKind::FilletConstantRadius,
        10_000_000,
    ));

    assert_eq!(
        report.classification.convexity,
        VertexConvexityClass::NearTangent
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

fn incident_face(brep: &BRep, vertex: VertexHandle) -> crate::topology::FaceHandle {
    brep.faces
        .iter()
        .find(|(_, face)| loop_vertices(brep, face.outer_loop).contains(&vertex))
        .map(|(handle, _)| handle)
        .expect("incident face")
}

fn edge_faces_for_test(brep: &BRep, edge: EdgeHandle) -> Vec<crate::topology::FaceHandle> {
    let mut faces = std::collections::BTreeSet::new();
    let first = brep.edges[edge].half_edge;
    for half_edge in [first, brep.half_edges[first].twin] {
        if let Some(loop_) = brep.half_edges.get(half_edge).map(|he| he.loop_) {
            faces.insert(brep.loops[loop_].face);
        }
    }
    faces.into_iter().collect()
}
