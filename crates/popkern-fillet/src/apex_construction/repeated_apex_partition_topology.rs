//! Low-level topology helpers for planned repeated-apex partitioning.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, TopologyRole,
    VertexHandle,
};

use super::error::ApexBuildError;

pub(super) fn oriented_edge_half_edges_on_face(
    brep: &BRep,
    face: FaceHandle,
    edge: EdgeHandle,
    expected_start: VertexHandle,
    expected_end: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    let edge_data = brep.edges.get(edge)?;
    for he in [
        edge_data.half_edge,
        brep.half_edges.get(edge_data.half_edge)?.twin,
    ] {
        let half_edge = brep.half_edges.get(he)?;
        let twin = half_edge.twin;
        let half_edge_face = brep.loops.get(half_edge.loop_)?.face;
        if half_edge_face != face {
            continue;
        }
        let end = brep.half_edges.get(twin)?.origin;
        if half_edge.origin == expected_start && end == expected_end {
            return Some((he, twin));
        }
    }
    None
}

pub(super) fn insert_unwired_edge_pair(
    brep: &mut BRep,
    start: VertexHandle,
    end: VertexHandle,
) -> (EdgeHandle, HalfEdgeHandle, HalfEdgeHandle) {
    let edge = brep.edges.insert(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: TopologyRole::ConstructionResidual,
    });
    let half_edge = brep.half_edges.insert(HalfEdge {
        origin: start,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge,
        loop_: Default::default(),
        pcurve: None,
    });
    let twin = brep.half_edges.insert(HalfEdge {
        origin: end,
        twin: half_edge,
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge,
        loop_: Default::default(),
        pcurve: None,
    });
    brep.half_edges[half_edge].twin = twin;
    brep.edges[edge].half_edge = half_edge;
    if brep.vertices[start].outgoing == HalfEdgeHandle::default() {
        brep.vertices[start].outgoing = half_edge;
    }
    if brep.vertices[end].outgoing == HalfEdgeHandle::default() {
        brep.vertices[end].outgoing = twin;
    }
    (edge, half_edge, twin)
}

pub(super) fn insert_boundary_half_edge(
    brep: &mut BRep,
    start: VertexHandle,
) -> (EdgeHandle, HalfEdgeHandle) {
    let edge = brep.edges.insert(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: TopologyRole::ConstructionResidual,
    });
    let half_edge = brep.half_edges.insert(HalfEdge {
        origin: start,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge,
        loop_: Default::default(),
        pcurve: None,
    });
    brep.edges[edge].half_edge = half_edge;
    brep.boundary_edges.insert(edge);
    if brep.vertices[start].outgoing == HalfEdgeHandle::default() {
        brep.vertices[start].outgoing = half_edge;
    }
    (edge, half_edge)
}

fn half_edge_destination_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let data = brep.half_edges.get(half_edge)?;
    if let Some(twin) = brep.half_edges.get(data.twin) {
        return Some(twin.origin);
    }
    if !brep.boundary_edges.contains(&data.edge) {
        return brep
            .half_edges
            .get(data.next)
            .map(|next_half_edge| next_half_edge.origin);
    }
    let edge = brep.edges.get(data.edge)?;
    let Some(curve) = edge.curve.as_ref() else {
        return brep
            .half_edges
            .get(data.next)
            .map(|next_half_edge| next_half_edge.origin);
    };
    let start_position = brep.vertices.get(data.origin)?.position;
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    let destination = if start_position.distance_to(&p0) <= start_position.distance_to(&p1) {
        p1
    } else {
        p0
    };
    let tolerance = brep.tolerance.linear().max(1.0e-9);
    brep.vertices
        .iter()
        .filter_map(|(vertex, data)| {
            let distance = data.position.distance_to(&destination);
            (distance <= tolerance).then_some((vertex, distance))
        })
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(vertex, _)| vertex)
        .or_else(|| {
            brep.half_edges
                .get(data.next)
                .map(|next_half_edge| next_half_edge.origin)
        })
}

pub(super) fn insert_face_with_outer_loop(
    brep: &mut BRep,
    template_face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
) -> Result<FaceHandle, ApexBuildError> {
    if half_edges.is_empty() {
        return Err(ApexBuildError::PlanMismatch {
            reason: "planned repeated-apex partition cannot create an empty loop".to_string(),
        });
    }
    let template =
        brep.faces
            .get(template_face)
            .cloned()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition template face {template_face:?} is dead"
                ),
            })?;
    let face = brep.faces.insert(Face {
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        shell: template.shell,
        surface: template.surface,
        name: None,
        tolerance: template.tolerance,
        role: Default::default(),
    });
    if brep
        .shells
        .get(template.shell)
        .is_some_and(|shell| !brep.faces.contains_key(shell.face))
    {
        brep.shells[template.shell].face = face;
    }
    let loop_handle = brep.loops.insert(Loop {
        half_edge: half_edges[0],
        face,
    });
    brep.faces[face].outer_loop = loop_handle;
    wire_existing_loop(brep, loop_handle, face, half_edges)?;
    Ok(face)
}

pub(crate) fn repair_stale_shell_anchors(brep: &mut BRep) -> Result<(), ApexBuildError> {
    let shells = brep.shells.keys().collect::<Vec<_>>();
    for shell in shells {
        let Some(anchor) = brep.shells.get(shell).map(|shell_data| shell_data.face) else {
            continue;
        };
        if brep.faces.contains_key(anchor) {
            continue;
        }
        let replacement = brep
            .faces
            .iter()
            .find_map(|(face, data)| (data.shell == shell).then_some(face))
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot repair stale shell anchor for \
                     {shell:?}: no live replacement face remains"
                ),
            })?;
        brep.shells[shell].face = replacement;
    }
    Ok(())
}

pub(super) fn wire_existing_loop(
    brep: &mut BRep,
    loop_handle: popkern_core::topology::LoopHandle,
    face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    if half_edges.is_empty() {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot wire empty loop {loop_handle:?}"
            ),
        });
    }
    let n = half_edges.len();
    for i in 0..n {
        let current = half_edges[i];
        let next = half_edges[(i + 1) % n];
        let prev = half_edges[(i + n - 1) % n];
        let next_start = brep
            .half_edges
            .get(next)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition loop {loop_handle:?} references stale \
                     next half-edge {next:?}"
                ),
            })?;
        let current_end = brep
            .half_edges
            .get(current)
            .and_then(|_| half_edge_destination_vertex(brep, current))
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition loop {loop_handle:?} references stale \
                     half-edge {current:?}"
                ),
            })?;
        let coincident_endpoint = brep
            .vertices
            .get(current_end)
            .zip(brep.vertices.get(next_start))
            .is_some_and(|(current_vertex, next_vertex)| {
                current_vertex.position.distance_to(&next_vertex.position)
                    <= brep.tolerance.linear().max(1.0e-9)
            });
        if current_end != next_start && !coincident_endpoint {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition loop {loop_handle:?} is not contiguous at \
                     {current:?}: edge ends at {current_end:?}, next starts at {next_start:?}"
                ),
            });
        }
        let half_edge =
            brep.half_edges
                .get_mut(current)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition loop {loop_handle:?} references stale \
                     half-edge {current:?}"
                    ),
                })?;
        half_edge.next = next;
        half_edge.prev = prev;
        half_edge.loop_ = loop_handle;
        let edge = half_edge.edge;
        if let Some(edge_data) = brep.edges.get_mut(edge) {
            edge_data.role = if brep
                .faces
                .get(face)
                .is_some_and(|face_data| face_data.role == TopologyRole::ConstructionResidual)
            {
                TopologyRole::ConstructionResidual
            } else {
                TopologyRole::Production
            };
        }
    }
    brep.loops[loop_handle].half_edge = half_edges[0];
    brep.loops[loop_handle].face = face;
    Ok(())
}

pub(super) fn validate_partition_cap_to_cylinder_seams(
    brep: &BRep,
    cap_face: FaceHandle,
    n_fillets: usize,
) -> Result<(), ApexBuildError> {
    let mut cylinder_faces = std::collections::BTreeSet::new();
    let mut cylinder_seams = 0usize;
    for (_, edge) in &brep.edges {
        let Some((face_a, face_b)) = edge_incident_faces(brep, edge.half_edge) else {
            continue;
        };
        if face_a == face_b {
            continue;
        }
        let neighbor_face = if face_a == cap_face {
            face_b
        } else if face_b == cap_face {
            face_a
        } else {
            continue;
        };
        let Some(surface) = brep
            .faces
            .get(neighbor_face)
            .and_then(|face| face.surface.as_ref())
        else {
            continue;
        };
        if matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }) {
            cylinder_seams += 1;
            cylinder_faces.insert(neighbor_face);
        }
    }

    if cylinder_faces.len() < n_fillets {
        count_materialized_exact_boundary_cap_cylinder_seams(
            brep,
            cap_face,
            &mut cylinder_faces,
            &mut cylinder_seams,
        );
    }

    if cylinder_faces.len() < n_fillets {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition exposed {cylinder_seams} cap-to-cylinder seam(s) \
                 across {} distinct cylinder face(s); expected at least {n_fillets}",
                cylinder_faces.len()
            ),
        });
    }
    Ok(())
}

fn count_materialized_exact_boundary_cap_cylinder_seams(
    brep: &BRep,
    cap_face: FaceHandle,
    cylinder_faces: &mut std::collections::BTreeSet<FaceHandle>,
    cylinder_seams: &mut usize,
) {
    let Some(cap_loop) = brep.faces.get(cap_face).map(|face| face.outer_loop) else {
        return;
    };
    let Ok(cap_sides) = collect_loop_half_edges(brep, cap_loop) else {
        return;
    };
    for (face, face_data) in &brep.faces {
        if face == cap_face || cylinder_faces.contains(&face) {
            continue;
        }
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if !matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        let Ok(cylinder_sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        if cap_sides.iter().any(|cap_side| {
            cylinder_sides.iter().any(|cylinder_side| {
                same_oriented_exact_boundary_side(brep, *cap_side, *cylinder_side)
            })
        }) {
            *cylinder_seams += 1;
            cylinder_faces.insert(face);
        }
    }
}

fn same_oriented_exact_boundary_side(
    brep: &BRep,
    cap_side: HalfEdgeHandle,
    cylinder_side: HalfEdgeHandle,
) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(cylinder_start) = brep.half_edges.get(cylinder_side).map(|data| data.origin) else {
        return false;
    };
    if cap_start != cylinder_start {
        return false;
    }
    let Some(cap_end) = half_edge_destination_vertex(brep, cap_side) else {
        return false;
    };
    let Some(cylinder_end) = half_edge_destination_vertex(brep, cylinder_side) else {
        return false;
    };
    cap_end == cylinder_end
}

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Result<Vec<HalfEdgeHandle>, ApexBuildError> {
    let loop_data = brep
        .loops
        .get(loop_handle)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot collect stale loop {loop_handle:?}"
            ),
        })?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition loop {loop_handle:?} starts at stale half-edge \
                 {start:?}"
            ),
        });
    }
    let mut cycle = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        cycle.push(cursor);
        let next = brep
            .half_edges
            .get(cursor)
            .map(|half_edge| half_edge.next)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition loop {loop_handle:?} references stale \
                     half-edge {cursor:?}"
                ),
            })?;
        cursor = next;
        if cursor == start {
            return Ok(cycle);
        }
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition loop {loop_handle:?} did not close within {} \
             half-edge steps",
            brep.half_edges.len()
        ),
    })
}

fn edge_incident_faces(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<(FaceHandle, FaceHandle)> {
    let side_a = brep.half_edges.get(half_edge)?;
    let side_b = brep.half_edges.get(side_a.twin)?;
    let face_a = brep
        .loops
        .get(side_a.loop_)
        .map(|loop_data| loop_data.face)?;
    let face_b = brep
        .loops
        .get(side_b.loop_)
        .map(|loop_data| loop_data.face)?;
    Some((face_a, face_b))
}
