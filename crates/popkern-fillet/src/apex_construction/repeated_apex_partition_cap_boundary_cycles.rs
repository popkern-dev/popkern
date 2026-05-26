use super::*;

pub(super) fn collect_loop_half_edges(
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

pub(super) fn intended_cap_boundary_cycle(
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
) -> Vec<HalfEdgeHandle> {
    let mut cycle = Vec::with_capacity(2 * fillet_ends.len());
    for (fillet_end, face_edge) in fillet_ends.iter().zip(face_edges) {
        cycle.push(fillet_end.cap_half_edge);
        cycle.push(face_edge.face_edge_cap_half_edge);
    }
    cycle
}

pub(super) fn refreshed_intended_cap_boundary_cycle(
    brep: &BRep,
    cap_loop: popkern_core::topology::LoopHandle,
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
) -> Option<Vec<HalfEdgeHandle>> {
    let expected_edges = intended_cap_boundary_edges(face_edges, fillet_ends);
    let actual_cycle = collect_loop_half_edges(brep, cap_loop).ok()?;
    if actual_cycle.len() != expected_edges.len() {
        return None;
    }
    let actual_edges = actual_cycle
        .iter()
        .map(|half_edge| brep.half_edges.get(*half_edge).map(|data| data.edge))
        .collect::<Option<Vec<_>>>()?;
    for offset in 0..actual_edges.len() {
        if expected_edges
            .iter()
            .enumerate()
            .all(|(i, expected)| actual_edges[(offset + i) % actual_edges.len()] == *expected)
        {
            return Some(
                (0..actual_cycle.len())
                    .map(|i| actual_cycle[(offset + i) % actual_cycle.len()])
                    .collect(),
            );
        }
    }
    None
}

fn intended_cap_boundary_edges(
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
) -> Vec<EdgeHandle> {
    let mut edges = Vec::with_capacity(2 * fillet_ends.len());
    for (fillet_end, face_edge) in fillet_ends.iter().zip(face_edges) {
        edges.push(fillet_end.edge);
        edges.push(face_edge.edge);
    }
    edges
}

pub(super) fn reversed_twin_cycle(
    brep: &BRep,
    current_cap_cycle: &[HalfEdgeHandle],
    cap_face: FaceHandle,
) -> Result<Vec<HalfEdgeHandle>, ApexBuildError> {
    current_cap_cycle
        .iter()
        .rev()
        .map(|half_edge| {
            brep.half_edges
                .get(*half_edge)
                .map(|half_edge_data| half_edge_data.twin)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition cap-boundary rewrite references stale \
                         half-edge {half_edge:?} on cap face {cap_face:?}"
                    ),
                })
        })
        .collect()
}

pub(super) fn validate_half_edge_cycle(
    brep: &BRep,
    cycle: &[HalfEdgeHandle],
    label: &'static str,
    cap_face: FaceHandle,
) -> Result<(), ApexBuildError> {
    if cycle.is_empty() {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition {label} for cap face {cap_face:?} is empty"
            ),
        });
    }
    for i in 0..cycle.len() {
        let current = cycle[i];
        let next = cycle[(i + 1) % cycle.len()];
        let current_end = if label.contains("support")
            || label.contains("strip")
            || label.contains("target-cylinder")
        {
            support_split_half_edge_end_vertex(brep, current)
        } else {
            half_edge_end_vertex(brep, current)
        }
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition {label} for cap face {cap_face:?} has stale \
                     current half-edge {current:?}"
            ),
        })?;
        let next_start = brep
            .half_edges
            .get(next)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition {label} for cap face {cap_face:?} has stale \
                     next half-edge {next:?}"
                ),
            })?;
        let coincident_support_endpoint = (label.contains("support") || label.contains("strip"))
            && vertices_are_coincident(brep, current_end, next_start);
        if current_end != next_start && !coincident_support_endpoint {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition {label} for cap face {cap_face:?} is not \
                     contiguous: {current:?} ends at {current_end:?}, but {next:?} starts at \
                     {next_start:?}"
                ),
            });
        }
    }
    Ok(())
}

pub(super) fn assert_cap_boundary_matches_selected_cycle(
    brep: &BRep,
    cap_face: FaceHandle,
    selected_cap_cycle: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    let cap_loop = brep
        .faces
        .get(cap_face)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition selected cap boundary references stale cap face \
                 {cap_face:?}"
            ),
        })?
        .outer_loop;
    let actual = collect_loop_half_edges(brep, cap_loop)?;
    if actual == selected_cap_cycle {
        return Ok(());
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition selected cap boundary mismatch: got {:?}, expected \
             selected cycle {:?}",
            actual, selected_cap_cycle
        ),
    })
}
