use super::partition_spilled_corner::normalize_spilled_fillet_corners;
use super::*;

#[path = "repeated_apex_partition_support_loop_residual.rs"]
mod support_loop_residual;

pub(super) use support_loop_residual::materialize_authored_support_boundary_residual_cycles;
use support_loop_residual::materialize_support_boundary_residual_cycles;

pub(super) fn materialize_closed_support_trim_cycles(
    brep: &mut BRep,
    face_edges: &[FaceEdgeOwnership],
    trim_sides: &[SupportOwnedTrimSide],
) -> Result<(), ApexBuildError> {
    let rollback = brep.clone();
    let result = materialize_closed_support_trim_cycles_candidate(brep, face_edges, trim_sides);
    if let Err(err) = result {
        *brep = rollback;
        return Err(err);
    }
    validate_live_half_edge_cycles_after_support_rewrite(brep)?;
    Ok(())
}

fn validate_live_half_edge_cycles_after_support_rewrite(brep: &BRep) -> Result<(), ApexBuildError> {
    for (loop_handle, loop_data) in &brep.loops {
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut visited = std::collections::BTreeSet::new();
        for _ in 0..=brep.half_edges.len() {
            let half_edge =
                brep.half_edges
                    .get(cursor)
                    .ok_or_else(|| ApexBuildError::PlanMismatch {
                        reason: format!(
                            "planned repeated-apex partition support ownership rewrite produced \
                         invalid live loop {loop_handle:?}: missing half-edge {cursor:?}"
                        ),
                    })?;
            if half_edge.loop_ != loop_handle {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition support ownership rewrite produced \
                         invalid live loop {loop_handle:?}: half-edge {cursor:?} belongs to \
                         {:?}",
                        half_edge.loop_
                    ),
                });
            }
            if !visited.insert(cursor) {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition support ownership rewrite produced \
                         invalid live loop {loop_handle:?}: repeated half-edge {cursor:?}"
                    ),
                });
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
        if cursor != start {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition support ownership rewrite produced invalid \
                     live loop {loop_handle:?}: cycle did not close from {start:?}"
                ),
            });
        }
    }
    Ok(())
}

pub(super) fn repair_support_source_loops_after_cylinder_rewire(
    brep: &mut BRep,
    trim_sides: &[SupportOwnedTrimSide],
) -> Result<Vec<SupportSourceLoopTransfer>, ApexBuildError> {
    let mut transfers = Vec::new();
    for trim_side in trim_sides {
        if !trim_side.requires_support_ownership {
            continue;
        }
        let Some(source) = trim_side.support_source else {
            continue;
        };
        if !brep.loops.contains_key(source.loop_) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot repair \
                     support loop after cylinder trim ownership transfer: stale source loop \
                     {loop_handle:?}",
                    fillet_index = trim_side.fillet_index,
                    loop_handle = source.loop_
                ),
            });
        }
        let mut source = source;
        if !brep.half_edges.contains_key(source.prev)
            || !brep.half_edges.contains_key(source.next)
            || !brep.half_edges.contains_key(trim_side.half_edge)
        {
            let Some(moved) = brep.half_edges.get(trim_side.half_edge) else {
                continue;
            };
            if moved.loop_ != source.loop_ {
                continue;
            }
            source.prev = moved.prev;
            source.next = moved.next;
        }
        let prev_still_points_to_moved = brep.half_edges[source.prev].next == trim_side.half_edge;
        let next_still_points_to_moved = brep.half_edges[source.next].prev == trim_side.half_edge;
        if !prev_still_points_to_moved || !next_still_points_to_moved {
            if !support_loop_contains_half_edge(brep, source.loop_, trim_side.half_edge) {
                continue;
            }
            let moved = &brep.half_edges[trim_side.half_edge];
            if moved.loop_ != source.loop_ {
                continue;
            }
            source.prev = moved.prev;
            source.next = moved.next;
        }
        let replacement = author_support_source_loop_replacement_span(brep, trim_side, source)?;
        brep.half_edges[source.prev].next = replacement;
        brep.half_edges[source.next].prev = replacement;
        brep.half_edges[replacement].prev = source.prev;
        brep.half_edges[replacement].next = source.next;
        brep.half_edges[replacement].loop_ = source.loop_;
        if brep.loops[source.loop_].half_edge == trim_side.half_edge {
            brep.loops[source.loop_].half_edge = replacement;
        }
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/support-loop-repair] FilletEnd[{}] replaced stolen cylinder \
                 side {:?} in source support loop {:?} with boundary span {:?} between {:?}->{:?}",
                trim_side.fillet_index,
                trim_side.half_edge,
                source.loop_,
                replacement,
                source.prev,
                source.next
            );
        }
        transfers.push(SupportSourceLoopTransfer {
            loop_: source.loop_,
            replacement,
        });
    }
    Ok(transfers)
}

pub(super) struct SupportSourceLoopTransfer {
    pub(super) loop_: LoopHandle,
    pub(super) replacement: HalfEdgeHandle,
}

fn support_loop_contains_half_edge(
    brep: &BRep,
    loop_handle: LoopHandle,
    half_edge: HalfEdgeHandle,
) -> bool {
    collect_loop_half_edges(brep, loop_handle)
        .ok()
        .is_some_and(|cycle| cycle.contains(&half_edge))
}

fn author_support_source_loop_replacement_span(
    brep: &mut BRep,
    trim_side: &SupportOwnedTrimSide,
    source: SupportOwnedTrimSource,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let start = brep
        .half_edges
        .get(trim_side.half_edge)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot replace stolen support \
                 loop trim {:?}: stale moved half-edge",
                trim_side.fillet_index, trim_side.half_edge
            ),
        })?;
    let end = brep
        .half_edges
        .get(source.next)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot replace stolen support \
                 loop trim {:?}: stale next half-edge {:?}",
                trim_side.fillet_index, trim_side.half_edge, source.next
            ),
        })?;
    let curve = brep
        .half_edges
        .get(trim_side.half_edge)
        .and_then(|half_edge| brep.edges.get(half_edge.edge))
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot replace stolen support \
                 loop trim {:?}: moved edge has no curve",
                trim_side.fillet_index, trim_side.half_edge
            ),
        })?;
    let (edge, replacement) = insert_boundary_half_edge(brep, start);
    brep.edges[edge].curve = Some(curve);
    brep.edges[edge].role = TopologyRole::Production;
    if support_split_half_edge_end_vertex(brep, replacement) != Some(end) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] authored support replacement \
                 span with wrong endpoint: expected {start:?}->{end:?}, got {replacement:?}",
                trim_side.fillet_index
            ),
        });
    }
    Ok(replacement)
}

fn materialize_closed_support_trim_cycles_candidate(
    brep: &mut BRep,
    face_edges: &[FaceEdgeOwnership],
    trim_sides: &[SupportOwnedTrimSide],
) -> Result<(), ApexBuildError> {
    trace_topology_counts(brep, "support-loop-entry");
    let source_loop_transfers =
        repair_support_source_loops_after_cylinder_rewire(brep, trim_sides)?;
    let mut authored_boundary_spans = source_loop_transfers
        .iter()
        .map(|transfer| transfer.replacement)
        .collect::<Vec<_>>();
    let mut by_support_face: std::collections::BTreeMap<FaceHandle, Vec<HalfEdgeHandle>> =
        std::collections::BTreeMap::new();
    let mut generated_trim_count_by_support_face =
        std::collections::BTreeMap::<FaceHandle, usize>::new();
    let mut support_faces_with_existing_boundary = std::collections::BTreeSet::new();
    for trim_side in trim_sides {
        if !trim_side.requires_support_ownership {
            continue;
        }
        if brep
            .half_edges
            .get(trim_side.half_edge)
            .and_then(|half_edge| brep.loops.get(half_edge.loop_))
            .is_none_or(|loop_data| loop_data.face != trim_side.cylinder_face)
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot plan \
                     support ownership for {cylinder_half_edge:?}: it is not owned by target \
                     cylinder face {cylinder_face:?}",
                    fillet_index = trim_side.fillet_index,
                    cylinder_half_edge = trim_side.half_edge,
                    cylinder_face = trim_side.cylinder_face
                ),
            });
        }
        let twin = brep
            .half_edges
            .get(trim_side.half_edge)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot plan \
                     support ownership for stale half-edge {cylinder_half_edge:?}",
                    fillet_index = trim_side.fillet_index,
                    cylinder_half_edge = trim_side.half_edge
                ),
            })?;
        if trim_side.support_side_owned {
            continue;
        }
        by_support_face
            .entry(trim_side.support_face)
            .or_default()
            .push(twin);
        *generated_trim_count_by_support_face
            .entry(trim_side.support_face)
            .or_default() += 1;
    }
    for face_edge in face_edges {
        if let Some(half_edges) = by_support_face.get_mut(&face_edge.support_face)
            && !half_edge_is_degenerate(brep, face_edge.support_half_edge)
        {
            half_edges.push(face_edge.support_half_edge);
            support_faces_with_existing_boundary.insert(face_edge.support_face);
        }
    }

    let mut support_groups = by_support_face.into_iter().collect::<Vec<_>>();
    support_groups.sort_by_key(|(support_face, _)| {
        let generated_count = generated_trim_count_by_support_face
            .get(support_face)
            .copied()
            .unwrap_or_default();
        (
            std::cmp::Reverse(generated_count),
            support_faces_with_existing_boundary.contains(support_face),
            *support_face,
        )
    });

    for (support_face, mut half_edges) in support_groups {
        trace_topology_counts(brep, "support-face-before");
        half_edges.retain(|half_edge| !half_edge_is_degenerate(brep, *half_edge));
        let support_split_plan = materialize_support_split_cycles(brep, support_face, &half_edges)?;
        authored_boundary_spans.extend(support_split_plan.authored_boundary_spans.iter().copied());
        if !support_split_plan.is_complete() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot materialize support ownership on \
                     {support_face:?}: generated trim sides {:?} do not form closed support \
                     cycle(s) yet. Interior support subdivision must add the missing support \
                     boundary spans before these cylinder trims can be exposed. {}; {}",
                    half_edges,
                    support_cycle_open_chain_diagnostic(brep, &half_edges),
                    support_split_plan.diagnostic()
                ),
            });
        }
        if support_split_plan.cycles.len() != 1 {
            let counts = topology_count_summary(brep);
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot materialize support ownership on \
                     {support_face:?}: got {} closed support cycles. Multi-cycle support \
                     subdivision requires an Euler-balanced residual face/ring operator; direct \
                     inner-loop wiring or duplicate-face insertion would violate shell topology. \
                     Current topology deficit before this support face is {}; one face-loop per \
                     closed support cycle would add only {}.",
                    support_split_plan.cycles.len(),
                    counts.deficit_to_genus_zero(),
                    support_split_plan.cycles.len()
                ),
            });
        }
        let cycle = &support_split_plan.cycles[0];
        validate_half_edge_cycle(brep, cycle, "materialized support trim cycle", support_face)?;
        let loop_handle = brep
            .faces
            .get(support_face)
            .map(|face| face.outer_loop)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot replace stale support face \
                     {support_face:?}"
                ),
            })?;
        wire_existing_loop(brep, loop_handle, support_face, cycle)?;
        trace_topology_counts(brep, "support-face-after");
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/support-loop-replace] support {support_face:?} -> \
                 rewired outer loop {:?} using {:?}",
                loop_handle, cycle
            );
        }
    }
    retire_consumed_support_source_loops(brep, &source_loop_transfers)?;
    normalize_spilled_fillet_corners(brep, trim_sides, &source_loop_transfers)?;
    materialize_support_boundary_residual_cycles(brep, &authored_boundary_spans)?;
    Ok(())
}

fn retire_consumed_support_source_loops(
    brep: &mut BRep,
    transfers: &[SupportSourceLoopTransfer],
) -> Result<(), ApexBuildError> {
    let loops = transfers
        .iter()
        .map(|transfer| transfer.loop_)
        .collect::<std::collections::BTreeSet<_>>();
    for loop_handle in loops {
        let Some(loop_data) = brep.loops.get(loop_handle).cloned() else {
            continue;
        };
        if brep
            .half_edges
            .iter()
            .any(|(_, half_edge)| half_edge.loop_ == loop_handle)
        {
            continue;
        }
        let Some(face_data) = brep.faces.get(loop_data.face) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot retire consumed support source loop \
                     {loop_handle:?}: owner face {:?} is stale",
                    loop_data.face
                ),
            });
        };
        if face_data.outer_loop != loop_handle || !face_data.inner_loops.is_empty() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot retire consumed support source loop \
                     {loop_handle:?}: owner face {:?} has non-simple loop ownership",
                    loop_data.face
                ),
            });
        }
        let shell = face_data.shell;
        brep.loops.remove(loop_handle);
        brep.faces.remove(loop_data.face);
        if brep
            .shells
            .get(shell)
            .is_some_and(|shell_data| shell_data.face == loop_data.face)
            && let Some(replacement_face) = brep
                .faces
                .iter()
                .find_map(|(face, data)| (data.shell == shell).then_some(face))
        {
            brep.shells[shell].face = replacement_face;
        }
    }
    Ok(())
}

fn trace_topology_counts(brep: &BRep, label: &'static str) {
    if !crate::direct_apex_trace::enabled() {
        return;
    }
    eprintln!("{}", topology_count_summary(brep).trace_label(label));
}
