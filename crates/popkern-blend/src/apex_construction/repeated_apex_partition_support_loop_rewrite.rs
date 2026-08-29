use super::partition_spilled_corner::normalize_spilled_fillet_corners;
use super::*;

#[path = "repeated_apex_partition_support_loop_assembly.rs"]
mod support_loop_assembly;
#[path = "repeated_apex_partition_support_loop_residual.rs"]
mod support_loop_residual;
#[path = "repeated_apex_partition_support_loop_routing.rs"]
mod support_loop_routing;

use support_loop_assembly::assemble_extended_facet_cycle;
pub(super) use support_loop_routing::{
    SupportSourceLoopTransfer, repair_support_source_loops_after_cylinder_rewire,
};
use support_loop_routing::{
    exchange_bridge_side_with_cap_chain, face_edge_endpoints_fit_face_plane,
    splice_side_out_of_donor_loop, winding_honest_bridge_side,
};

pub(super) use support_loop_residual::materialize_authored_support_boundary_residual_cycles;
use support_loop_residual::materialize_support_boundary_residual_cycles;

pub(super) fn materialize_closed_support_trim_cycles(
    brep: &mut BRep,
    face_edges: &[FaceEdgeOwnership],
    trim_sides: &[SupportOwnedTrimSide],
    route_face_edges_to_coplanar_facets: bool,
) -> Result<(), ApexBuildError> {
    let rollback = brep.clone();
    let result = materialize_closed_support_trim_cycles_candidate(
        brep,
        face_edges,
        trim_sides,
        route_face_edges_to_coplanar_facets,
    );
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
    route_face_edges_to_coplanar_facets: bool,
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
        // Strict V-free routing: the FaceEdge's neighbor side belongs to
        // the support group whose facet plane carries both endpoints —
        // the retiring K-P fragment it currently sits in is merely
        // coplanar scaffolding. Index arithmetic cannot express this
        // pairing (K-P corners and facets adjoin fillets differently).
        let group_key = if route_face_edges_to_coplanar_facets {
            by_support_face
                .keys()
                .copied()
                .find(|support_face| {
                    face_edge_endpoints_fit_face_plane(
                        brep,
                        *support_face,
                        face_edge.support_half_edge,
                    )
                })
                .unwrap_or(face_edge.support_face)
        } else {
            face_edge.support_face
        };
        if let Some(half_edges) = by_support_face.get_mut(&group_key)
            && !half_edge_is_degenerate(brep, face_edge.support_half_edge)
        {
            let selected =
                winding_honest_bridge_side(brep, face_edge.support_half_edge, half_edges);
            if selected == face_edge.support_half_edge {
                if group_key != face_edge.support_face {
                    splice_side_out_of_donor_loop(brep, face_edge.support_half_edge);
                }
            } else {
                exchange_bridge_side_with_cap_chain(brep, selected, face_edge.support_half_edge);
            }
            if group_key != face_edge.support_face
                && let Some(face) = brep.faces.get_mut(face_edge.support_face)
            {
                face.role = TopologyRole::ConstructionResidual;
                face.surface = None;
            }
            half_edges.push(selected);
            support_faces_with_existing_boundary.insert(group_key);
        }
    }

    if crate::direct_apex_trace::enabled() {
        let transfer_loops: std::collections::BTreeSet<LoopHandle> = source_loop_transfers
            .iter()
            .map(|transfer| transfer.loop_)
            .collect();
        for (support_face, sides) in &by_support_face {
            for side in sides {
                let Some(data) = brep.half_edges.get(*side) else {
                    continue;
                };
                let end = brep.half_edges.get(data.twin).map(|twin| twin.origin);
                let twin_owner = brep
                    .half_edges
                    .get(data.twin)
                    .map(|twin| twin.loop_)
                    .and_then(|loop_handle| {
                        brep.loops.get(loop_handle).map(|loop_data| {
                            (
                                loop_handle,
                                loop_data.face,
                                brep.faces.get(loop_data.face).map(|face| face.role),
                                transfer_loops.contains(&loop_handle),
                            )
                        })
                    });
                eprintln!(
                    "[direct-apex-trace §14.9/group-ownership] support {support_face:?} side \
                     {side:?} {:?}->{end:?} twin {:?} owned by {twin_owner:?} (loop, face, role, \
                     consumed-source-transfer)",
                    data.origin, data.twin
                );
            }
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
        if let Some(spliced) = assemble_extended_facet_cycle(brep, support_face, &half_edges)? {
            let loop_handle = brep
                .faces
                .get(support_face)
                .map(|face| face.outer_loop)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition cannot replace stale support face                          {support_face:?}"
                    ),
                })?;
            let displaced_old: Vec<HalfEdgeHandle> = collect_loop_half_edges(brep, loop_handle)
                .unwrap_or_default()
                .into_iter()
                .filter(|side| !spliced.contains(side))
                .collect();
            wire_existing_loop(brep, loop_handle, support_face, &spliced)?;
            for side in &displaced_old {
                let Some(data) = brep.half_edges.get(*side) else {
                    continue;
                };
                let side_start = data.origin;
                let Some(side_end) = brep.half_edges.get(data.twin).map(|twin| twin.origin) else {
                    continue;
                };
                let mut holding = vec![*side];
                if side_end != side_start {
                    let (seam_edge, seam_side, _) =
                        insert_unwired_edge_pair(brep, side_end, side_start);
                    brep.edges[seam_edge].role = TopologyRole::ConstructionResidual;
                    holding.push(seam_side);
                }
                let remainder_face = insert_face_with_outer_loop(brep, support_face, &holding)?;
                brep.faces[remainder_face].role = TopologyRole::ConstructionResidual;
                brep.faces[remainder_face].surface = None;
            }
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/support-loop-replace] support {support_face:?}                      explicit-assembled ({} sides; {} displaced retired)",
                    spliced.len(),
                    displaced_old.len()
                );
            }
            trace_topology_counts(brep, "support-face-after");
            continue;
        }
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
