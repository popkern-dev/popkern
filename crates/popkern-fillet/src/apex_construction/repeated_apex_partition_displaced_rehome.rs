//! Handover completion at end of partition: every half-edge that
//! claims a loop it is no longer reachable from gets a closed
//! `ConstructionResidual` home.
//!
//! Mid-transaction, displaced claims are deliberate hand-offs consumed
//! by later passes (measured: rehoming inside every `wire_existing_loop`
//! broke strict acceptance 8→4 cap/support seams), so this sweep runs
//! ONCE, after the partition's last loop rewrite and residual demotion.
//! What is still displaced then is genuinely orphaned; parking each
//! maximal chain in its own closed residual loop keeps the scaffold
//! loop-partitionable (the PK-092 mechanism-broken proof showed later
//! weld/retire passes cannot repair a non-partitionable scaffold).

use popkern_core::topology::TopologyRole;

use super::*;

/// Sweep the whole model for half-edges unreachable from their claimed
/// loop's anchor and re-home them into closed residual holding faces.
/// Returns how many holding faces were created.
pub(super) fn rehome_all_displaced_claimers(brep: &mut BRep) -> Result<usize, ApexBuildError> {
    let mut reachable = std::collections::BTreeSet::new();
    for (loop_handle, loop_data) in &brep.loops {
        let start = loop_data.half_edge;
        let mut cursor = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(data) = brep.half_edges.get(cursor) else {
                break;
            };
            if data.loop_ == loop_handle {
                reachable.insert(cursor);
            }
            cursor = data.next;
            if cursor == start {
                break;
            }
        }
    }
    let mut displaced: std::collections::BTreeSet<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter_map(|(heh, data)| {
            (brep.loops.contains_key(data.loop_) && !reachable.contains(&heh)).then_some(heh)
        })
        .collect();
    if displaced.is_empty() {
        return Ok(0);
    }
    let mut created = 0usize;
    while let Some(&seed) = displaced.iter().next() {
        // Walk back to the chain head (predecessor also displaced).
        let mut head = seed;
        for _ in 0..displaced.len() {
            let prev = brep.half_edges[head].prev;
            if prev != head && displaced.contains(&prev) && prev != seed {
                head = prev;
            } else {
                break;
            }
        }
        // Maximal forward chain.
        let mut chain = vec![head];
        let mut cursor = head;
        for _ in 0..displaced.len() {
            let next = brep.half_edges[cursor].next;
            if next != head && displaced.contains(&next) {
                chain.push(next);
                cursor = next;
            } else {
                break;
            }
        }
        for he in &chain {
            displaced.remove(he);
        }
        let chain_start = brep.half_edges[chain[0]].origin;
        let tail = *chain.last().ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: "displaced chain lost its tail".to_owned(),
        })?;
        let chain_end = half_edge_chain_end(brep, tail)?;
        let template_face = brep
            .loops
            .get(brep.half_edges[chain[0]].loop_)
            .map(|loop_data| loop_data.face)
            .filter(|face| brep.faces.contains_key(*face))
            .or_else(|| brep.faces.keys().next())
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: "displaced-claimer rehome found no template face".to_owned(),
            })?;
        let mut holding_cycle = chain.clone();
        if chain_end != chain_start {
            let (seam_edge, seam_side, _) = insert_unwired_edge_pair(brep, chain_end, chain_start);
            brep.edges[seam_edge].role = TopologyRole::ConstructionResidual;
            holding_cycle.push(seam_side);
        }
        let holding_face = insert_face_with_outer_loop(brep, template_face, &holding_cycle)?;
        brep.faces[holding_face].role = TopologyRole::ConstructionResidual;
        brep.faces[holding_face].surface = None;
        created += 1;
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/displaced-rehome] {} orphaned side(s) re-homed into \
                 residual {holding_face:?}",
                chain.len()
            );
        }
    }
    Ok(created)
}

fn half_edge_chain_end(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Result<VertexHandle, ApexBuildError> {
    let data = brep
        .half_edges
        .get(half_edge)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("displaced chain tail {half_edge:?} is stale"),
        })?;
    if let Some(twin) = brep.half_edges.get(data.twin) {
        return Ok(twin.origin);
    }
    brep.half_edges
        .get(data.next)
        .map(|next| next.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("displaced chain tail {half_edge:?} has no recoverable endpoint"),
        })
}
