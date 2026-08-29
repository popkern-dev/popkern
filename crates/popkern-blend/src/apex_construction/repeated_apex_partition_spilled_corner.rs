//! Spilled-fillet base-corner normalization.
//!
//! When a fillet strip ends on the body's base, the far-side trim mef
//! splits the support face along the exact cylinder∩support arc and the
//! cylinder loop replacement steals the arc's cylinder side, leaving the
//! main support fragment rimmed by a twinless curve-clone while the
//! corner fragment (the "pocket") keeps the arc's real support side.
//! The OCC reference for this configuration (BRepFilletAPI on a diamond
//! pyramid slant edge) shows the corner fragment is removed material:
//! the support face's boundary runs polygon-edge → anchor → arc →
//! anchor → polygon-edge and shares the arc edge directly with the
//! strip; no pocket face exists.
//!
//! This pass performs the Euler-neutral half of that contract: it swaps
//! the arc's real support side into the main support loop (replacing
//! the clone) and hands the clone to the corner fragment, then retires
//! the corner fragment as `ConstructionResidual`. Production scans and
//! the defect census exclude residual records, so the base↔strip seam
//! becomes a true twin pair and the pocket stops shadowing the support
//! plane.

use popkern_core::topology::TopologyRole;

use super::partition_support_loop_rewrite::SupportSourceLoopTransfer;
use super::*;

/// Swap real arc support sides into the support loops rimmed by clones
/// and retire the corner fragments; returns how many corners were
/// normalized. Every step is structurally guarded; corners that don't
/// match the spilled signature are left untouched.
pub(super) fn normalize_spilled_fillet_corners(
    brep: &mut BRep,
    trim_sides: &[SupportOwnedTrimSide],
    transfers: &[SupportSourceLoopTransfer],
) -> Result<usize, ApexBuildError> {
    let mut normalized = 0usize;
    // Arrangement A — the mef chain kept the real arc support side in
    // the main support face and parked the twinless clone in the corner
    // fragment: the fragment only needs retirement.
    for transfer in transfers {
        let clone = transfer.replacement;
        let Some(clone_side) = brep.half_edges.get(clone) else {
            trace_decline(clone, "clone side is stale");
            continue;
        };
        if !half_edge_is_unpaired_boundary(brep, clone_side.twin) {
            trace_decline(clone, "clone twin is paired");
            continue;
        }
        let pocket_loop = clone_side.loop_;
        let Some(pocket_face) = brep.loops.get(pocket_loop).map(|l| l.face) else {
            trace_decline(clone, "clone loop is stale");
            continue;
        };
        if brep
            .faces
            .get(pocket_face)
            .is_none_or(|face| face.role == TopologyRole::ConstructionResidual)
        {
            trace_decline(clone, "owner face already residual or stale");
            continue;
        }
        if !pocket_signature_matches(brep, pocket_loop, clone) {
            let cycle_len = collect_loop_half_edges(brep, pocket_loop)
                .map(|cycle| cycle.len())
                .unwrap_or(0);
            trace_decline(
                clone,
                &format!(
                    "owner loop lacks the corner signature (face {pocket_face:?}, {cycle_len} \
                     sides)"
                ),
            );
            continue;
        }
        if let Some(face) = brep.faces.get_mut(pocket_face) {
            face.role = TopologyRole::ConstructionResidual;
            face.surface = None;
        }
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/spilled-corner] corner fragment {pocket_face:?} \
                 (rimmed by clone {clone:?}) retired; support face keeps the true arc seam"
            );
        }
        normalized += 1;
    }
    // Arrangement B — the clone landed in the main support loop and the
    // corner fragment holds the real arc side: swap them, then retire.
    for trim_side in trim_sides {
        if !trim_side.requires_support_ownership {
            continue;
        }
        let Some(stolen) = brep.half_edges.get(trim_side.half_edge) else {
            continue;
        };
        let arc_support_side = stolen.twin;
        let Some(support_side) = brep.half_edges.get(arc_support_side) else {
            continue;
        };
        let pocket_loop = support_side.loop_;
        let Some(pocket_face) = brep.loops.get(pocket_loop).map(|l| l.face) else {
            continue;
        };
        if brep
            .faces
            .get(pocket_face)
            .is_none_or(|face| face.role == TopologyRole::ConstructionResidual)
        {
            continue;
        }
        if !pocket_signature_matches(brep, pocket_loop, arc_support_side) {
            continue;
        }
        // The clone rim authored for this trim's source support loop:
        // a same-span boundary side produced by the source-loop repair.
        let Some(clone) = transfers.iter().find_map(|transfer| {
            let candidate = transfer.replacement;
            let side = brep.half_edges.get(candidate)?;
            if half_edge_is_unpaired_boundary(brep, side.twin)
                && half_edge_spans_same_vertices(brep, candidate, arc_support_side)
            {
                Some(candidate)
            } else {
                None
            }
        }) else {
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/spilled-corner] FilletEnd[{}]: no twinless \
                     same-span clone among {} transfer(s) for arc side {arc_support_side:?}",
                    trim_side.fillet_index,
                    transfers.len()
                );
            }
            continue;
        };
        let same_direction = brep
            .half_edges
            .get(arc_support_side)
            .zip(brep.half_edges.get(clone))
            .is_some_and(|(a, b)| a.origin == b.origin);
        if !same_direction {
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/spilled-corner] FilletEnd[{}]: clone {clone:?} \
                     traverses opposite to arc side {arc_support_side:?}; corner left as-is",
                    trim_side.fillet_index
                );
            }
            continue;
        }
        swap_loop_positions(brep, arc_support_side, clone)?;
        if let Some(face) = brep.faces.get_mut(pocket_face) {
            face.role = TopologyRole::ConstructionResidual;
            face.surface = None;
        }
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/spilled-corner] FilletEnd[{}]: arc support side \
                 {arc_support_side:?} took over the support rim from clone {clone:?}; corner \
                 fragment {pocket_face:?} retired",
                trim_side.fillet_index
            );
        }
        normalized += 1;
    }
    Ok(normalized)
}

/// The spilled-corner fragment: a loop of exactly three sides — the arc
/// support side plus two sides meeting at a single interior corner
/// vertex that lies on neither arc endpoint.
fn pocket_signature_matches(
    brep: &BRep,
    pocket_loop: popkern_core::topology::LoopHandle,
    arc_support_side: HalfEdgeHandle,
) -> bool {
    let Ok(cycle) = collect_loop_half_edges(brep, pocket_loop) else {
        return false;
    };
    if cycle.len() != 3 || !cycle.contains(&arc_support_side) {
        return false;
    }
    let arc_vertices: std::collections::BTreeSet<_> = half_edge_vertices(brep, arc_support_side)
        .into_iter()
        .collect();
    if arc_vertices.len() != 2 {
        return false;
    }
    let corner_vertices: std::collections::BTreeSet<_> = cycle
        .iter()
        .filter(|he| **he != arc_support_side)
        .flat_map(|he| half_edge_vertices(brep, *he))
        .filter(|vertex| !arc_vertices.contains(vertex))
        .collect();
    corner_vertices.len() == 1
}

fn trace_decline(clone: HalfEdgeHandle, reason: &str) {
    if crate::direct_apex_trace::enabled() {
        eprintln!("[direct-apex-trace §14.9/spilled-corner] clone {clone:?}: {reason}");
    }
}

/// A boundary side's opposite: either absent from the map outright or
/// present but never wired into a loop.
fn half_edge_is_unpaired_boundary(brep: &BRep, twin: HalfEdgeHandle) -> bool {
    match brep.half_edges.get(twin) {
        None => true,
        Some(twin_side) => !brep.loops.contains_key(twin_side.loop_),
    }
}

fn half_edge_vertices(brep: &BRep, half_edge: HalfEdgeHandle) -> Vec<VertexHandle> {
    let Some(side) = brep.half_edges.get(half_edge) else {
        return Vec::new();
    };
    let mut vertices = vec![side.origin];
    // Derive the end from the loop successor, not the twin: boundary
    // clones carry an unwired twin whose origin is meaningless.
    if let Some(next) = brep.half_edges.get(side.next) {
        vertices.push(next.origin);
    }
    vertices
}

fn half_edge_spans_same_vertices(brep: &BRep, a: HalfEdgeHandle, b: HalfEdgeHandle) -> bool {
    let va: std::collections::BTreeSet<_> = half_edge_vertices(brep, a).into_iter().collect();
    let vb: std::collections::BTreeSet<_> = half_edge_vertices(brep, b).into_iter().collect();
    va.len() == 2 && va == vb
}

/// Exchange the loop positions (prev/next/loop_ and loop anchors) of two
/// half-edges that span the same vertex pair. Euler-neutral.
fn swap_loop_positions(
    brep: &mut BRep,
    first: HalfEdgeHandle,
    second: HalfEdgeHandle,
) -> Result<(), ApexBuildError> {
    let first_data =
        brep.half_edges
            .get(first)
            .cloned()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("spilled-corner normalization: stale half-edge {first:?}"),
            })?;
    let second_data =
        brep.half_edges
            .get(second)
            .cloned()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("spilled-corner normalization: stale half-edge {second:?}"),
            })?;
    for (handle, take_from) in [(first, &second_data), (second, &first_data)] {
        if take_from.prev == handle || take_from.next == handle {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "spilled-corner normalization: half-edges {first:?}/{second:?} are loop \
                     neighbors; swap would corrupt the cycle"
                ),
            });
        }
    }
    let apply = |brep: &mut BRep, handle: HalfEdgeHandle, prev, next, loop_| {
        brep.half_edges[handle].prev = prev;
        brep.half_edges[handle].next = next;
        brep.half_edges[handle].loop_ = loop_;
        brep.half_edges[prev].next = handle;
        brep.half_edges[next].prev = handle;
    };
    apply(
        brep,
        first,
        second_data.prev,
        second_data.next,
        second_data.loop_,
    );
    apply(
        brep,
        second,
        first_data.prev,
        first_data.next,
        first_data.loop_,
    );
    if brep.loops[second_data.loop_].half_edge == second {
        brep.loops[second_data.loop_].half_edge = first;
    }
    if brep.loops[first_data.loop_].half_edge == first {
        brep.loops[first_data.loop_].half_edge = second;
    }
    Ok(())
}
