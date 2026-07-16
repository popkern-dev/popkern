//! Multi-outside-owner (cross-outside-owner) support continuation
//! (PK-092).
//!
//! The single-owner cylinder-mediated machinery requires each
//! generator trim's endpoints to live on ONE outside loop, and the
//! core preflight classifies the repeated-apex trim rewrite's trims
//! as `CrossOutsideOwners` — each trim individually runs from an
//! outside loop to an anchor on the cylinder's own loop (measured on
//! pyramid4_zoo h = 1.2 and the strict V-free gate), so the plan
//! died with the missing-primitive refusal and the main strips were
//! left untrimmed under the partition's wedges.
//!
//! The insight: the PAIRED view is single-owner. Each fillet end
//! carries exactly two cross-owner trims, and their two OUTSIDE
//! endpoints share one loop (the cap loop) while their two ANCHOR
//! endpoints share the cylinder's loop. Pairing them synthesizes a
//! valid single-owner `CylinderMediatedSupportContinuationRequest`
//! (outside span on the cap loop, cylinder interval from the outside
//! adjacency to the anchors' support adjacency), which the EXISTING
//! transactional batch operator executes — including the owned-span
//! exposure phases a raw `face_family_bridge` lacks (its preflight
//! rejects owned boundary sides, measured on this exact fixture).
//!
//! All-or-nothing: one batch preflight + apply for every fillet end,
//! snapshot rollback on any refusal, and the caller's decline path
//! is unchanged when the selection is not in this class.

use std::collections::BTreeMap;

use popkern_core::euler::CylinderMediatedSupportContinuationRequest;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_side_support_plan::{
    OutsideOutsideSupportTrimPlan, PlannedTopologyOperation, SideSupportRegionPlanSelection,
};
use super::repeated_apex_side_support_result::SideSupportRegionContinuationResult;
use super::repeated_apex_topology_knowledge::TopologyKnowledgeGraph;

/// One fillet end's paired cross-owner trims, normalized so that
/// `outside_*` operands live on the shared outside loop and
/// `anchor_*` operands on the shared cylinder loop.
struct PairedCrossOwnerTrims {
    fillet_index: usize,
    cylinder_face: FaceHandle,
    outside_face: FaceHandle,
    outside_loop: LoopHandle,
    /// (vertex, half-edge-at-vertex) for the paired generator ends on
    /// the outside loop, in (start, end) request order.
    outside_start: (VertexHandle, HalfEdgeHandle),
    outside_end: (VertexHandle, HalfEdgeHandle),
    /// Anchor half-edges on the cylinder loop (origins are the
    /// materialized anchors).
    anchor_half_edges: [HalfEdgeHandle; 2],
    cylinder_loop: LoopHandle,
}

/// Attempt the cross-owner continuation for a Blocked selection.
/// `Ok(Some(..))` = committed; `Ok(None)` = not this class;
/// `Err(reason)` = right class, transaction declined (BRep restored).
pub(super) fn apply_cross_owner_continuation(
    selection: &SideSupportRegionPlanSelection,
    knowledge: &TopologyKnowledgeGraph,
    brep: &mut BRep,
) -> Result<Option<SideSupportRegionContinuationResult>, String> {
    let Some(pairs) = paired_cross_owner_trims(selection, brep)? else {
        return Ok(None);
    };
    let mut requests = Vec::with_capacity(pairs.len());
    for pair in &pairs {
        requests.push(paired_request(pair, knowledge, brep)?);
    }
    let rollback = brep.clone();
    let outcome =
        (|| -> Result<SideSupportRegionContinuationResult, String> {
            let batch = brep
                .preflight_cylinder_mediated_support_continuation_batch(requests.clone())
                .map_err(|e| format!("paired batch preflight: {e}"))?;
            let result = brep
                .cylinder_mediated_support_continuation_batch(&batch)
                .map_err(|e| format!("paired batch apply: {e}"))?;
            match result.postconditions {
            popkern_core::euler::CylinderMediatedSupportContinuationPostconditions::Satisfied => {
                Ok(
                    SideSupportRegionContinuationResult::from_outside_outside_cylinder_payload(
                        result.bridge_faces,
                        result.bridge_loops,
                    ),
                )
            }
            popkern_core::euler::CylinderMediatedSupportContinuationPostconditions::Refused {
                reason,
            } => Err(format!("paired batch refused after mutation boundary: {reason}")),
        }
        })();
    match outcome {
        Ok(result) => Ok(Some(result)),
        Err(reason) => {
            *brep = rollback;
            Err(reason)
        }
    }
}

/// Build the single-owner request for one paired fillet end. The
/// cylinder interval runs on the cylinder loop from the outside
/// adjacency (the FilletEnd's cylinder-side half-edge) to the
/// support adjacency at the anchors; the support face is the common
/// twin owner of the anchor half-edges (the far-side support the
/// generators run to).
fn paired_request(
    pair: &PairedCrossOwnerTrims,
    knowledge: &TopologyKnowledgeGraph,
    brep: &BRep,
) -> Result<CylinderMediatedSupportContinuationRequest, String> {
    let twin_owner = |he: HalfEdgeHandle| -> Option<(FaceHandle, LoopHandle)> {
        let twin = brep.half_edges.get(he)?.twin;
        let loop_ = brep.half_edges.get(twin)?.loop_;
        Some((brep.loops.get(loop_)?.face, loop_))
    };
    let (support_a, _) = twin_owner(pair.anchor_half_edges[0]).ok_or_else(|| {
        format!(
            "FilletEnd[{}]: anchor half-edge {:?} has no live twin owner",
            pair.fillet_index, pair.anchor_half_edges[0]
        )
    })?;
    let (support_b, _) = twin_owner(pair.anchor_half_edges[1]).ok_or_else(|| {
        format!(
            "FilletEnd[{}]: anchor half-edge {:?} has no live twin owner",
            pair.fillet_index, pair.anchor_half_edges[1]
        )
    })?;
    if support_a != support_b {
        return Err(format!(
            "FilletEnd[{}]: anchors owned by different supports ({support_a:?} vs {support_b:?})",
            pair.fillet_index
        ));
    }
    // Outside/cylinder adjacency: the half-edge on the OUTSIDE loop
    // whose twin is THIS pair's cylinder, restricted to the arc
    // between the paired endpoints — the cap loop borders every
    // fillet's cylinder, so an unrestricted pick selects a foreign
    // FilletEnd and the core derives an interval sweeping other
    // fillets' territory (measured: exposure rejected on a foreign
    // multi-owner span). The arc that contains the adjacency IS the
    // strip span; if it lies on the complementary arc, the request's
    // endpoint order is swapped to match the core's convention (the
    // selected interval must contain the adjacency).
    let arc_with_adjacency =
        |from: HalfEdgeHandle, to: HalfEdgeHandle| -> Option<(usize, HalfEdgeHandle)> {
            let arc = directed_arc(brep, pair.outside_loop, from, to)?;
            let len = arc.len();
            arc.into_iter()
                .find(|&he| twin_owner(he).map(|(f, _)| f) == Some(pair.cylinder_face))
                .map(|he| (len, he))
        };
    let forward = arc_with_adjacency(pair.outside_start.1, pair.outside_end.1);
    let backward = arc_with_adjacency(pair.outside_end.1, pair.outside_start.1);
    // Both arcs can border this cylinder (the wedge context edges twin
    // to it as well), and picking the long arc makes the core derive a
    // span sweeping every other fillet's territory (measured: a 16-
    // half-edge segment with self-owned cap pairs that the exposure
    // rightly rejects). The strip span is the SHORTER qualifying arc.
    let (outside_start, outside_end, outside_cylinder_half_edge) = match (forward, backward) {
        (Some((lf, hf)), Some((lb, hb))) => {
            if lf <= lb {
                (pair.outside_start, pair.outside_end, hf)
            } else {
                (pair.outside_end, pair.outside_start, hb)
            }
        }
        (Some((_, hf)), None) => (pair.outside_start, pair.outside_end, hf),
        (None, Some((_, hb))) => (pair.outside_end, pair.outside_start, hb),
        (None, None) => {
            return Err(format!(
                "FilletEnd[{}]: neither outside arc between {:?} and {:?} borders cylinder {:?}",
                pair.fillet_index, pair.outside_start.1, pair.outside_end.1, pair.cylinder_face
            ));
        }
    };
    // Support/cylinder adjacency half-edge on the SUPPORT loop whose
    // twin is the cylinder, chosen at one of the anchors: the twin of
    // the anchor half-edge itself.
    let support_cylinder_half_edge = brep
        .half_edges
        .get(pair.anchor_half_edges[0])
        .map(|d| d.twin)
        .ok_or_else(|| format!("FilletEnd[{}]: stale anchor twin", pair.fillet_index))?;
    let support_loop = brep
        .half_edges
        .get(support_cylinder_half_edge)
        .map(|d| d.loop_)
        .ok_or_else(|| format!("FilletEnd[{}]: stale support loop", pair.fillet_index))?;
    // Directed cylinder interval from the OUTSIDE adjacency's
    // cylinder-side half-edge to the anchor.
    let outside_adjacency_on_cylinder = brep
        .half_edges
        .get(outside_cylinder_half_edge)
        .map(|d| d.twin)
        .ok_or_else(|| format!("FilletEnd[{}]: stale outside twin", pair.fillet_index))?;
    let cylinder_interval = directed_arc(
        brep,
        pair.cylinder_loop,
        outside_adjacency_on_cylinder,
        pair.anchor_half_edges[0],
    )
    .or_else(|| {
        directed_arc(
            brep,
            pair.cylinder_loop,
            outside_adjacency_on_cylinder,
            pair.anchor_half_edges[1],
        )
    })
    .ok_or_else(|| {
        format!(
            "FilletEnd[{}]: no directed cylinder interval from {:?} to an anchor",
            pair.fillet_index, outside_adjacency_on_cylinder
        )
    })?;
    Ok(CylinderMediatedSupportContinuationRequest {
        support_face: support_a,
        outside_face: pair.outside_face,
        cylinder_face: pair.cylinder_face,
        start: outside_start.0,
        end: outside_end.0,
        outside_loop: pair.outside_loop,
        outside_start_half_edge: outside_start.1,
        outside_end_half_edge: outside_end.1,
        support_loop,
        support_cylinder_half_edge,
        outside_cylinder_half_edge,
        cylinder_interval,
        context_owner_faces: knowledge.context_faces().to_vec(),
    })
}

/// Group and normalize the plan's cross-owner trims per fillet end,
/// or `None` when the selection is not in this class.
fn paired_cross_owner_trims(
    selection: &SideSupportRegionPlanSelection,
    brep: &BRep,
) -> Result<Option<Vec<PairedCrossOwnerTrims>>, String> {
    let plan = match selection {
        SideSupportRegionPlanSelection::Selected(plan)
        | SideSupportRegionPlanSelection::Blocked(plan, _) => plan,
        SideSupportRegionPlanSelection::Refused(_) => return Ok(None),
    };
    let mut trims_by_fillet: BTreeMap<usize, Vec<&OutsideOutsideSupportTrimPlan>> = BTreeMap::new();
    for operation in plan.operations() {
        if let PlannedTopologyOperation::OutsideOutsideSupportContinuation(region) = operation {
            for trim in region.trims() {
                trims_by_fillet
                    .entry(trim.fillet_index)
                    .or_default()
                    .push(trim);
            }
        }
    }
    if trims_by_fillet.is_empty() {
        return Ok(None);
    }
    let mut pairs = Vec::with_capacity(trims_by_fillet.len());
    for (fillet_index, trims) in &trims_by_fillet {
        let [ta, tb] = trims.as_slice() else {
            return Err(format!(
                "cross-owner pairing needs exactly 2 trims for FilletEnd[{fillet_index}]; got {}",
                trims.len()
            ));
        };
        if ta.start_loop == ta.end_loop && tb.start_loop == tb.end_loop {
            // Single-owner trims — not this class.
            return Ok(None);
        }
        if ta.cylinder_face != tb.cylinder_face {
            return Err(format!(
                "cross-owner FilletEnd[{fillet_index}]: trims name different cylinders \
                 ({:?} vs {:?})",
                ta.cylinder_face, tb.cylinder_face
            ));
        }
        // Normalize: ta runs outside->anchor, tb anchor->outside; the
        // outside side is whichever shared loop is NOT the cylinder's.
        if !(ta.start_loop == tb.end_loop && ta.end_loop == tb.start_loop) {
            return Err(format!(
                "cross-owner FilletEnd[{fillet_index}] does not straddle a shared loop pair: \
                 ta {:?}->{:?}, tb {:?}->{:?}",
                ta.start_loop, ta.end_loop, tb.start_loop, tb.end_loop
            ));
        }
        let cylinder_owns = |l: LoopHandle| {
            brep.loops
                .get(l)
                .is_some_and(|d| d.face == ta.cylinder_face)
        };
        let (outside_face, outside_loop, outside_start, outside_end, anchors, cylinder_loop) =
            if cylinder_owns(ta.end_loop) {
                (
                    ta.start_face,
                    ta.start_loop,
                    (ta.start, ta.start_half_edge),
                    (tb.end, tb.end_half_edge),
                    [ta.end_half_edge, tb.start_half_edge],
                    ta.end_loop,
                )
            } else if cylinder_owns(ta.start_loop) {
                (
                    ta.end_face,
                    ta.end_loop,
                    (tb.start, tb.start_half_edge),
                    (ta.end, ta.end_half_edge),
                    [ta.start_half_edge, tb.end_half_edge],
                    ta.start_loop,
                )
            } else {
                return Err(format!(
                    "cross-owner FilletEnd[{fillet_index}]: neither shared loop belongs to the \
                     cylinder {:?}",
                    ta.cylinder_face
                ));
            };
        pairs.push(PairedCrossOwnerTrims {
            fillet_index: *fillet_index,
            cylinder_face: ta.cylinder_face,
            outside_face,
            outside_loop,
            outside_start,
            outside_end,
            anchor_half_edges: anchors,
            cylinder_loop,
        });
    }
    Ok(Some(pairs))
}

fn directed_arc(
    brep: &BRep,
    loop_handle: LoopHandle,
    from: HalfEdgeHandle,
    to: HalfEdgeHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let mut cursor = from;
    let mut arc = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        arc.push(cursor);
        if cursor == to {
            return Some(arc);
        }
        cursor = half_edge.next;
        if cursor == from {
            break;
        }
    }
    None
}
