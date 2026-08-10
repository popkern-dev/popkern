//! Support-source loop repair and coplanar FaceEdge routing for the
//! repeated-apex support-loop rewrite: replaces stolen cylinder sides
//! with boundary spans, classifies FaceEdge sides onto the support
//! group whose facet plane carries both endpoints, and splices routed
//! sides out of their donor loops.

use super::super::*;

use super::{author_support_source_loop_replacement_span, support_loop_contains_half_edge};

pub(in crate::apex_construction) fn repair_support_source_loops_after_cylinder_rewire(
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

pub(in crate::apex_construction) struct SupportSourceLoopTransfer {
    pub(in crate::apex_construction) loop_: LoopHandle,
    pub(in crate::apex_construction) replacement: HalfEdgeHandle,
}

/// The bridge side a support group can actually consume: the one
/// whose direction chains with the group's already-fixed sides (the
/// generator directions are twins of the wired cylinder loops, so
/// they encode the shell winding). At mirrored corners the final cap
/// winding is the REVERSE of the §14.7 mef authoring parity
/// (reauthored FilletEnd seams twin the cylinder strips' retained
/// windings), which puts the chainable direction on the recorded CAP
/// side of the bridge — take it only when the recorded support side
/// cannot chain and the twin chains at both endpoints.
pub(super) fn winding_honest_bridge_side(
    brep: &BRep,
    recorded: HalfEdgeHandle,
    group: &[HalfEdgeHandle],
) -> HalfEdgeHandle {
    let chains = |side: HalfEdgeHandle| -> bool {
        let (Some(start), Some(end)) = (
            brep.half_edges.get(side).map(|data| data.origin),
            support_split_half_edge_end_vertex(brep, side),
        ) else {
            return false;
        };
        let feeds_start = group
            .iter()
            .any(|other| support_split_half_edge_end_vertex(brep, *other) == Some(start));
        let continues_end = group.iter().any(|other| {
            brep.half_edges
                .get(*other)
                .is_some_and(|data| data.origin == end)
        });
        feeds_start && continues_end
    };
    let twin = brep.half_edges.get(recorded).map(|data| data.twin);
    match twin {
        Some(twin) if !chains(recorded) && chains(twin) => {
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/bridge-side] recorded support side {recorded:?} \
                     opposes the group winding; taking cap-side twin {twin:?}"
                );
            }
            twin
        }
        _ => recorded,
    }
}

/// Exchange bridge-edge side ownership between the retiring cap chain
/// and the support group: the group takes the cap-owned side, and the
/// recorded support side (spliced out of its retiring donor loop)
/// takes the taken side's chain slot so later cap-boundary scans find
/// the honest-direction side on the cap loop. The taken side is left
/// unwired until the support cycle claims it.
pub(super) fn exchange_bridge_side_with_cap_chain(
    brep: &mut BRep,
    taken_for_support: HalfEdgeHandle,
    recorded_support_side: HalfEdgeHandle,
) {
    let Some(slot) = brep.half_edges.get(taken_for_support).cloned() else {
        return;
    };
    splice_side_out_of_donor_loop(brep, recorded_support_side);
    if slot.prev != taken_for_support
        && slot.next != taken_for_support
        && brep.half_edges.contains_key(slot.prev)
        && brep.half_edges.contains_key(slot.next)
    {
        brep.half_edges[slot.prev].next = recorded_support_side;
        brep.half_edges[slot.next].prev = recorded_support_side;
        brep.half_edges[recorded_support_side].prev = slot.prev;
        brep.half_edges[recorded_support_side].next = slot.next;
    }
    brep.half_edges[recorded_support_side].loop_ = slot.loop_;
    if brep
        .loops
        .get(slot.loop_)
        .is_some_and(|loop_data| loop_data.half_edge == taken_for_support)
    {
        brep.loops[slot.loop_].half_edge = recorded_support_side;
    }
    brep.half_edges[taken_for_support].loop_ = LoopHandle::default();
}

/// Remove a side from its current loop's chain before another loop
/// adopts it: predecessor and successor join, and the loop re-anchors
/// if the side was its anchor. The donor is retiring scaffolding.
pub(super) fn splice_side_out_of_donor_loop(brep: &mut BRep, side: HalfEdgeHandle) {
    let Some(data) = brep.half_edges.get(side).cloned() else {
        return;
    };
    if data.prev != side
        && data.next != side
        && brep.half_edges.contains_key(data.prev)
        && brep.half_edges.contains_key(data.next)
    {
        brep.half_edges[data.prev].next = data.next;
        brep.half_edges[data.next].prev = data.prev;
    }
    if brep
        .loops
        .get(data.loop_)
        .is_some_and(|loop_data| loop_data.half_edge == side)
        && data.next != side
    {
        brep.loops[data.loop_].half_edge = data.next;
    }
}

pub(super) fn face_edge_endpoints_fit_face_plane(
    brep: &BRep,
    face: FaceHandle,
    side: HalfEdgeHandle,
) -> bool {
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
    let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    if !matches!(surface.geometry, SurfaceGeometry::Plane { .. }) {
        return false;
    }
    let Some(data) = brep.half_edges.get(side) else {
        return false;
    };
    let Some(end) = brep
        .half_edges
        .get(data.twin)
        .map(|twin| twin.origin)
        .or_else(|| brep.half_edges.get(data.next).map(|next| next.origin))
    else {
        return false;
    };
    let fit = brep.tolerance.linear().max(1.0e-9) * 10.0;
    [data.origin, end].into_iter().all(|vertex| {
        brep.vertices.get(vertex).is_some_and(|v| {
            surface
                .geometry
                .closest_point(&v.position)
                .map(|(u, w)| surface.geometry.point_at(u, w).distance_to(&v.position) <= fit)
                .unwrap_or(false)
        })
    })
}
