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
