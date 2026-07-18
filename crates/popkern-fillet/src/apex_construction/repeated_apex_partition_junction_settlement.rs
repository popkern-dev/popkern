//! Central pocket-junction settlement for planned partition cycles
//! (PK-129). Extracted from the topology module to hold it under the
//! 500-line budget; see [`settle_planned_loop_junction_gaps`].

use popkern_core::geometry::CurveEval;
use popkern_core::topology::{BRep, HalfEdgeHandle, TopologyRole, VertexHandle};

/// Settle pocket-junction gaps in a PLANNED half-edge cycle (the
/// PK-129 design fork). At N = 3 the trim-rewrite execution order
/// leaves live split/pocket half-edges BETWEEN consecutive planned
/// sides, so a planned cycle is not contiguous at those junctions.
/// Instead of each consumer bridging its own cycle (the cap boundary
/// did; the strip loop's mirror gap then failed one stage later),
/// every planned cycle is settled HERE, before wiring: junction gaps
/// are spliced with the live connecting chain found by bounded
/// vertex-chain search (deterministic smallest-handle choice;
/// `.next` chains are not authoritative on freshly cloned sides).
/// Cycles without gaps (all N >= 4 plans) pass through unchanged.
pub(crate) fn settle_planned_loop_junction_gaps(
    brep: &BRep,
    cycle: &[HalfEdgeHandle],
) -> Vec<HalfEdgeHandle> {
    const MAX_BRIDGE: usize = 8;
    let mut settled = Vec::with_capacity(cycle.len());
    for (index, &current) in cycle.iter().enumerate() {
        settled.push(current);
        let next = cycle[(index + 1) % cycle.len()];
        let Some(current_end) = half_edge_destination_vertex(brep, current) else {
            continue;
        };
        let Some(next_start) = brep.half_edges.get(next).map(|half_edge| half_edge.origin) else {
            continue;
        };
        if current_end == next_start {
            continue;
        }
        let mut cursor_vertex = current_end;
        let mut chain = Vec::new();
        'bridge: for _ in 0..MAX_BRIDGE {
            // Candidate rank: connecting material must not be stolen
            // from a LIVE production loop (measured: the greedy
            // smallest-handle choice took the target cylinder strip's
            // own seam side into the cap cycle, and the strip loop
            // then failed to close). Rank 0 = unwired fresh clones
            // (the reauthored FilletEnd pieces), rank 1 = residual
            // scaffolding to be consumed; live production sides are
            // excluded.
            let claim_rank = |handle: HalfEdgeHandle| -> Option<u8> {
                let data = brep.half_edges.get(handle)?;
                match brep.loops.get(data.loop_) {
                    None => Some(0),
                    Some(loop_data) => match brep.faces.get(loop_data.face) {
                        None => Some(0),
                        Some(face_data) => {
                            if face_data.role == TopologyRole::ConstructionResidual {
                                Some(1)
                            } else {
                                None
                            }
                        }
                    },
                }
            };
            let mut candidates: Vec<(u8, HalfEdgeHandle, VertexHandle)> = brep
                .half_edges
                .iter()
                .filter(|(handle, data)| {
                    data.origin == cursor_vertex
                        && *handle != next
                        && !cycle.contains(handle)
                        && !chain.contains(handle)
                })
                .filter_map(|(handle, _)| {
                    let rank = claim_rank(handle)?;
                    half_edge_destination_vertex(brep, handle).map(|end| (rank, handle, end))
                })
                .collect();
            candidates.sort_by_key(|&(rank, handle, _)| (rank, handle));
            // Prefer a candidate landing exactly on the next start.
            if let Some(&(_, handle, _)) = candidates.iter().find(|(_, _, end)| *end == next_start)
            {
                chain.push(handle);
                if crate::direct_apex_trace::enabled() {
                    let pos = |v: VertexHandle| {
                        brep.vertices
                            .get(v)
                            .map(|d| {
                                format!(
                                    "({:.4},{:.4},{:.4})",
                                    d.position.x, d.position.y, d.position.z
                                )
                            })
                            .unwrap_or_else(|| "stale".into())
                    };
                    let extent: f64 = chain
                        .iter()
                        .filter_map(|&h| {
                            let d = brep.half_edges.get(h)?;
                            let a = brep.vertices.get(d.origin)?.position;
                            let e = half_edge_destination_vertex(brep, h)?;
                            let b = brep.vertices.get(e)?.position;
                            Some(a.distance_to(&b))
                        })
                        .sum();
                    let loop_info = |h: HalfEdgeHandle| {
                        let Some(d) = brep.half_edges.get(h) else {
                            return "stale".to_string();
                        };
                        let Some(l) = brep.loops.get(d.loop_) else {
                            return format!("{:?}:no-loop", d.loop_);
                        };
                        let face_role = brep
                            .faces
                            .get(l.face)
                            .map(|f| format!("{:?}", f.role))
                            .unwrap_or_else(|| "stale-face".into());
                        format!("{:?} face {:?} role {face_role}", d.loop_, l.face)
                    };
                    eprintln!(
                        "[junction-settle] gap {current:?}->{next:?}: {} chain hes, extent {extent:.3e}, {} -> {}, donors [{}]",
                        chain.len(),
                        pos(current_end),
                        pos(next_start),
                        chain
                            .iter()
                            .map(|&h| loop_info(h))
                            .collect::<Vec<_>>()
                            .join("; ")
                    );
                }
                settled.append(&mut chain);
                break 'bridge;
            }
            let Some(&(_, handle, end)) = candidates.first() else {
                break;
            };
            chain.push(handle);
            cursor_vertex = end;
        }
    }
    settled
}

pub(in super::super) fn half_edge_destination_vertex(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<VertexHandle> {
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
