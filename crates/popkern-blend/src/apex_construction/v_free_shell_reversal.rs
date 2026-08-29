//! Twin-exchange loop reversal for winding reconciliation (PK-153):
//! every loop of the mirror class re-wires to the reversed traversal.
//! Interior seams swap halves between the two adjacent loops (which
//! keeps twins antiparallel by construction), dead sides flip origin
//! in place, and each side carries the reversed pcurve of the side
//! its loop previously held over the same edge — preserving periodic
//! branch choices a fresh projection can wrap.

use std::collections::BTreeSet;

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

/// Reverse every loop of the given faces by twin-exchange. Interior
/// seams swap halves between the two adjacent loops; dead sides flip
/// origin in place. Returns the number of faces reversed.
pub(super) fn reverse_class_loops(brep: &mut BRep, faces: &BTreeSet<FaceHandle>) -> usize {
    struct LoopReversal {
        loop_handle: LoopHandle,
        new_cycle: Vec<HalfEdgeHandle>,
        /// Traversal-aligned pcurves for the new cycle, derived by
        /// reversing the pcurve of the side this loop previously
        /// carried over the same edge — this preserves periodic
        /// branch choices that a fresh projection can wrap.
        new_pcurves: Vec<Option<popkern_core::geometry::pcurve::Pcurve>>,
        flipped_origins: Vec<(HalfEdgeHandle, VertexHandle)>,
    }
    let mut reversals: Vec<LoopReversal> = Vec::new();
    for &face in faces {
        let Some(face_data) = brep.faces.get(face) else {
            continue;
        };
        let mut loops = vec![face_data.outer_loop];
        loops.extend(face_data.inner_loops.iter().copied());
        for loop_handle in loops {
            let Some(loop_data) = brep.loops.get(loop_handle) else {
                continue;
            };
            let mut cycle = Vec::new();
            let start = loop_data.half_edge;
            let mut current = start;
            for _ in 0..brep.half_edges.len() {
                cycle.push(current);
                let Some(data) = brep.half_edges.get(current) else {
                    break;
                };
                current = data.next;
                if current == start {
                    break;
                }
            }
            let ends: Vec<VertexHandle> = cycle
                .iter()
                .map(|side| {
                    let next = brep.half_edges[*side].next;
                    brep.half_edges[next].origin
                })
                .collect();
            let mut new_cycle = Vec::with_capacity(cycle.len());
            let mut new_pcurves = Vec::with_capacity(cycle.len());
            let mut flipped_origins = Vec::new();
            for (side, end) in cycle.iter().zip(&ends).rev() {
                match brep.half_edges.get(brep.half_edges[*side].twin) {
                    Some(_) => new_cycle.push(brep.half_edges[*side].twin),
                    None => {
                        new_cycle.push(*side);
                        flipped_origins.push((*side, *end));
                    }
                }
                new_pcurves.push(
                    brep.half_edges[*side]
                        .pcurve
                        .as_ref()
                        .map(|pcurve| pcurve.reversed()),
                );
            }
            reversals.push(LoopReversal {
                loop_handle,
                new_cycle,
                new_pcurves,
                flipped_origins,
            });
        }
    }
    for reversal in &reversals {
        for (side, origin) in &reversal.flipped_origins {
            brep.half_edges[*side].origin = *origin;
        }
        let n = reversal.new_cycle.len();
        for k in 0..n {
            let side = reversal.new_cycle[k];
            let next = reversal.new_cycle[(k + 1) % n];
            brep.half_edges[side].next = next;
            brep.half_edges[next].prev = side;
            brep.half_edges[side].loop_ = reversal.loop_handle;
            brep.half_edges[side].pcurve = reversal.new_pcurves[k].clone();
        }
        if let Some(first) = reversal.new_cycle.first() {
            brep.loops[reversal.loop_handle].half_edge = *first;
        }
    }
    faces.len()
}
