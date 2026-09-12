//! Pinched keyhole loop normalization (PK-468).
//!
//! A face whose inner loop SHARES A VERTEX with its outer loop is a
//! keyhole the fragment walk left split: at the pinch vertex the walk
//! closed the outer cycle straight through (`in → out`) and emitted
//! the remaining boundary — typically a load-bearing seam slit plus
//! the ring it charts — as a separate inner loop touching the outer
//! loop at that vertex. The canonical form threads ONE boundary cycle
//! through the pinch (measured: the M4 socket-cap tip-chamfer cone
//! serves with its seam slit and tip rim inside the outer loop, while
//! the M5 countersunk tip cone carried the identical geometry as the
//! split form — the mesher's period-window fold then collapsed the
//! slit into a zero-width spur, the face declined to a fallback that
//! produced no triangles, and every neighbour tore along the shared
//! boundary: 159 open mesh edges).
//!
//! The merge is a pure half-edge splice at the shared vertex — no
//! vertex, edge, or face changes, loop count −1 — and is transactional:
//! the merged loop must be one cycle carrying every member of both
//! originals, and structural validation must not regress, else the
//! face keeps its split form. Geometry (curves, pcurves, positions) is
//! untouched, so semantic ownership and stable names are preserved by
//! construction.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

/// Report from the pinched-keyhole merge pass.
#[derive(Debug, Clone, Copy, Default)]
pub struct KeyholeLoopReport {
    /// Inner loops spliced into their face's outer loop.
    pub merged: usize,
    /// Candidates declined (cycle or validation regressed) and rolled
    /// back.
    pub declined: usize,
}

/// Merge inner loops that share a vertex with their face's outer loop.
///
/// Iterates to a fixed point: each merge removes one inner loop, and a
/// face may carry several pinches.
pub fn merge_pinched_keyhole_loops(brep: &mut BRep) -> KeyholeLoopReport {
    let mut report = KeyholeLoopReport::default();
    while let Some((face, inner, he_outer, he_inner)) = find_pinched_pair(brep) {
        if splice_at_pinch(brep, face, inner, he_outer, he_inner) {
            report.merged += 1;
        } else {
            report.declined += 1;
            break;
        }
    }
    report
}

/// Find one (face, inner loop, outer half-edge, inner half-edge) pair
/// whose two half-edges originate at the same vertex. Handle identity,
/// not position: a positional match between distinct vertices is a
/// different class (the near-tangent weld family), never a keyhole.
fn find_pinched_pair(
    brep: &BRep,
) -> Option<(FaceHandle, LoopHandle, HalfEdgeHandle, HalfEdgeHandle)> {
    for (fh, face) in brep.faces.iter() {
        if face.inner_loops.is_empty() {
            continue;
        }
        let outer_members = loop_members(brep, face.outer_loop)?;
        for &inner in &face.inner_loops {
            let inner_members = loop_members(brep, inner)?;
            for &hi in &inner_members {
                let vi = brep.half_edges.get(hi)?.origin;
                for &ho in &outer_members {
                    if brep.half_edges.get(ho)?.origin == vi {
                        return Some((fh, inner, ho, hi));
                    }
                }
            }
        }
    }
    None
}

/// Collect a loop's half-edges by walking `next` from its anchor.
/// `None` when the cycle is broken (defensive: a broken loop is louder
/// upstream residue this pass must not touch).
fn loop_members(brep: &BRep, lh: LoopHandle) -> Option<Vec<HalfEdgeHandle>> {
    let start = brep.loops.get(lh)?.half_edge;
    let mut members = Vec::new();
    let mut cursor = start;
    for _ in 0..brep.half_edges.len() + 1 {
        members.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(members);
        }
    }
    None
}

/// Splice `inner` into the face's outer loop at the shared vertex:
/// `prev(he_outer) → he_inner … prev(he_inner) → he_outer`. Rolls the
/// whole splice back when the merged cycle does not carry every member
/// of both loops or structural validation regresses.
fn splice_at_pinch(
    brep: &mut BRep,
    face: FaceHandle,
    inner: LoopHandle,
    he_outer: HalfEdgeHandle,
    he_inner: HalfEdgeHandle,
) -> bool {
    let snapshot = brep.clone();
    let Some(outer_loop) = brep.faces.get(face).map(|f| f.outer_loop) else {
        return false;
    };
    let (Some(outer_count), Some(inner_members)) = (
        loop_members(brep, outer_loop).map(|m| m.len()),
        loop_members(brep, inner),
    ) else {
        return false;
    };
    let (Some(io), Some(ii)) = (
        brep.half_edges.get(he_outer).map(|h| h.prev),
        brep.half_edges.get(he_inner).map(|h| h.prev),
    ) else {
        return false;
    };
    // The splice: incoming-outer continues into the inner cycle, and
    // incoming-inner continues into the outer cycle.
    brep.half_edges[io].next = he_inner;
    brep.half_edges[he_inner].prev = io;
    brep.half_edges[ii].next = he_outer;
    brep.half_edges[he_outer].prev = ii;
    for h in &inner_members {
        brep.half_edges[*h].loop_ = outer_loop;
    }
    if let Some(f) = brep.faces.get_mut(face) {
        f.inner_loops.retain(|l| *l != inner);
    }
    brep.loops.remove(inner);

    // The merged loop must be one cycle carrying every member.
    let merged_ok = loop_members(brep, outer_loop)
        .is_some_and(|m| m.len() == outer_count + inner_members.len());
    let validation_ok = merged_ok && {
        let before = snapshot.validate_collect_structural().len();
        let after = brep.validate_collect_structural().len();
        after <= before
    };
    if !validation_ok {
        *brep = snapshot;
        return false;
    }
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [keyhole] face {face:?}: inner loop {inner:?} ({} half-edge(s)) spliced into the outer loop at the pinch vertex",
            inner_members.len()
        );
    }
    true
}

#[cfg(test)]
#[path = "keyhole_loop_tests.rs"]
mod tests;
