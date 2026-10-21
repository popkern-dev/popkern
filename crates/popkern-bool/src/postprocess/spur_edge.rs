//! Post-merge spur collapse: remove degenerate spur edges and
//! degenerate self-twinned stub half-edges.
//!
//! Two related residue classes from the merge/walk phase (PK-293,
//! measured on a rim-filleted cylinder minus a cuboid):
//!
//! 1. **MEV spurs** — a "there-and-back" half-edge pair whose tip
//!    vertex has valence 1. Retired transactionally with the `kev`
//!    Euler operator.
//! 2. **Self-twinned stubs** — the assembly's dangling-half-edge
//!    cleanup points a half-wired survivor's twin at ITSELF (see
//!    `merge/assemble.rs`); when later dissolution passes shrink that
//!    survivor's edge to a single position, the result carries a
//!    zero-length self-twinned half-edge wired into a real loop, and
//!    `validate_all` reports its vertex as `LowValence`. `kev`
//!    cannot retire these (its preconditions assume healthy
//!    pairing), so they are excised directly — the stub spans one
//!    position, making the loop stitch a geometric no-op — under a
//!    snapshot/rollback guard like [`super::loop_repair`].
//!
//! Like `loop_repair`, this is detection-side repair: the walk-level
//! half-wiring that strands these halves is the recorded follow-up.
//! Only degenerate (≤ 10× linear tolerance) residue is touched; a
//! long spur is a louder upstream bug and stays visible.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::serve_region::ServeRegion;

/// Report from a spur-collapse pass.
#[derive(Debug, Clone, Copy, Default)]
pub struct SpurCollapseReport {
    /// MEV spur edges removed (with their tip vertices) via `kev`.
    pub spurs_removed: usize,
    /// Degenerate self-twinned stub half-edges excised.
    pub self_twins_removed: usize,
    /// Ghost faces (every boundary half-edge a degenerate self-twin)
    /// retired whole.
    pub ghost_faces_removed: usize,
}

/// Remove degenerate spur edges, ghost faces, and self-twinned
/// stubs. Iterates to a fixed point so chains retire fully.
pub fn collapse_degenerate_spur_edges(
    brep: &mut BRep,
    tolerance: &Tolerance,
) -> SpurCollapseReport {
    collapse_degenerate_spur_edges_in(brep, tolerance, &ServeRegion::Whole)
}

/// [`collapse_degenerate_spur_edges`] on the faces `region` admits: a
/// ghost face or a stub's face in it, a spur edge with both sides in it
/// (removals create nothing, so the region stays exact).
pub(crate) fn collapse_degenerate_spur_edges_in(
    brep: &mut BRep,
    tolerance: &Tolerance,
    region: &ServeRegion,
) -> SpurCollapseReport {
    let mut report = SpurCollapseReport::default();
    let length_gate = tolerance.linear() * 10.0;
    loop {
        if retire_one_ghost_face(brep, length_gate, region) {
            report.ghost_faces_removed += 1;
            continue;
        }
        if retire_one_mev_spur(brep, length_gate, region) {
            report.spurs_removed += 1;
            continue;
        }
        if excise_one_self_twin_stub(brep, length_gate, region) {
            report.self_twins_removed += 1;
            continue;
        }
        break;
    }
    report
}

/// True when the half-edge is a degenerate self-twinned stub: its
/// twin is itself and its edge has no curve (the assembly cleanup's
/// signature) or a curve spanning less than `length_gate`.
fn is_degenerate_self_twin(
    brep: &BRep,
    heh: popkern_core::topology::HalfEdgeHandle,
    length_gate: f64,
) -> bool {
    let Some(he) = brep.half_edges.get(heh) else {
        return false;
    };
    if he.twin != heh {
        return false;
    }
    let Some(edge) = brep.edges.get(he.edge) else {
        return true; // dangling edge reference: pure residue
    };
    match edge.curve.as_ref() {
        None => true,
        Some(c) => {
            use popkern_core::geometry::curve::CurveEval;
            c.geometry
                .point_at(c.t_start)
                .distance_to(&c.geometry.point_at(c.t_end))
                <= length_gate
        }
    }
}

/// Retire ONE ghost face — a face whose EVERY loop half-edge is a
/// degenerate self-twinned stub (walk residue with no real
/// boundary). Removes the face, its loops, half-edges, edges, and
/// orphaned vertices, with snapshot rollback if structural
/// validation regresses.
fn retire_one_ghost_face(brep: &mut BRep, length_gate: f64, region: &ServeRegion) -> bool {
    let target = brep.faces.iter().find_map(|(fh, f)| {
        if !region.admits_face(fh) {
            return None;
        }
        let mut loop_handles = vec![f.outer_loop];
        loop_handles.extend(f.inner_loops.iter().copied());
        // Walked with an early exit: a real face fails on its first
        // half-edge (PK-815: the 64-hole caps were walked whole).
        let mut hes = Vec::new();
        for lh in &loop_handles {
            let l = brep.loops.get(*lh)?;
            let start = l.half_edge;
            let mut cur = start;
            let mut steps = 0usize;
            loop {
                if !is_degenerate_self_twin(brep, cur, length_gate) {
                    return None;
                }
                hes.push(cur);
                cur = brep.half_edges.get(cur)?.next;
                steps += 1;
                if cur == start {
                    break;
                }
                if steps > 100_000 {
                    return None;
                }
            }
        }
        (!hes.is_empty()).then_some((fh, loop_handles, hes))
    });
    let Some((face_h, loop_handles, hes)) = target else {
        return false;
    };
    let snapshot = brep.clone();
    let mut edge_handles = Vec::new();
    let mut vertex_handles = Vec::new();
    for &heh in &hes {
        if let Some(he) = brep.half_edges.get(heh) {
            edge_handles.push(he.edge);
            vertex_handles.push(he.origin);
        }
        brep.remove_half_edge_reanchored(heh);
    }
    for eh in edge_handles {
        if !brep.half_edges.iter().any(|(_, he)| he.edge == eh) {
            brep.edges.remove(eh);
        }
    }
    for lh in loop_handles {
        brep.loops.remove(lh);
    }
    brep.faces.remove(face_h);
    for vh in vertex_handles {
        match brep
            .half_edges
            .iter()
            .find(|(_, he)| he.origin == vh)
            .map(|(heh, _)| heh)
        {
            Some(out_he) => {
                if let Some(v) = brep.vertices.get_mut(vh) {
                    v.outgoing = out_he;
                }
            }
            None => {
                brep.vertices.remove(vh);
            }
        }
    }
    let before = snapshot.validate_collect_structural().len();
    let after = brep.validate_collect_structural().len();
    if after > before {
        if popkern_core::trace::trace_enabled("spur") {
            eprintln!(
                "[spur] ghost face {face_h:?}: rollback — structural issues {before} -> {after}"
            );
        }
        *brep = snapshot;
        return false;
    }
    true
}

/// Retire ONE short MEV spur (tip valence 1) via the transactional
/// `kev` operator. Returns whether anything was removed.
fn retire_one_mev_spur(brep: &mut BRep, length_gate: f64, region: &ServeRegion) -> bool {
    // The short edges first (each as its half-edge and twin): the
    // valence census over every half-edge is paid only when one exists
    // (PK-815: every boolean's result).
    let short: Vec<_> = brep
        .edges
        .iter()
        .filter_map(|(eh, d)| {
            let he = brep.half_edges.get(d.half_edge)?;
            if he.twin == d.half_edge {
                return None; // self-twin: handled by the stub pass
            }
            let twin = brep.half_edges.get(he.twin)?;
            let a = brep.vertices.get(he.origin)?.position;
            let b = brep.vertices.get(twin.origin)?.position;
            (a.distance_to(&b) <= length_gate && region.admits_edge(brep, eh))
                .then_some((d.half_edge, he.twin))
        })
        .collect();
    if short.is_empty() {
        return false;
    }
    // Dense on the vertex handles.
    let mut valence: slotmap::SecondaryMap<popkern_core::topology::VertexHandle, usize> =
        slotmap::SecondaryMap::with_capacity(brep.vertices.len());
    for (_, he) in brep.half_edges.iter() {
        if let Some(n) = valence.get_mut(he.origin) {
            *n += 1;
        } else {
            valence.insert(he.origin, 1);
        }
    }
    let origin_valence = |h| {
        brep.half_edges
            .get(h)
            .and_then(|he| valence.get(he.origin).copied())
            .unwrap_or(0)
    };
    // `kev` removes the edge plus the TWIN's origin (the tip); pick
    // whichever half-edge points INTO a valence-1 tip.
    let candidates: Vec<_> = short
        .into_iter()
        .filter_map(|(he, twin)| {
            if origin_valence(twin) == 1 {
                Some(he)
            } else if origin_valence(he) == 1 {
                Some(twin)
            } else {
                None
            }
        })
        .collect();
    for he_into_tip in candidates {
        if brep.half_edges.contains_key(he_into_tip) && brep.kev(he_into_tip).is_ok() {
            return true;
        }
    }
    false
}

/// Excise ONE degenerate self-twinned stub half-edge from its loop
/// chain, with snapshot rollback if structural validation regresses.
/// Returns whether anything was removed.
fn excise_one_self_twin_stub(brep: &mut BRep, length_gate: f64, region: &ServeRegion) -> bool {
    let candidates: Vec<_> = brep
        .half_edges
        .iter()
        .filter_map(|(heh, he)| {
            if he.twin != heh || !region.admits_face(brep.loops.get(he.loop_)?.face) {
                return None;
            }
            // Degenerate: the stub's edge spans (numerically) one
            // point AND its loop-chain hop is a no-op.
            let edge = brep.edges.get(he.edge)?;
            let span = edge
                .curve
                .as_ref()
                .map(|c| {
                    use popkern_core::geometry::curve::CurveEval;
                    c.geometry
                        .point_at(c.t_start)
                        .distance_to(&c.geometry.point_at(c.t_end))
                })
                .unwrap_or(0.0);
            let origin = brep.vertices.get(he.origin)?.position;
            let next_origin = brep
                .half_edges
                .get(he.next)
                .and_then(|n| brep.vertices.get(n.origin))
                .map(|v| v.position)?;
            if span > length_gate || origin.distance_to(&next_origin) > length_gate {
                return None;
            }
            Some(heh)
        })
        .collect();
    for target in candidates {
        if brep.half_edges.contains_key(target) && excise_stub(brep, target) {
            return true;
        }
    }
    false
}

/// Remove `stub_h` (a self-twinned half-edge) from its loop and the
/// model. Rolls back on structural validation regression.
fn excise_stub(brep: &mut BRep, stub_h: popkern_core::topology::HalfEdgeHandle) -> bool {
    let snapshot = brep.clone();
    let stub = &brep.half_edges[stub_h];
    let loop_h = stub.loop_;
    let next_h = stub.next;
    let edge_h = stub.edge;
    let tip_v = stub.origin;

    if next_h == stub_h {
        // The stub is its entire loop: retire the loop when it is an
        // inner loop of its face; a degenerate OUTER loop is a louder
        // bug and stays visible.
        let Some(loop_data) = brep.loops.get(loop_h) else {
            return false;
        };
        let face_h = loop_data.face;
        let Some(face) = brep.faces.get_mut(face_h) else {
            return false;
        };
        if face.outer_loop == loop_h {
            if popkern_core::trace::trace_enabled("spur") {
                eprintln!("[spur] {stub_h:?}: 1-element OUTER loop — left visible");
            }
            return false;
        }
        face.inner_loops.retain(|&l| l != loop_h);
        brep.loops.remove(loop_h);
        brep.remove_half_edge_reanchored(stub_h);
    } else {
        // Stitch the loop chain past the stub.
        let mut prev_h = next_h;
        let mut steps = 0usize;
        loop {
            let Some(prev) = brep.half_edges.get(prev_h) else {
                return false;
            };
            if prev.next == stub_h {
                break;
            }
            prev_h = prev.next;
            steps += 1;
            if steps > 100_000 {
                return false; // unterminated chain — leave visible
            }
        }
        brep.half_edges[prev_h].next = next_h;
        brep.half_edges[next_h].prev = prev_h;
        if let Some(loop_data) = brep.loops.get_mut(loop_h)
            && loop_data.half_edge == stub_h
        {
            loop_data.half_edge = next_h;
        }
        brep.remove_half_edge_reanchored(stub_h);
    }

    // Retire the edge if no other half-edge references it, and the
    // vertex if nothing originates there any more.
    if !brep.half_edges.iter().any(|(_, he)| he.edge == edge_h) {
        brep.edges.remove(edge_h);
    }
    match brep
        .half_edges
        .iter()
        .find(|(_, he)| he.origin == tip_v)
        .map(|(heh, _)| heh)
    {
        Some(out_he) => {
            if let Some(v) = brep.vertices.get_mut(tip_v) {
                v.outgoing = out_he;
            }
        }
        None => {
            brep.vertices.remove(tip_v);
        }
    }

    // Structural no-regression guard: never trade one defect class
    // for another.
    let before = snapshot.validate_collect_structural().len();
    let after = brep.validate_collect_structural().len();
    if after > before {
        if popkern_core::trace::trace_enabled("spur") {
            eprintln!(
                "[spur] {stub_h:?}: rollback — structural issues {before} -> {after}: {:?}",
                brep.validate_collect_structural().first()
            );
        }
        *brep = snapshot;
        return false;
    }
    true
}

#[cfg(test)]
#[path = "spur_edge_tests.rs"]
mod tests;
