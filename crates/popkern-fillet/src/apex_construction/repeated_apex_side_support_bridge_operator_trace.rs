//! PK-167 stage 1b observability: dump the face-family bridge cycle
//! structure at execution time (post-exposure) under
//! `POPKERN_TRACE=support-loop-rewrite`. The in-place support-loop rewrite derivation
//! needs the measured near/far classification of the source loop's
//! sides — this trace is how that structure was captured on the
//! V-free mixed-cube fixture (each shared-start source is a 3-sided
//! corner fragment: two dead sides where the corner material was
//! removed plus one rail side already twinning its strip).

use popkern_core::topology::BRep;

use super::super::repeated_apex_side_support_bridge::FaceFamilyBridgePrimitivePlan;

pub(super) fn trace_cycle_structure(primitive: &FaceFamilyBridgePrimitivePlan, candidate: &BRep) {
    if !popkern_core::trace::trace_enabled("support-loop-rewrite") {
        return;
    }
    for region in primitive.regions() {
        for (ci, cycle) in region.boundary_cycles().iter().enumerate() {
            for (label, sides) in [
                ("support", cycle.support_interval()),
                ("outside", cycle.outside_interval()),
            ] {
                for &h in sides {
                    let Some(data) = candidate.half_edges.get(h) else {
                        eprintln!("[slir] c{ci} {label} {h:?} STALE");
                        continue;
                    };
                    let Some(v) = candidate.vertices.get(data.origin) else {
                        eprintln!("[slir] c{ci} {label} {h:?} stale-origin");
                        continue;
                    };
                    let p = v.position;
                    let twin = candidate.half_edges.get(data.twin).map(|t| {
                        let f = candidate.loops.get(t.loop_).map(|l| l.face);
                        let kind = f.and_then(|f| candidate.faces.get(f)).map(|fd| {
                            (
                                fd.role,
                                fd.surface
                                    .as_ref()
                                    .map(|s| std::mem::discriminant(&s.geometry)),
                            )
                        });
                        format!("{f:?} {kind:?}")
                    });
                    eprintln!(
                        "[slir] c{ci} {label} {h:?} o=({:.3},{:.3},{:.3}) twin {}",
                        p.x,
                        p.y,
                        p.z,
                        twin.unwrap_or_else(|| "DEAD".into())
                    );
                }
            }
            trace_gap_chain(candidate, cycle, ci);
        }
    }
}

/// Walk the gap's outer boundary by VERTEX IDENTITY between the
/// removed interval's endpoints — the open verification for the
/// in-place rewrite plan: the chain must close on shared vertex
/// handles (dead sides chained destination->origin), or the rewrite
/// needs a pre-weld/connector step first.
fn trace_gap_chain(
    candidate: &BRep,
    cycle: &super::super::repeated_apex_side_support_bridge::FaceFamilyBridgeBoundaryCycle,
    ci: usize,
) {
    let support = cycle.support_interval();
    if support.len() < 2 {
        return;
    }
    let removed = &support[..support.len() - 1];
    let Some(first) = candidate.half_edges.get(removed[0]) else {
        return;
    };
    let start_vertex = first.origin;
    let Some(&last) = removed.last() else { return };
    let end_vertex = candidate
        .half_edges
        .get(last)
        .and_then(|h| candidate.half_edges.get(h.next))
        .map(|n| n.origin);
    let Some(end_vertex) = end_vertex else { return };
    // Dead sides indexed by destination vertex.
    let dead: Vec<_> = candidate
        .half_edges
        .iter()
        .filter(|(h, data)| !candidate.half_edges.contains_key(data.twin) && !removed.contains(h))
        .map(|(h, _)| h)
        .collect();
    let dest_of = |h| {
        candidate
            .half_edges
            .get(h)
            .and_then(|d: &popkern_core::topology::HalfEdge| candidate.half_edges.get(d.next))
            .map(|n| n.origin)
    };
    let mut chain = Vec::new();
    let mut want = start_vertex;
    for _ in 0..16 {
        let Some(&h) = dead
            .iter()
            .find(|&&h| dest_of(h) == Some(want) && !chain.contains(&h))
        else {
            break;
        };
        chain.push(h);
        want = candidate.half_edges[h].origin;
        if want == end_vertex {
            break;
        }
    }
    let closed = want == end_vertex;
    eprintln!(
        "[slir] c{ci} gap chain: {} sides, vertex-identity {} (start {start_vertex:?} end {end_vertex:?})",
        chain.len(),
        if closed { "CLOSES" } else { "BROKEN" }
    );
    for &h in &chain {
        let d = &candidate.half_edges[h];
        let p = candidate.vertices[d.origin].position;
        eprintln!(
            "[slir]   gap side {h:?} o=({:.3},{:.3},{:.3})",
            p.x, p.y, p.z
        );
    }
}
