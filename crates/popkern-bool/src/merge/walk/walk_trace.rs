//! Env-gated election tracing (the loop-dump trace topic / the loop-dump-all trace topic),
//! extracted under the file-size ratchet.

/// Per-face pre-walk dump (the `fragment-walk` trace topic): the face's
/// half-edge population and each half-edge's endpoints, twin, source, and
/// curve parameter range. Extracted from the walk body under the
/// file-size ratchet.
pub(super) fn dump_face_pre(
    graph: &super::super::graph::IntersectionGraph,
    face_ref: &super::super::graph::FaceRef,
    he_indices: &[usize],
) {
    if !popkern_core::trace::trace_enabled("fragment-walk") {
        return;
    }
    let n_int = he_indices
        .iter()
        .filter(|&&h| {
            matches!(
                graph.half_edges[h].source,
                super::super::graph::EdgeSource::Shared { .. }
            )
        })
        .count();
    eprintln!(
        "  [walk-pre] {face_ref:?}: {} hes ({n_int} intersection)",
        he_indices.len()
    );
    for &h in he_indices.iter() {
        let he = &graph.half_edges[h];
        let o = he.origin;
        let e = graph.half_edges[he.twin].origin;
        let po = graph.nodes[o].position;
        let pe = graph.nodes[e].position;
        eprintln!(
            "      [he] {h} n{o}({:.9},{:.9},{:.9}) -> n{e}({:.9},{:.9},{:.9}) twin {} src {:?} curve {}",
            po.x,
            po.y,
            po.z,
            pe.x,
            pe.y,
            pe.z,
            he.twin,
            he.source,
            he.curve
                .as_ref()
                .or(graph.half_edges[he.twin].curve.as_ref())
                .map(|c| format!("[{:.2},{:.2}]", c.t_start, c.t_end))
                .unwrap_or_else(|| "none".into())
        );
    }
}

/// the loop-dump trace topic / the loop-dump-all trace topic election tracing (extracted from
/// the walk body under the file-size ratchet).
pub(super) fn dump_election_loops(
    graph: &super::super::graph::IntersectionGraph,
    face_ref: &super::super::graph::FaceRef,
    outer_loop: &[usize],
    inner_loops_out: &[Vec<usize>],
    separate_outers: &[Vec<usize>],
) {
    if popkern_core::trace::trace_enabled("loop-dump-all") {
        for (tag, lp) in std::iter::once(("outer", outer_loop))
            .chain(inner_loops_out.iter().map(|l| ("inner", l.as_slice())))
            .chain(separate_outers.iter().map(|l| ("separate", l.as_slice())))
        {
            let hes: Vec<String> = lp
                .iter()
                .map(|&h| {
                    let n = graph.half_edges[h].origin;
                    let e = graph.half_edges[graph.half_edges[h].twin].origin;
                    let p = graph.nodes[n].position;
                    let shared = matches!(
                        graph.half_edges[h].source,
                        super::super::graph::EdgeSource::Shared { .. }
                    );
                    format!(
                        "{h}:n{n}->n{e}({:.2},{:.2},{:.2}){}",
                        p.x,
                        p.y,
                        p.z,
                        if shared { "S" } else { "" }
                    )
                })
                .collect();
            eprintln!("  [loopdump-all] {face_ref:?} {tag} [{}]", hes.join(" "));
        }
    }
    if popkern_core::trace::trace_enabled("loop-dump") && outer_loop.len() >= 20 {
        eprintln!("  [loopdump] {face_ref:?} outer {} hes:", outer_loop.len());
        for &h in outer_loop {
            let n = graph.half_edges[h].origin;
            let p = graph.nodes[n].position;
            let shared = matches!(
                graph.half_edges[h].source,
                super::super::graph::EdgeSource::Shared { .. }
            );
            eprintln!(
                "    he {h} ({:.1},{:.1},{:.1}) shared={shared}",
                p.x, p.y, p.z
            );
        }
    }
}
