//! Env-gated election tracing (PK_LOOP_DUMP / PK_LOOP_DUMP_ALL),
//! extracted under the file-size ratchet.

/// PK_LOOP_DUMP / PK_LOOP_DUMP_ALL election tracing (extracted from
/// the walk body under the file-size ratchet).
pub(super) fn dump_election_loops(
    graph: &super::super::graph::IntersectionGraph,
    face_ref: &super::super::graph::FaceRef,
    outer_loop: &[usize],
    inner_loops_out: &[Vec<usize>],
    separate_outers: &[Vec<usize>],
) {
    if std::env::var_os("PK_LOOP_DUMP_ALL").is_some() {
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
    if std::env::var_os("PK_LOOP_DUMP").is_some() && outer_loop.len() >= 20 {
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
