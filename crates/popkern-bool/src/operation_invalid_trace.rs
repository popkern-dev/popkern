//! Invalid-result trace dump for the boolean pipeline (loops, edges
//! with curve endpoints/midpoints, face->surface bindings) — split
//! from `operation_helpers.rs` under the 500-line file budget.

/// Boolean-trace dump of a validation-declined result: every loop's
/// half-edge chain with positions, and every edge's curve endpoints
/// against its vertex pair — the wiring detail post-check findings
/// cannot show.
pub(crate) fn trace_invalid_result(result: &popkern_core::topology::BRep) {
    if !popkern_core::trace::trace_enabled("boolean") {
        return;
    }
    for (lh, l) in result.loops.iter() {
        let mut chain = Vec::new();
        let mut he = l.half_edge;
        for _ in 0..result.half_edges.len() + 1 {
            let Some(h) = result.half_edges.get(he) else {
                break;
            };
            let p = result.vertices.get(h.origin).map(|v| v.position);
            chain.push((he, h.edge, p));
            he = h.next;
            if he == l.half_edge {
                break;
            }
        }
        eprintln!("  [invalid-loop] {lh:?} {} hes:", chain.len());
        for (heh, eh, p) in &chain {
            if let Some(p) = p {
                eprintln!(
                    "    he {heh:?} edge {eh:?} origin ({:.4},{:.4},{:.4})",
                    p.x, p.y, p.z
                );
            }
        }
    }
    use popkern_core::geometry::curve::CurveEval as _;
    for (eh, e) in result.edges.iter() {
        let Some(c) = e.curve.as_ref() else { continue };
        let p0 = c.geometry.point_at(c.t_start);
        let p1 = c.geometry.point_at(c.t_end);
        let pm = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
        let he = result.half_edges.get(e.half_edge);
        let v0 = he
            .and_then(|h| result.vertices.get(h.origin))
            .map(|v| v.position);
        let v1 = he
            .and_then(|h| result.half_edges.get(h.twin))
            .and_then(|t| result.vertices.get(t.origin))
            .map(|v| v.position);
        eprintln!(
            "  [invalid-edge] {eh:?} curve ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) mid ({:.4},{:.4},{:.4}) t [{:.6},{:.6}] v0 {:?} v1 {:?}",
            p0.x,
            p0.y,
            p0.z,
            p1.x,
            p1.y,
            p1.z,
            pm.x,
            pm.y,
            pm.z,
            c.t_start,
            c.t_end,
            v0.map(|p| format!("({:.6},{:.6},{:.6})", p.x, p.y, p.z)),
            v1.map(|p| format!("({:.6},{:.6},{:.6})", p.x, p.y, p.z)),
        );
    }
    // Loop → face → surface bindings: an invalid edge's loop dump is
    // unreadable without knowing WHICH surface each loop trims (a
    // phantom edge class is diagnosed by the face it claims to lie
    // on).
    for (fh, f) in result.faces.iter() {
        let kind = f
            .surface
            .as_ref()
            .map(|s| s.geometry.kind_name())
            .unwrap_or("<none>");
        let loops: Vec<String> = [f.outer_loop]
            .into_iter()
            .chain(f.inner_loops.iter().copied())
            .map(|lh| format!("{lh:?}"))
            .collect();
        eprintln!("  [invalid-face] {fh:?} surface {kind} loops {loops:?}");
    }
}
