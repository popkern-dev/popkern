//! A deterministic JSON dump of a cell complex for diagnostics.
//!
//! Emitted with a decline (and under `POPKERN_TRACE=complex`), so a
//! failed boolean can be inspected and reproduced without rerunning it:
//! vertices, edges with their radial fin order, fins, faces with owners
//! and side shells, and shells with membership. Hand-written (no serde
//! dependency in this crate); key order and number formatting are
//! fixed, so equal complexes dump byte-identically.

use std::fmt::Write;

use super::model::{CellComplex, SideIndex};

/// The complex as a JSON object.
pub(crate) fn complex_artifact_json(c: &CellComplex) -> String {
    let mut s = String::new();
    let _ = write!(s, "{{\"operands\":{},", c.operand_count);
    s.push_str("\"vertices\":[");
    for (i, v) in c.vertices.iter().enumerate() {
        sep(&mut s, i);
        let _ = write!(
            s,
            "{{\"id\":{i},\"p\":[{:e},{:e},{:e}],\"tol\":{:e},\"node\":{}}}",
            v.position.x,
            v.position.y,
            v.position.z,
            v.tolerance,
            opt(v.source_node)
        );
    }
    s.push_str("],\"edges\":[");
    for (i, e) in c.edges.iter().enumerate() {
        sep(&mut s, i);
        let radial: Vec<String> = e.radial.iter().map(|f| f.0.to_string()).collect();
        let _ = write!(
            s,
            "{{\"id\":{i},\"v\":[{},{}],\"curve\":{},\"tol\":{:e},\"radial\":[{}],\"ordered\":{},\"key\":{},\"frontier\":{}}}",
            e.start.0,
            e.end.0,
            e.curve.is_some(),
            e.tolerance,
            radial.join(","),
            e.radial_ordered,
            opt(e.source_key),
            e.frontier
        );
    }
    s.push_str("],\"fins\":[");
    for (i, f) in c.fins.iter().enumerate() {
        sep(&mut s, i);
        let _ = write!(
            s,
            "{{\"id\":{i},\"edge\":{},\"face\":{},\"loop\":{},\"fwd\":{},\"next\":{},\"he\":{}}}",
            f.edge.0,
            f.face.0,
            f.loop_.0,
            f.forward,
            f.next.0,
            opt(f.source_half_edge)
        );
    }
    s.push_str("],\"faces\":[");
    for (i, f) in c.faces.iter().enumerate() {
        sep(&mut s, i);
        let owners: Vec<String> = f
            .owners
            .iter()
            .map(|o| {
                format!(
                    "{{\"operand\":{},\"face\":\"{:?}\",\"same\":{}}}",
                    o.operand, o.face, o.same_sense
                )
            })
            .collect();
        let loops: Vec<String> = f.loops.iter().map(|l| l.0.to_string()).collect();
        let frags: Vec<String> = f.source_fragments.iter().map(|x| x.to_string()).collect();
        let _ = write!(
            s,
            "{{\"id\":{i},\"owners\":[{}],\"loops\":[{}],\"ahead\":{},\"behind\":{},\"folded_out\":{},\"fragments\":[{}]}}",
            owners.join(","),
            loops.join(","),
            opt(f.sides[SideIndex::Ahead as usize].shell.map(|x| x.index())),
            opt(f.sides[SideIndex::Behind as usize].shell.map(|x| x.index())),
            f.folded_is_outward,
            frags.join(",")
        );
    }
    s.push_str("],\"shells\":[");
    for (i, sh) in c.shells.iter().enumerate() {
        sep(&mut s, i);
        let sides: Vec<String> = sh
            .sides
            .iter()
            .map(|(f, side)| format!("[{},{}]", f.0, *side as u8))
            .collect();
        let membership = sh.membership.as_ref().map_or("null".to_string(), |m| {
            let ks: Vec<String> = m.iter().map(|k| k.to_string()).collect();
            format!("[{}]", ks.join(","))
        });
        let _ = write!(
            s,
            "{{\"id\":{i},\"sides\":[{}],\"membership\":{membership}}}",
            sides.join(",")
        );
    }
    s.push_str("]}");
    s
}

fn sep(s: &mut String, i: usize) {
    if i > 0 {
        s.push(',');
    }
}

fn opt(v: Option<usize>) -> String {
    v.map_or("null".to_string(), |x| x.to_string())
}
