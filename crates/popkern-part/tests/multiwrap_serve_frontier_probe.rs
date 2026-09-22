//! Multi-wrap serve frontier probe (PK-461 -> PK-494): measure where
//! thread-cut serves break as tool turn count grows, with per-face
//! open-edge attribution and served-topology dumps. Run it at the
//! lengths under repair and watch for torn serves (open_edges>0 is
//! WORSE than a decline).
//!
//! Run by hand:
//! `PK461_PROBE_LEN=<len> cargo test -p popkern-part --release --test multiwrap_serve_frontier_probe -- --ignored --nocapture`
//!
//! Frontier at PK-494 close: M6 rods serve watertight through len 80
//! (81 tool turns) with the record-side classification witness; at
//! len 100 the merge's strand gate declines honestly (a structured
//! half-claimed-records decline — never a torn serve).

use popkern_core::tolerance::Tolerance;
use popkern_part::fastener::{ThreadSpec, ThreadedRod};

#[test]
#[ignore = "PK-494 multi-wrap frontier probe — run by hand with PK461_PROBE_LEN"]
fn probe_turn_scale() {
    let len: f64 = std::env::var("PK461_PROBE_LEN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(20.0);
    let spec = ThreadSpec::parse("M6").expect("M6");
    match ThreadedRod::new(spec, len)
        .tip_chamfer(false)
        .build(Tolerance::default())
    {
        Ok(rod) => {
            let open = rod.open_edge_count().expect("mesh");
            println!("PK461_PROBE len={len} SERVED open_edges={open}");
            {
                let b = rod.brep();
                println!(
                    "PK461_TOPO faces={} edges={} loops={}",
                    b.faces.len(),
                    b.edges.len(),
                    b.loops.len()
                );
                for (fh, f) in b.faces.iter() {
                    let mut n = 0;
                    let mut he = b.loops[f.outer_loop].half_edge;
                    for _ in 0..b.half_edges.len() + 1 {
                        n += 1;
                        he = b.half_edges[he].next;
                        if he == b.loops[f.outer_loop].half_edge {
                            break;
                        }
                    }
                    let kind = f
                        .surface
                        .as_ref()
                        .map(|s| s.geometry.kind_name())
                        .unwrap_or("<none>");
                    if n > 6 || !f.inner_loops.is_empty() {
                        println!(
                            "PK461_FACE {fh:?} {kind} outer_hes={n} inners={}",
                            f.inner_loops.len()
                        );
                    }
                }
                if open > 0 {
                    report_open_edges_per_face(b);
                }
            }
        }
        Err(e) => println!("PK461_PROBE len={len} DECLINED: {e:?}"),
    }
}

/// Attribute welded open mesh edges to their originating faces (same
/// welding as popkern-part's watertightness gate).
fn report_open_edges_per_face(brep: &popkern_core::topology::BRep) {
    use std::collections::HashMap;
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut id = |p: popkern_core::geometry::Point3| -> usize {
        let k = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        let next = key.len();
        *key.entry(k).or_insert(next)
    };
    let mut counts: HashMap<(usize, usize), (usize, Vec<usize>)> = HashMap::new();
    for (ti, t) in mesh.triangles.iter().enumerate() {
        let w = [
            id(mesh.vertices[t[0]].position),
            id(mesh.vertices[t[1]].position),
            id(mesh.vertices[t[2]].position),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            let e = counts
                .entry((a.min(b), a.max(b)))
                .or_insert((0, Vec::new()));
            e.0 += 1;
            e.1.push(ti);
        }
    }
    let mut per_face: HashMap<String, usize> = HashMap::new();
    let mut printed = 0;
    for ((a, _), (c, tris)) in counts.iter() {
        if *c == 1 {
            let fh = mesh.face_map[tris[0]];
            let kind = brep.faces[fh]
                .surface
                .as_ref()
                .map(|s| s.geometry.kind_name())
                .unwrap_or("<none>");
            *per_face.entry(format!("{fh:?} {kind}")).or_insert(0) += 1;
            if printed < 30 && a % 7 == 0 {
                let t = &mesh.triangles[tris[0]];
                let p = mesh.vertices[t[0]].position;
                println!("PK461_TEAR {fh:?} near ({:.3},{:.3},{:.3})", p.x, p.y, p.z);
                printed += 1;
            }
        }
    }
    let mut rows: Vec<(String, usize)> = per_face.into_iter().collect();
    rows.sort_by_key(|r| std::cmp::Reverse(r.1));
    for (face, n) in rows.iter().take(12) {
        println!("PK461_OPEN {face}: {n}");
    }
}
