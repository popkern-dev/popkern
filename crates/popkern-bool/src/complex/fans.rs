//! Vertex fans of a selected boundary (PK-975).
//!
//! At a vertex, every kept face passing through it contributes a
//! corner (the fin arriving at the vertex and its successor leaving
//! it). Two corners are in one fan when they share an edge the kept
//! faces use exactly twice — a manifold edge of the result. A vertex
//! whose corners form more than one fan is a non-manifold vertex of
//! the result (Lee & Lee 2001: one partial vertex per fan): two
//! bodies touching at a point, or along an edge the graph carries as
//! two coincident edges, show up here.

use std::collections::BTreeMap;

use super::ids::{EdgeId, FinId, VertexId};
use super::model::CellComplex;

/// Vertices where the kept faces (`kept_face[face]`) form more than
/// one fan, in vertex order.
pub(crate) fn split_vertices(c: &CellComplex, kept_face: &[bool]) -> Vec<VertexId> {
    let kept_fin = |f: FinId| kept_face[c.fins[f.index()].face.index()];
    let kept_uses: Vec<usize> = c
        .edges
        .iter()
        .map(|e| e.radial.iter().filter(|&&f| kept_fin(f)).count())
        .collect();
    // Corners per vertex: (arriving fin, leaving fin).
    let mut corners: BTreeMap<VertexId, Vec<(FinId, FinId)>> = BTreeMap::new();
    for (fi, fin) in c.fins.iter().enumerate() {
        let f = FinId::from_index(fi);
        if !kept_fin(f) {
            continue;
        }
        if let Some(v) = c.fin_end(f) {
            corners.entry(v).or_default().push((f, fin.next));
        }
    }
    let mut out = Vec::new();
    for (v, list) in corners {
        if list.len() < 2 {
            continue;
        }
        // Union corners through shared manifold edges.
        let mut parent: Vec<usize> = (0..list.len()).collect();
        fn find(p: &mut [usize], mut x: usize) -> usize {
            while p[x] != x {
                p[x] = p[p[x]];
                x = p[x];
            }
            x
        }
        let mut by_edge: BTreeMap<EdgeId, Vec<usize>> = BTreeMap::new();
        for (k, &(arrive, leave)) in list.iter().enumerate() {
            for f in [arrive, leave] {
                let e = c.fins[f.index()].edge;
                if kept_uses[e.index()] == 2 {
                    by_edge.entry(e).or_default().push(k);
                }
            }
        }
        for ks in by_edge.values() {
            for w in ks.windows(2) {
                let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
                if a != b {
                    parent[a.max(b)] = a.min(b);
                }
            }
        }
        let fans = (0..list.len())
            .filter(|&k| find(&mut parent, k) == k)
            .count();
        if fans > 1 {
            if popkern_core::trace::trace_enabled("contacts") {
                eprintln!(
                    "  [contacts] complex vertex {v:?} at {:?}: {fans} fans of {} corners",
                    c.vertices[v.index()].position,
                    list.len()
                );
                for (k, &(arrive, leave)) in list.iter().enumerate() {
                    let edge = |f: FinId| c.fins[f.index()].edge;
                    eprintln!(
                        "    corner {k} (fan {}) face {:?}: arrive {:?} ({} kept uses), leave {:?} ({} kept uses)",
                        find(&mut parent, k),
                        c.fins[arrive.index()].face,
                        edge(arrive),
                        kept_uses[edge(arrive).index()],
                        edge(leave),
                        kept_uses[edge(leave).index()]
                    );
                }
            }
            out.push(v);
        }
    }
    out
}
