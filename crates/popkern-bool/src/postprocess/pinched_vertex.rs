//! Pinched vertices of a served boolean result (PK-992).
//!
//! The half-edges leaving a manifold vertex form ONE cycle under
//! `h -> twin(h).next`: the faces around it close up into a single
//! fan, the vertex's link is a disk. A vertex where they form two or
//! more cycles is shared by separate fans: two cones of material
//! touching at a point that the result carries as ONE vertex, which no
//! coincident-copy test sees. Measured: main − branch of an
//! equal-radius tee whose branch passes through the main cylinder
//! leaves the main body's material in two lobes, |z − 60| > |x|,
//! meeting only at the singular points (0, ±40, 60), and the cellular
//! route served that body as two-manifold edges around a pinched
//! vertex. Under decision D1 it declines with the vertex as a point
//! contact.

use popkern_core::topology::{BRep, ContactEntity, HalfEdgeHandle, VertexHandle};

use crate::error::{ContactKind, NonManifoldContact};

/// Vertices of `brep` whose leaving half-edges form more than one
/// fan, in vertex order, each as a point contact.
pub(crate) fn pinched_vertices(brep: &BRep) -> Vec<(NonManifoldContact, Vec<ContactEntity>)> {
    let mut leaving: Vec<(VertexHandle, HalfEdgeHandle)> = brep
        .half_edges
        .iter()
        .map(|(h, he)| (he.origin, h))
        .collect();
    leaving.sort();
    let mut out = Vec::new();
    for run in leaving.chunk_by(|a, b| a.0 == b.0) {
        let v = run[0].0;
        let hs: Vec<HalfEdgeHandle> = run.iter().map(|&(_, h)| h).collect();
        if hs.len() < 2 || fan_count(brep, v, &hs) < 2 {
            continue;
        }
        let Some(vertex) = brep.vertices.get(v) else {
            continue;
        };
        out.push((
            NonManifoldContact {
                kind: ContactKind::Vertex,
                location: vertex.position,
            },
            vec![ContactEntity::Vertex(v)],
        ));
    }
    out
}

/// Cycles of `h -> twin(h).next` among the half-edges `hs` leaving `v`.
/// A walk that leaves `v` or meets a missing twin (an open sheet's
/// rim) ends its fan there.
fn fan_count(brep: &BRep, v: VertexHandle, hs: &[HalfEdgeHandle]) -> usize {
    // A vertex has a handful of half-edges: a list beats a tree.
    let mut seen: Vec<HalfEdgeHandle> = Vec::with_capacity(hs.len());
    let mut fans = 0;
    for &start in hs {
        if seen.contains(&start) {
            continue;
        }
        fans += 1;
        let mut cur = start;
        for _ in 0..=hs.len() {
            seen.push(cur);
            let Some(next) = brep
                .half_edges
                .get(cur)
                .and_then(|h| h.mate().and_then(|m| brep.half_edges.get(m)))
                .map(|tw| tw.next)
            else {
                break;
            };
            if next == start
                || seen.contains(&next)
                || brep.half_edges.get(next).is_none_or(|h| h.origin != v)
            {
                break;
            }
            cur = next;
        }
    }
    fans
}
