//! Every shared intersection record must split both faces it lies on
//! (plan Stage 3, "phantom" findings).
//!
//! A record is the curve where a face of A meets a face of B; the walk
//! splits both faces along it. When one face's walk splits along it and
//! the other's does not, the complex gets an edge whose fins come from
//! one side only, so the regions on either side of the unsplit face
//! leak into each other through it — the shells come out with a face's
//! two sides joined, and the decline names the symptom, not the cause
//! (measured: a cuboid bite through a rim-filleted cylinder, the cutter
//! wall split along its torus record while the torus was not). This
//! check names the cause. Trim-residue records partition one face by
//! design (both halves registered to their owner) and are exempt.

use crate::merge::graph::{EdgeSource, IntersectionGraph};

use super::ids::EdgeId;
use super::model::CellComplex;

/// Complex edges lifted from a two-sided shared record that only one of
/// the record's faces uses.
pub(crate) fn half_registered_records(c: &CellComplex, graph: &IntersectionGraph) -> Vec<EdgeId> {
    let mut out = Vec::new();
    for (ei, edge) in c.edges.iter().enumerate() {
        let Some(key) = edge.source_key else {
            continue;
        };
        let he = &graph.half_edges[key];
        let EdgeSource::Shared { face_a, face_b } = &he.source else {
            continue;
        };
        if he.face == graph.half_edges[he.twin].face {
            continue;
        }
        let uses = |operand: usize, face| {
            edge.radial.iter().any(|fin| {
                c.faces[c.fins[fin.index()].face.index()]
                    .owners
                    .iter()
                    .any(|o| o.operand == operand && o.face == face)
            })
        };
        if uses(0, *face_a) != uses(1, *face_b) {
            out.push(EdgeId::from_index(ei));
        }
    }
    out
}
