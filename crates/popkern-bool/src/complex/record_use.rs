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
//!
//! The leak needs the unsplit face to pass THROUGH the edge. A record
//! piece outside that face's trim (in a hole, past its rim) splits only
//! the other face, and both sides of that split lie on the same side of
//! the unsplit face's operand: a phantom split, harmless to membership
//! (measured: a keyway slot's wall split along the block's bottom plane
//! across the bore, where the bottom face has its hole). Only a piece
//! outside the unused face's trim is exempt; one running along its rim
//! should have merged with the rim edge and is still reported. The
//! piece's ends may lie on the rim: a record leaving the face does
//! (measured: a tangent boss's records through a web's notch in a ring
//! wall, ending on the notch rim at the tangent point).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

use crate::merge::graph::{EdgeSource, IntersectionGraph};

use super::ids::EdgeId;
use super::model::{CellComplex, ComplexEdge};

/// Complex edges lifted from a two-sided shared record that only one of
/// the record's faces uses, unless the record lies wholly outside the
/// unused face's trim.
pub(crate) fn half_registered_records(
    c: &CellComplex,
    graph: &IntersectionGraph,
    operands: [&BRep; 2],
    linear: f64,
) -> Vec<EdgeId> {
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
        let (used_a, used_b) = (uses(0, *face_a), uses(1, *face_b));
        let unused = match (used_a, used_b) {
            (true, false) => (1, *face_b),
            (false, true) => (0, *face_a),
            _ => continue,
        };
        // An unanswerable trim query keeps the decline: no evidence
        // the record misses the face.
        if passes_through(c, edge, operands[unused.0], unused.1, linear) != Some(false) {
            out.push(EdgeId::from_index(ei));
        } else if popkern_core::trace::trace_enabled("complex") {
            eprintln!(
                "[complex] record edge {ei} splits one face only; it misses operand {} {:?}",
                unused.0, unused.1
            );
        }
    }
    out
}

/// Whether `edge` passes through `face`'s trim: a sample strictly
/// inside it, or every quarter point inside or on its rim (a piece
/// riding the rim). `None` when the face or a sample cannot be
/// classified.
fn passes_through(
    c: &CellComplex,
    edge: &ComplexEdge,
    brep: &BRep,
    face: FaceHandle,
    linear: f64,
) -> Option<bool> {
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let (pa, pb) = (
        c.vertices[edge.start.index()].position,
        c.vertices[edge.end.index()].position,
    );
    let mut along_rim = true;
    for s in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let p = match &edge.curve {
            Some(cv) => cv
                .geometry
                .point_at(cv.t_start + (cv.t_end - cv.t_start) * s),
            None => pa + (pb - pa) * s,
        };
        let (u, v) = surface.closest_point(&p)?;
        // A record of this face lies on its surface; one that does not
        // is no evidence either way.
        if surface.point_at(u, v).distance_to(&p) > linear {
            return None;
        }
        let at = brep.face_uv_trim_containment(face, surface, u, v)?;
        if at == FaceUvContainment::Interior {
            return Some(true);
        }
        if s != 0.0 && s != 1.0 && at == FaceUvContainment::Outside {
            along_rim = false;
        }
    }
    // Along the face's rim the record should have merged with the rim
    // edge: no evidence it misses the face. Its ends on the rim are where
    // any record leaving the face lies, and a piece outside the trim
    // that only touches the rim between them is tangent to it there;
    // neither says anything.
    Some(along_rim)
}

#[cfg(test)]
#[path = "record_use_tests.rs"]
mod tests;
