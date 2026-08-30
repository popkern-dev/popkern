//! Outer-loop half-edge walk shared by the setback vertex topology
//! and the boundary-ring construction.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// The half-edges of `face`'s outer loop in cycle order.
pub(crate) fn face_loop_half_edges(
    brep: &BRep,
    face: FaceHandle,
) -> Result<Vec<HalfEdgeHandle>, String> {
    let f = brep.faces.get(face).ok_or("face handle invalid")?;
    let start = brep
        .loops
        .get(f.outer_loop)
        .ok_or("outer loop invalid")?
        .half_edge;
    let mut out = Vec::new();
    let mut he = start;
    for _ in 0..10_000 {
        out.push(he);
        he = brep.half_edges[he].next;
        if he == start {
            return Ok(out);
        }
    }
    Err("unterminated face loop".to_string())
}
