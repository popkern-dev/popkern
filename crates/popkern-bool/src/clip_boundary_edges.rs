//! Face boundary-edge inventories for the clip's junction and
//! crossing certificates (split from `clip_endpoint_refine.rs` under
//! the 500-line file budget).

use popkern_core::topology::{BRep, FaceHandle};

/// The OUTER loop's edge handles only (the refinement's lateral
/// junction search walks this; hole rims go through
/// [`boundary_edges`]).
pub(crate) fn boundary_edges_outer(
    brep: &BRep,
    face: FaceHandle,
) -> Vec<popkern_core::topology::EdgeHandle> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let Some(lp) = brep.loops.get(face_data.outer_loop) else {
        return out;
    };
    let start = lp.half_edge;
    let mut he = start;
    for _ in 0..brep.half_edges.len().max(1) {
        let Some(h) = brep.half_edges.get(he) else {
            break;
        };
        out.push(h.edge);
        he = h.next;
        if he == start {
            break;
        }
    }
    out
}

/// [`boundary_edges`] restricted to the face's TRIM boundary: SEAM
/// edges — both half-edge uses on this same face (a periodic wall's
/// closure line) — are excluded. A record crossing the seam stays ON
/// the face, so the seam is never a trim-crossing event; feeding it
/// to the crossing/junction certificates truncated a groove-runout
/// record at its seam hit (the far retrim dragged the
/// junction-exact end n86 back 0.59 onto the seam junction while
/// every removed-span sample abstained by "riding" the seam), and
/// the missing tail stranded the LH two-start's cellular assembly
/// (PK-491, 24 once-used sides).
pub(crate) fn trim_boundary_edges(
    brep: &BRep,
    face: FaceHandle,
) -> Vec<popkern_core::topology::EdgeHandle> {
    let is_seam = |eh: popkern_core::topology::EdgeHandle| -> bool {
        let Some(edge) = brep.edges.get(eh) else {
            return false;
        };
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            return false;
        };
        let Some(tw) = brep.half_edges.get(he.twin) else {
            return false;
        };
        let face_of = |lh: popkern_core::topology::LoopHandle| brep.loops.get(lh).map(|l| l.face);
        match (face_of(he.loop_), face_of(tw.loop_)) {
            (Some(fa), Some(fb)) => fa == fb && fa == face,
            _ => false,
        }
    };
    let mut out = boundary_edges(brep, face);
    out.retain(|&eh| !is_seam(eh));
    out
}

/// Every boundary edge handle of `face` — outer AND inner loops. The
/// refinement's junction search must see hole rims too: a record
/// exiting through a bore mouth crosses the top face's INNER loop
/// (measured, PK-465 — the angled tap's exit records cross the mouth
/// rim, which the outer-loop-only walk never visited).
pub(crate) fn boundary_edges(
    brep: &BRep,
    face: FaceHandle,
) -> Vec<popkern_core::topology::EdgeHandle> {
    let mut out = Vec::new();
    let Some(face_data) = brep.faces.get(face) else {
        return out;
    };
    let loops = std::iter::once(face_data.outer_loop).chain(face_data.inner_loops.iter().copied());
    for lh in loops {
        let Some(lp) = brep.loops.get(lh) else {
            continue;
        };
        let start = lp.half_edge;
        let mut he = start;
        for _ in 0..brep.half_edges.len().max(1) {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            out.push(h.edge);
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    out
}
