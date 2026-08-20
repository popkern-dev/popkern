//! Boundary measurement helpers for the sliver-absorption
//! classifiers (split from `sliver_absorb_classifiers.rs` for the
//! 500-line structural cap).
//!
//! Measurement contract (PK-296): vertex positions are the
//! pipeline's live truth — mid-pipeline surgery moves vertices
//! without re-authoring curves, so a genuinely collapsed sliver can
//! carry stale macroscopic curves. Stored curves arbitrate only
//! where the vertex skeleton is STRUCTURALLY silent: topologically
//! closed edges (hole rims, closed-fillet contact circles), whose
//! endpoint span is 0 by construction rather than by collapse.

use popkern_core::topology::{BRep, FaceHandle};

/// True when the outer loop of `face` walks a topologically CLOSED
/// edge (half-edge origin == twin origin — a full circle such as a
/// hole rim or a closed-fillet contact circle). Faces whose
/// degenerate vertex measure involves a closed edge get a second
/// look through their curves before being classified as slivers.
pub(super) fn outer_loop_has_closed_edge(brep: &BRep, face: FaceHandle) -> bool {
    let Some(fd) = brep.faces.get(face) else {
        return false;
    };
    let Some(ld) = brep.loops.get(fd.outer_loop) else {
        return false;
    };
    let start = ld.half_edge;
    if !brep.half_edges.contains_key(start) {
        return false;
    }
    let mut he = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let hd = &brep.half_edges[he];
        if let Some(twin) = brep.half_edges.get(hd.twin)
            && twin.origin == hd.origin
        {
            return true;
        }
        he = hd.next;
        if he == start {
            break;
        }
    }
    false
}

/// Curve-sampled outer-loop boundary polygon. Half-edges with a
/// stored curve contribute 9 stations oriented along the half-edge's
/// walk direction; curve-less half-edges contribute their origin
/// vertex. Used only for the closed-edge second look (see
/// [`outer_loop_has_closed_edge`]) — same defect class as the
/// popkern-bool sliver_removal bigon fix (PK-295); the helper is
/// intentionally the same shape, duplicated because bool and fillet
/// cannot share non-core code.
fn sampled_outer_loop_polygon(
    brep: &BRep,
    face: FaceHandle,
) -> Vec<popkern_core::geometry::Point3> {
    use popkern_core::geometry::curve::CurveEval;
    let Some(fd) = brep.faces.get(face) else {
        return Vec::new();
    };
    let Some(ld) = brep.loops.get(fd.outer_loop) else {
        return Vec::new();
    };
    let start = ld.half_edge;
    if !brep.half_edges.contains_key(start) {
        return Vec::new();
    }
    let mut poly = Vec::new();
    let mut he = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let hd = &brep.half_edges[he];
        let origin = brep.vertices[hd.origin].position;
        if let Some(curve) = brep.edges.get(hd.edge).and_then(|e| e.curve.as_ref()) {
            let eval: &dyn CurveEval = &curve.geometry;
            // Orient samples along THIS half-edge: the stored curve
            // runs the edge's canonical direction, which is the
            // twin's direction half the time.
            let forward = eval.point_at(curve.t_start).distance_to(&origin)
                <= eval.point_at(curve.t_end).distance_to(&origin);
            for i in 0..9 {
                let f = i as f64 / 9.0;
                let f = if forward { f } else { 1.0 - f };
                poly.push(eval.point_at(curve.t_start + (curve.t_end - curve.t_start) * f));
            }
        } else {
            poly.push(origin);
        }
        he = hd.next;
        if he == start {
            break;
        }
    }
    poly
}

/// Magnitude-sum fan area of the curve-sampled boundary polygon.
/// Magnitude sum, not vector sum: signed vector areas of
/// opposite-sense rings cancel exactly on closed bands (a full
/// cylinder wall's top and bottom rims), which would re-introduce
/// the zero-area false positive the second look exists to remove.
/// The magnitude sum can only over-estimate, so it errs toward NOT
/// classifying a face as a sliver — the safe direction for a pass
/// whose action is to kef the face away.
pub(super) fn sampled_boundary_area(brep: &BRep, face: FaceHandle) -> f64 {
    let positions = sampled_outer_loop_polygon(brep, face);
    if positions.len() < 3 {
        return 0.0;
    }
    let p0 = positions[0];
    let mut area = 0.0_f64;
    for i in 1..positions.len() - 1 {
        let v1 = positions[i] - p0;
        let v2 = positions[i + 1] - p0;
        area += v1.cross(&v2).length() * 0.5;
    }
    area
}
