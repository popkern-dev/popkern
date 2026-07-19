//! PK-152 pair-corner repair: shared Euler surgery and the stage-1
//! equal-radius analytic path (see `blend_seam_corner` for the
//! recognizer and `blend_seam_corner_march` for the mismatched-radius
//! marching path).

use popkern_core::geometry::Point3;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::BlendSeamCornerDeclineReason as Reason;

use super::blend_seam_corner::PairCornerCandidate;
use super::blend_seam_corner_geometry as geometry;

/// Retire the flap and stub chains with Euler operators only, leaving
/// the seam edge running `anchor_near` → `anchor_far`. `kef` keeps
/// the given half-edge's face, so the flap (the twin side) dies into
/// `strip_n`; the remaining `chain_n` stubs become spurs dissolved
/// tip-first; `chain_p` contracts onto the far anchor, whose position
/// `kve` preserves (origin side).
pub(super) fn perform_flap_surgery(brep: &mut BRep, c: &PairCornerCandidate) -> Result<(), Reason> {
    let last_n = *c.chain_n.last().ok_or(Reason::EulerRejected)?;
    brep.kef(last_n).map_err(|_| Reason::EulerRejected)?;
    for he in c.chain_n[..c.chain_n.len() - 1].iter().rev() {
        brep.kev(*he).map_err(|_| Reason::EulerRejected)?;
    }
    for he in &c.chain_p {
        brep.kve(*he).map_err(|_| Reason::EulerRejected)?;
    }
    // The seam edge must now run between the two anchors.
    let he_n_data = brep.half_edges.get(c.he_n).ok_or(Reason::EulerRejected)?;
    let he_p_data = brep
        .half_edges
        .get(he_n_data.twin)
        .ok_or(Reason::EulerRejected)?;
    if he_n_data.origin != c.anchor_near || he_p_data.origin != c.anchor_far {
        return Err(Reason::EulerRejected);
    }
    Ok(())
}

/// Commit one edge's re-authored curve and per-side pcurves.
pub(super) fn commit_edge_authoring(
    brep: &mut BRep,
    edge: popkern_core::topology::EdgeHandle,
    curve: popkern_core::geometry::CurveOnEdge,
    staged: Vec<(FaceHandle, HalfEdgeHandle, Pcurve)>,
) {
    if let Some(e) = brep.edges.get_mut(edge) {
        e.curve = Some(curve);
    }
    for (face, he, pcurve) in staged {
        if let Some(h) = brep.half_edges.get_mut(he) {
            h.pcurve = Some(pcurve.clone());
        }
        brep.set_face_use_pcurve(face, he, pcurve);
    }
}

/// Stage-1 equal-radius repair: analytic bisector-plane ellipse seam.
pub(super) fn repair_pair_corner_equal(
    brep: &mut BRep,
    c: &PairCornerCandidate,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let cyl_n =
        super::blend_seam_corner::strip_cylinder(brep, c.strip_n).ok_or(Reason::EulerRejected)?;
    let cyl_p =
        super::blend_seam_corner::strip_cylinder(brep, c.strip_p).ok_or(Reason::EulerRejected)?;
    let pos_of = |v: VertexHandle| -> Option<Point3> { Some(brep.vertices.get(v)?.position) };
    let corner_pos = pos_of(c.corner).ok_or(Reason::EulerRejected)?;
    let a_pos = pos_of(c.anchor_near).ok_or(Reason::EulerRejected)?;
    let b_pos = pos_of(c.anchor_far).ok_or(Reason::EulerRejected)?;

    // Author the seam geometry FIRST — every analytic decline exits
    // before any topology is touched.
    let frame = geometry::solve_equal_radius_seam_frame(&cyl_n, &cyl_p, a_pos, b_pos, tolerance)?;
    let primary = brep
        .edges
        .get(c.seam_edge)
        .ok_or(Reason::EulerRejected)?
        .half_edge;
    // After surgery `he_n` runs anchor_near → anchor_far and `he_p`
    // the reverse; the stored curve follows the primary half-edge.
    let (p_start, p_end) = if primary == c.he_n {
        (a_pos, b_pos)
    } else {
        (b_pos, a_pos)
    };
    let seam_curve = geometry::author_seam_curve(&frame, p_start, p_end, corner_pos, tolerance)?;

    perform_flap_surgery(brep, c)?;

    // Strip-side pcurves, anchored to the loop neighbours' charts.
    let staged: Vec<(FaceHandle, HalfEdgeHandle, Pcurve)> = [
        (c.strip_n, c.he_n, a_pos, b_pos),
        (c.strip_p, c.he_p, b_pos, a_pos),
    ]
    .into_iter()
    .map(|(face, he, origin_pos, dest_pos)| {
        let surface = brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
            .ok_or(Reason::EulerRejected)?;
        let h = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?;
        let origin_anchor = geometry::ChartAnchor {
            uv: neighbour_chart_point(brep, face, h.prev, true),
        };
        let dest_anchor = geometry::ChartAnchor {
            uv: neighbour_chart_point(brep, face, h.next, false),
        };
        let pcurve = geometry::author_side_pcurve(
            &surface,
            &frame,
            origin_pos,
            dest_pos,
            &origin_anchor,
            &dest_anchor,
            tolerance,
        )?;
        Ok((face, he, pcurve))
    })
    .collect::<Result<_, Reason>>()?;

    commit_edge_authoring(brep, c.seam_edge, seam_curve, staged);

    // Postconditions: the seam must sit on both strips and the body
    // must still validate.
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if brep.measure_edge_deviation(c.seam_edge) > pos_tol {
        return Err(Reason::ValidationFailed);
    }
    brep.validate().map_err(|_| Reason::ValidationFailed)?;
    Ok(())
}

/// Authoritative chart position of a seam endpoint on one face: the
/// adjacent half-edge's pcurve endpoint (its end when the neighbour
/// arrives at the vertex, its start when it leaves).
pub(super) fn neighbour_chart_point(
    brep: &BRep,
    face: FaceHandle,
    neighbour: HalfEdgeHandle,
    at_end: bool,
) -> Option<(f64, f64)> {
    let pc = brep.face_use_pcurve(face, neighbour).or(brep
        .half_edges
        .get(neighbour)?
        .pcurve
        .as_ref())?;
    let (t0, t1) = pc.parameter_range();
    let p = pc.point_at(if at_end { t1 } else { t0 });
    Some((p.x, p.y))
}
