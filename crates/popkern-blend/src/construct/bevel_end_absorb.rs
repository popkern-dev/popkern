//! Planar bevel reflex ends: the terminal triangle joins the end face
//! (PK-994).
//!
//! **The class.** A chamfer (planar bevel) whose end vertex is a reflex
//! corner of one support (a pocket or hole rim, a tunnel mouth, a
//! re-entrant L corner) with a sharp transversal edge. The anchor-split
//! stage lands that support's rail end on the transversal edge's line
//! BEYOND the corner ([`super::perp_anchor_rail_proximity::bevel_reflex_end`]),
//! and the terminal-disc stage carves the triangle (corner, rail end,
//! stub anchor) off the strip with one `mef`.
//!
//! **Why it merges.** All three triangle vertices lie on the end face's
//! plane: the corner and the stub anchor are on the transversal edges,
//! and the rail end is where the bevel crosses the transversal edge's
//! line, inside that plane. The exact result is the end face extended
//! by the triangle, not a second coplanar face: the removed wedge ends
//! flush with the end face, whose boundary now runs past the old
//! corner to the bevel. One `kef` across the stub edge merges the
//! triangle into the end face (the end face's half-edge survives, so it
//! keeps its handle, surface and name). The corner is then left between
//! two collinear edges with the same face pair, the transversal edge
//! and the spur; one `kve` contracts the spur into the rail end and the
//! transversal edge, re-authored as an exact segment, reaches it.
//!
//! **Invariants:** Euler operators only (χ-neutral: kef F−1 E−1, kve
//! V−1 E−1); the strip, the support and the end face keep their
//! handles. **Fallback:** a triangle that is not coplanar with a face
//! across one of its corner edges, or whose corner keeps a third edge,
//! is left as it was: still an exact, valid face, reported by the
//! `bevel-end` trace. A failed operator restores the snapshot and
//! returns the reason; the caller declines the operation.
//! **Tolerance:** the caller's blend tolerance, coplanarity and
//! collinearity gated at `linear × 10` (the terminal-disc gate).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveOnEdge, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

/// Merge the terminal triangle `disc` (holding `corner`) into the end
/// face it is coplanar with, then contract the corner out of the end
/// face's boundary. `Ok(true)` when the triangle was absorbed,
/// `Ok(false)` when the configuration does not qualify (the triangle
/// stays), `Err` when an operator failed (the body is restored).
pub(super) fn absorb_bevel_end_triangle(
    brep: &mut BRep,
    disc: FaceHandle,
    corner: VertexHandle,
    tolerance: &Tolerance,
) -> Result<bool, String> {
    let trace = popkern_core::trace::trace_enabled("bevel-end");
    let gate = (tolerance.linear() * 10.0).max(1e-9);
    let Some((h_out, h_in)) = corner_half_edges(brep, disc, corner) else {
        return Ok(false);
    };
    let disc_points = loop_points(brep, disc);
    // The stub is the corner edge whose other face holds the whole
    // triangle in its plane; the other corner edge is the spur.
    let mut stub: Option<(HalfEdgeHandle, FaceHandle, EdgeHandle)> = None;
    for (he, other) in [(h_out, h_in), (h_in, h_out)] {
        let twin = brep.half_edges[he].twin;
        let face = brep.loops[brep.half_edges[twin].loop_].face;
        if face_holds(brep, face, &disc_points, gate) {
            stub = Some((twin, face, brep.half_edges[other].edge));
            break;
        }
    }
    let Some((stub_on_end, end_face, spur_edge)) = stub else {
        if trace {
            eprintln!(
                "[bevel-end] {disc:?}: no coplanar end face across the corner; triangle kept"
            );
        }
        return Ok(false);
    };
    let snapshot = brep.clone();
    let restore = |brep: &mut BRep, reason: String| -> Result<bool, String> {
        *brep = snapshot.clone();
        Err(reason)
    };
    if let Err(e) = brep.kef(stub_on_end) {
        return restore(brep, format!("bevel end absorb: kef failed: {e}"));
    }
    match join_through_corner(brep, corner, spur_edge, end_face, gate) {
        Ok(()) => {}
        Err(reason) => return restore(brep, reason),
    }
    let _ = brep.recompute_face_pcurves(end_face);
    if trace {
        eprintln!("[bevel-end] triangle {disc:?} merged into end face {end_face:?}");
    }
    Ok(true)
}

/// The disc loop's half-edges leaving and entering `corner`.
fn corner_half_edges(
    brep: &BRep,
    face: FaceHandle,
    corner: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    let start = brep.loops.get(brep.faces.get(face)?.outer_loop)?.half_edge;
    let mut cur = start;
    for _ in 0..=brep.half_edges.len() {
        let he = brep.half_edges.get(cur)?;
        if he.origin == corner {
            return Some((cur, he.prev));
        }
        cur = he.next;
        if cur == start {
            break;
        }
    }
    None
}

fn loop_points(brep: &BRep, face: FaceHandle) -> Vec<Point3> {
    brep.face_outer_loop_positions(face)
}

fn face_holds(brep: &BRep, face: FaceHandle, points: &[Point3], gate: f64) -> bool {
    let Some(SurfaceGeometry::Plane { origin, normal, .. }) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    else {
        return false;
    };
    let Some(n) = normal.normalized() else {
        return false;
    };
    points.iter().all(|p| (*p - *origin).dot(&n).abs() <= gate)
}

/// After the merge the corner sits between the spur and the transversal
/// edge, both now bounding the support and `end_face`. Contract the
/// spur into its tip and re-author the transversal edge as the exact
/// segment to the tip. A corner that keeps a third edge, or whose two
/// edges are not one straight line, stays (still valid).
fn join_through_corner(
    brep: &mut BRep,
    corner: VertexHandle,
    spur_edge: EdgeHandle,
    end_face: FaceHandle,
    gate: f64,
) -> Result<(), String> {
    let outgoing: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == corner)
        .map(|(h, _)| h)
        .collect();
    if outgoing.len() != 2 {
        return Ok(());
    }
    let Some(&spur) = outgoing
        .iter()
        .find(|&&h| brep.half_edges[h].edge == spur_edge)
    else {
        return Ok(());
    };
    let Some(&keep) = outgoing.iter().find(|&&h| h != spur) else {
        return Ok(());
    };
    let dest = |brep: &BRep, h: HalfEdgeHandle| brep.half_edges[brep.half_edges[h].twin].origin;
    let faces_of = |brep: &BRep, h: HalfEdgeHandle| {
        let f1 = brep.loops[brep.half_edges[h].loop_].face;
        let f2 = brep.loops[brep.half_edges[brep.half_edges[h].twin].loop_].face;
        (f1.min(f2), f1.max(f2))
    };
    let pair = faces_of(brep, spur);
    if pair != faces_of(brep, keep) || (pair.0 != end_face && pair.1 != end_face) {
        return Ok(());
    }
    let (tip, far) = (dest(brep, spur), dest(brep, keep));
    let p_c = brep.vertices[corner].position;
    let (d_tip, d_far) = (
        brep.vertices[tip].position - p_c,
        brep.vertices[far].position - p_c,
    );
    // Collinear through the corner, on opposite sides of it.
    let (l_tip, l_far) = (d_tip.length(), d_far.length());
    if l_tip <= gate
        || l_far <= gate
        || d_tip.cross(&d_far).length() > gate * l_tip.max(l_far)
        || d_tip.dot(&d_far) >= 0.0
    {
        return Ok(());
    }
    let kept_edge = brep.half_edges[keep].edge;
    // Straight: its curve's interior stations lie on the chord line.
    let straight = brep.edges[kept_edge].curve.as_ref().is_none_or(|c| {
        [0.25, 0.5, 0.75].iter().all(|f| {
            let m = c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * f);
            (m - p_c).cross(&d_far).length() <= gate * l_far
        })
    });
    if !straight {
        return Ok(());
    }
    // kve on the half-edge tip → corner merges the corner into the tip.
    let tip_to_corner = brep.half_edges[spur].twin;
    brep.kve(tip_to_corner)
        .map_err(|e| format!("bevel end absorb: kve failed: {e}"))?;
    let (p_far, p_tip) = (brep.vertices[far].position, brep.vertices[tip].position);
    let he = brep
        .edges
        .get(kept_edge)
        .ok_or_else(|| "bevel end absorb: transversal edge vanished in kve".to_string())?
        .half_edge;
    let segment = if brep.half_edges[he].origin == far {
        CurveOnEdge::line_segment(p_far, p_tip, gate / 10.0)
    } else {
        CurveOnEdge::line_segment(p_tip, p_far, gate / 10.0)
    };
    if let Some(edge) = brep.edges.get_mut(kept_edge) {
        edge.curve = segment;
    }
    Ok(())
}
