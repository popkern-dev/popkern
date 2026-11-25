//! Face boundary-edge inventories for the clip's junction and
//! crossing certificates (split from `clip_endpoint_refine.rs` under
//! the 500-line file budget).

use popkern_core::geometry::curve::{CurveEval as _, CurveGeometry};
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
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
        let Some(tw) = he.mate().and_then(|m| brep.half_edges.get(m)) else {
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

/// Parameter window of `curve` around `tc` holding every point within
/// `reach` (chord) of `curve(tc)`. Every junction the acceptance
/// admits lies within 16 bands of the endpoint, which is within 8
/// bands of `curve(tc)`, so the windowed solve finds the same
/// candidates as a whole-edge solve; on a thread flank's full-length
/// boundary helix the whole-edge solve was 16% of a lone screw's
/// build (PK-799). The window is sized from the local speed with 4x
/// headroom (arc length at this reach equals chord to well under
/// that). It falls back to the whole range where the speed is
/// degenerate, and where it would clamp at an end of a closed edge:
/// the junction may then lie across the seam, at the range's other
/// end (a rim circle's start vertex).
pub(crate) fn junction_window(
    curve: &CurveGeometry,
    tc: f64,
    full: (f64, f64),
    reach: f64,
) -> (f64, f64) {
    let h = (full.1 - full.0) * 1e-6;
    let speed = curve
        .point_at((tc + h).min(full.1))
        .distance_to(&curve.point_at((tc - h).max(full.0)))
        / ((tc + h).min(full.1) - (tc - h).max(full.0)).max(f64::MIN_POSITIVE);
    if !(speed.is_finite() && speed > 0.0) {
        return full;
    }
    let half = 4.0 * reach / speed;
    let clamps = tc - half < full.0 || tc + half > full.1;
    if clamps && curve.point_at(full.0).distance_to(&curve.point_at(full.1)) <= reach {
        return full;
    }
    ((tc - half).max(full.0), (tc + half).min(full.1))
}

/// Whether `piece` is its closed `source` whole, its ends the source's
/// parameter closure rather than window ends: no boundary crossing
/// made them, so there is nothing to refine. The far retrim once moved
/// such an end onto a junction the loop passes mid-curve and kept a
/// sliver of the loop (PK-1143: a dimple's rim circle touching the
/// plate's edge, its closure 0.24 short of the touch, came back as the
/// 0.24 between them and the cut declined).
pub(crate) fn spans_closed_source(piece: &CurveOnSurface, source: &CurveOnSurface) -> bool {
    let (p0, p1) = piece.curve_3d.parameter_range();
    let (s0, s1) = source.curve_3d.parameter_range();
    let slack = 1e-9 * (s1 - s0).abs();
    source.is_closed && (p0 - s0).abs() <= slack && (p1 - s1).abs() <= slack
}

/// Whether `p`, the piece's `end`, is the open `source`'s own end there.
/// The intersector stopped the curve at it; it is no clip crossing, and
/// the far re-trim, which corrects a crossing's parameter, has nothing
/// to correct (measured, CC-V 68: an arm of a tangent crossing, its own
/// end at the crossing 4.4e-4 inside B's rim, was re-trimmed onto the
/// rim's junction 6.8e-4 along it, and the arm became a stub along the
/// rim).
pub(crate) fn is_open_source_end(
    source: &CurveOnSurface,
    end: bool,
    p: &popkern_core::geometry::Point3,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    let (s0, s1) = source.curve_3d.parameter_range();
    !source.is_closed
        && source
            .curve_3d
            .point_at(if end { s1 } else { s0 })
            .distance_to(p)
            <= tolerance.linear()
}

/// Whether `p`, the piece's `end`, is the open `source`'s own end at an
/// exact touch of the two surfaces (their normals parallel within the
/// angular tolerance): a tangent crossing the intersector refined and
/// split its branches at. No refinement moves it: its branches pass
/// within the tolerance of everything near it, a face boundary inside
/// the touch's patch among them (measured, CC-V 64: an arm's own end at
/// the crossing was dragged 3.1e-4 onto B's rim junction, and the arm's
/// piece there closed on itself at the junction).
pub(crate) fn is_exact_touch_end(
    source: &CurveOnSurface,
    end: bool,
    p: &popkern_core::geometry::Point3,
    (surf_a, surf_b): (
        &popkern_core::geometry::surface::SurfaceGeometry,
        &popkern_core::geometry::surface::SurfaceGeometry,
    ),
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval as _;
    let normal = |s: &popkern_core::geometry::surface::SurfaceGeometry| {
        let (u, v) = s.closest_point(p)?;
        s.normal_at(u, v).normalized()
    };
    is_open_source_end(source, end, p, tolerance)
        && normal(surf_a)
            .zip(normal(surf_b))
            .is_some_and(|(na, nb)| na.cross(&nb).length() <= tolerance.angular())
}

/// Whether `p` lies on a seam of `face` (an edge both of whose sides
/// bound it) within ten tolerances, away from the seam's end vertices.
/// An exact touch's end moves onto such a junction only: a seam is no
/// boundary, and a crossing within the tolerance of its face's seam is
/// on it, the seam's junction its chart copy (SW-V 836: a crossing
/// 1e-5 off a sphere's seam joins the seam's vertex there). A pole, a
/// seam's end, and another face's boundary are other places its arms
/// reach as themselves (YS-V 1302, CC-V 64).
pub(crate) fn on_seam_interior(
    brep: &BRep,
    face: FaceHandle,
    p: &popkern_core::geometry::Point3,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    use popkern_core::geometry::curve::closest_point_in_window;
    let near = 10.0 * tolerance.linear();
    brep.edges.iter().any(|(_, e)| {
        let Some(h) = brep.half_edges.get(e.half_edge) else {
            return false;
        };
        let Some(m) = h.mate().and_then(|m| brep.half_edges.get(m)) else {
            return false;
        };
        let face_of = |l| {
            brep.loops
                .get(l)
                .map(|l: &popkern_core::topology::Loop| l.face)
        };
        if face_of(h.loop_) != Some(face) || face_of(m.loop_) != Some(face) {
            return false;
        }
        let at_vertex = [h.origin, m.origin].into_iter().any(|v| {
            brep.vertices
                .get(v)
                .is_some_and(|v| v.position.distance_to(p) <= near)
        });
        let Some(c) = e.curve.as_ref() else {
            return false;
        };
        let (a, b) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
        let t = closest_point_in_window(&c.geometry, p, a, b);
        !at_vertex && c.geometry.point_at(t).distance_to(p) <= near
    })
}

/// Unit tangent at the piece's endpoint pointing OUTWARD (past the
/// end), by an inside finite difference. Zero vector on degenerate
/// spans — the caller's lateral test then reduces to the plain
/// distance gate.
pub(crate) fn endpoint_outward_tangent(
    piece: &CurveOnSurface,
    end: bool,
) -> popkern_core::geometry::Vec3 {
    let (t0, t1) = piece.curve_3d.parameter_range();
    let span = t1 - t0;
    let (t_at, t_in) = if end {
        (t1, t1 - span * 1e-3)
    } else {
        (t0, t0 + span * 1e-3)
    };
    let d = piece.curve_3d.point_at(t_at) - piece.curve_3d.point_at(t_in);
    d.normalized()
        .unwrap_or(popkern_core::geometry::Vec3::zero())
}
