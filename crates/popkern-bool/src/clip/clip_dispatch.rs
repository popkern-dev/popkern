//! Clip entry dispatch — route a record between the fast Cyrus–Beck
//! line clipper and the UV-space clipper, and hold the cached entry
//! that shares one pierce junction inventory across a face pair's
//! records (PK-511). Split from `clip.rs` under the 500-line budget.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::topology::{BRep, FaceHandle};

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use super::clip_pierce_windows;
use super::clip_windows::clip_curve_uv_space;
use super::{TrimResidue, clip_line_to_face, face_has_curved_boundary, face_polygon_is_reflex};

/// Clip an intersection curve to both face boundaries and produce the
/// trimmed `CurveOnSurface` pieces — one per in-face span (a record
/// can enter and leave a trimmed domain several times; every span is
/// a real intersection edge). Empty when the curve does not pass
/// through both faces.
pub(crate) fn clip_curve_to_faces(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
    residues: Option<&mut Vec<TrimResidue>>,
) -> Vec<CurveOnSurface> {
    // One-shot cache: single-record callers (section, imprint,
    // self-intersect, tests) pay exactly the old cost. The per-pair
    // record loop uses [`clip_curve_to_faces_cached`], which shares
    // the pierce junction inventory across a pair's records
    // (PK-511 — the inventory is record-independent, and recomputing
    // it per wrap record was the intersection stage's ~N^1.8 term).
    let mut cache = clip_pierce_windows::PierceInventoryCache::default();
    clip_curve_to_faces_cached(brep_a, face_a, brep_b, face_b, curve, residues, &mut cache)
}

/// [`clip_curve_to_faces`] with a caller-held pierce junction
/// inventory cache — pass one cache per face pair so every record of
/// the pair reuses the pair-level junction computation.
pub(crate) fn clip_curve_to_faces_cached(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
    residues: Option<&mut Vec<TrimResidue>>,
    cache: &mut clip_pierce_windows::PierceInventoryCache,
) -> Vec<CurveOnSurface> {
    // For non-linear curves, use UV-space clipping.
    if curve.curve_3d.degree != 1 {
        return clip_curve_uv_space(brep_a, face_a, brep_b, face_b, curve, residues, cache);
    }

    // The fast Cyrus–Beck line clipper works on the face's VERTEX
    // polygon — meaningless for faces bounded by closed rims (a
    // cylinder cap has ONE vertex; the periodic wall has two). Such
    // faces route to UV-space clipping, which samples the boundary
    // PCURVES and handles the true trimmed domain (measured: every
    // wall×cap intersection segment of a pocket cut into a swept
    // cylinder was silently discarded here, leaving the tool
    // unsplit and the boolean output corrupt).
    let degenerate_vertex_polygon = |brep: &BRep, face: popkern_core::topology::FaceHandle| {
        brep.face_outer_loop_positions(face).len() < 3
    };
    // A face with INNER loops routes to UV clipping too: the vertex
    // polygon cannot represent holes, so a line lying wholly inside a
    // hole clipped "on-face" and a slot chain imprinted into the hole
    // region (the keyed-bore / vent-slot class, PK-427).
    let has_holes = |brep: &BRep, face: popkern_core::topology::FaceHandle| {
        brep.faces
            .get(face)
            .is_some_and(|f| !f.inner_loops.is_empty())
    };
    if degenerate_vertex_polygon(brep_a, face_a)
        || degenerate_vertex_polygon(brep_b, face_b)
        || has_holes(brep_a, face_a)
        || has_holes(brep_b, face_b)
        || face_polygon_is_reflex(brep_a, face_a)
        || face_polygon_is_reflex(brep_b, face_b)
        || face_has_curved_boundary(brep_a, face_a)
        || face_has_curved_boundary(brep_b, face_b)
    {
        return clip_curve_uv_space(brep_a, face_a, brep_b, face_b, curve, residues, cache);
    }

    let (Some(range_a), Some(range_b)) = (
        clip_line_to_face(brep_a, face_a, &curve.curve_3d),
        clip_line_to_face(brep_b, face_b, &curve.curve_3d),
    ) else {
        return Vec::new();
    };

    // Intersect the two ranges.
    let t_start = range_a.0.max(range_b.0);
    let t_end = range_a.1.min(range_b.1);

    if t_start >= t_end {
        return Vec::new();
    }

    // Build a new degree-1 NURBS for the clipped segment.
    let p_start = curve.curve_3d.point_at(t_start);
    let p_end = curve.curve_3d.point_at(t_end);

    let clipped_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p_start, p_end],
        weights: vec![1.0, 1.0],
    };

    // Rebuild UV curves by evaluating the originals at the clipped parameters.
    let uv1_start = curve.curve_uv1.point_at(t_start);
    let uv1_end = curve.curve_uv1.point_at(t_end);
    let clipped_uv1 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![uv1_start, uv1_end],
        weights: vec![1.0, 1.0],
    };

    let uv2_start = curve.curve_uv2.point_at(t_start);
    let uv2_end = curve.curve_uv2.point_at(t_end);
    let clipped_uv2 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![uv2_start, uv2_end],
        weights: vec![1.0, 1.0],
    };

    vec![CurveOnSurface {
        curve_3d: clipped_3d,
        curve_uv1: clipped_uv1,
        curve_uv2: clipped_uv2,
        is_closed: false,
    }]
}
