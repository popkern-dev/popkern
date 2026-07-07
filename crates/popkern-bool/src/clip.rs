//! Curve clipping to face boundaries.
//!
//! Clips intersection curves to the polygonal boundaries of BRep faces.
//! Used by the intersection stage to trim curves that extend beyond face
//! boundaries due to the UV-box clipping limitations for non-orthogonal
//! plane parameterizations.
//!
//! Line (degree-1) clipping is in the [`clip_line`] submodule.

mod clip_line;
mod clip_uv;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

pub(crate) use clip_line::{clip_curve_to_single_face, clip_line_to_face};
pub(crate) use clip_uv::{
    UvCrossing, clip_curve_to_crossings, clip_curve_uv_space, find_uv_crossings,
};
#[cfg(test)]
use clip_uv::{is_near_seam, shift_curve_u, uv_range_from_crossings};

/// Stage 2b.58: extract the underlying `CurveOnSurface` (NURBS
/// triple) from an `EdgeTriple` produced by
/// `compute_edge_geometry_unanchored`. Returns `None` if the
/// triple's curves aren't NURBS-backed — shouldn't happen for the
/// marching path the unanchored API uses, but we guard against
/// future analytic-curve paths defensively.
///
/// Bridges the EdgeTriple-shape returned by the new factory with
/// the `CurveOnSurface`-shape that `clip_curve_to_faces` and the
/// rest of `popkern-bool`'s pipeline still consume.
pub(crate) fn edge_triple_to_curve_on_surface(
    triple: &popkern_core::geometry::edge_triple::EdgeTriple,
) -> Option<CurveOnSurface> {
    use popkern_core::geometry::pcurve::PcurveGeometry;
    let CurveGeometry::Nurbs(curve_3d) = &triple.curve_3d else {
        return None;
    };
    let PcurveGeometry::Nurbs(curve_uv1) = &triple.pcurve_a.geometry else {
        return None;
    };
    let PcurveGeometry::Nurbs(curve_uv2) = &triple.pcurve_b.geometry else {
        return None;
    };
    let is_closed = matches!(
        triple.topology,
        popkern_core::geometry::edge_triple::EdgeTopology::Closed { .. }
    );
    Some(CurveOnSurface {
        curve_3d: curve_3d.clone(),
        curve_uv1: curve_uv1.clone(),
        curve_uv2: curve_uv2.clone(),
        is_closed,
    })
}

/// Collect vertex positions from the outer loop of a face by walking half-edges.
///
/// Returns an empty `Vec` if the face's loop has no valid half-edges.
/// Uses a bounded iteration limit to prevent infinite loops on broken topology.
pub(crate) fn face_vertex_positions(brep: &BRep, face_handle: FaceHandle) -> Vec<Point3> {
    let face = &brep.faces[face_handle];
    let loop_handle = face.outer_loop;
    let loop_data = match brep.loops.get(loop_handle) {
        Some(l) => l,
        None => return Vec::new(),
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut positions = Vec::new();
    let mut current = start_he;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = &brep.half_edges[current];
        positions.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}

/// Extract the UV boundary of a face from its half-edge pcurves.
///
/// Returns a list of `(half_edge_handle, pcurve)` forming the outer loop
/// boundary in UV space. Returns `None` if any half-edge lacks a pcurve.
pub(crate) fn face_uv_boundary(
    brep: &BRep,
    face: FaceHandle,
) -> Option<Vec<(HalfEdgeHandle, popkern_core::geometry::pcurve::Pcurve)>> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut boundary = Vec::new();
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = brep.half_edges.get(cursor)?;
        let pcurve = he.pcurve.as_ref()?.clone();
        boundary.push((cursor, pcurve));
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    Some(boundary)
}

/// Clip an intersection curve to both face boundaries and produce a trimmed
/// `CurveOnSurface`. Returns `None` if the curve does not pass through both faces.
pub(crate) fn clip_curve_to_faces(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
) -> Option<CurveOnSurface> {
    // For non-linear curves, use UV-space clipping.
    if curve.curve_3d.degree != 1 {
        return clip_curve_uv_space(brep_a, face_a, brep_b, face_b, curve);
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
    if degenerate_vertex_polygon(brep_a, face_a) || degenerate_vertex_polygon(brep_b, face_b) {
        return clip_curve_uv_space(brep_a, face_a, brep_b, face_b, curve);
    }

    let range_a = clip_line_to_face(brep_a, face_a, &curve.curve_3d)?;
    let range_b = clip_line_to_face(brep_b, face_b, &curve.curve_3d)?;

    // Intersect the two ranges.
    let t_start = range_a.0.max(range_b.0);
    let t_end = range_a.1.min(range_b.1);

    if t_start >= t_end {
        return None;
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

    Some(CurveOnSurface {
        curve_3d: clipped_3d,
        curve_uv1: clipped_uv1,
        curve_uv2: clipped_uv2,
        is_closed: false,
    })
}

#[cfg(test)]
#[path = "clip_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "clip_uv_tests.rs"]
mod uv_tests;
