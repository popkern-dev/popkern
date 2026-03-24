//! Helper functions for face splitting along intersection curves.
//!
//! Extracted from `split.rs` to keep file sizes under the 500-line limit.

use popkern_core::{
    geometry::curve::CurveEval,
    geometry::curve_on_surface::CurveOnSurface,
    geometry::nurbs_curve::NurbsCurve,
    geometry::surface::SurfaceGeometry,
    geometry::{Point3, Vec3},
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle},
};

use crate::clip::{
    UvCrossing, clip_curve_to_crossings, clip_curve_to_single_face, face_uv_boundary,
    find_uv_crossings,
};
use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::intersect::IntersectionRecord;
use crate::split_inner_loop::point_inside_face_uv;

use super::FragmentMap;

/// Walk the outer loop of a face, collecting `(half_edge, origin_vertex)` pairs.
pub(super) fn walk_face_loop(
    brep: &BRep,
    face: FaceHandle,
) -> Result<Vec<(HalfEdgeHandle, VertexHandle)>, popkern_core::error::CoreError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| popkern_core::error::CoreError::InvalidHandle("face".into()))?;
    let loop_data = brep
        .loops
        .get(face_data.outer_loop)
        .ok_or_else(|| popkern_core::error::CoreError::InvalidHandle("loop".into()))?;
    let start = loop_data.half_edge;
    let mut pairs = Vec::new();
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = &brep.half_edges[cursor];
        pairs.push((cursor, he.origin));
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    Ok(pairs)
}

/// Compute the parameter `t` at which `point` projects onto the line segment
/// from `a` to `b`. Returns `None` if the point is not within `tol` of the
/// segment or if `t` is not in `(eps, 1-eps)` (i.e. already at an endpoint).
pub(super) fn point_on_segment(a: &Point3, b: &Point3, point: &Point3, tol: f64) -> Option<f64> {
    let ab = *b - *a;
    let len_sq = ab.length_squared();
    if len_sq < tol * tol {
        return None; // degenerate edge
    }
    let ap = *point - *a;
    let t = ap.dot(&ab) / len_sq;
    // Must be strictly interior — if close to 0 or 1, the endpoint already exists.
    let eps = 1e-8;
    if t <= eps || t >= 1.0 - eps {
        return None;
    }
    let closest = *a + Vec3::new(ab.x * t, ab.y * t, ab.z * t);
    if point.distance_to(&closest) < tol {
        Some(t)
    } else {
        None
    }
}

/// Ensure vertices exist on the face boundary at 3D points (planar fallback).
pub(super) fn ensure_endpoint_vertices(
    brep: &mut BRep,
    face: FaceHandle,
    endpoints: &[Point3],
    tolerance: &Tolerance,
) -> Result<bool, popkern_core::error::CoreError> {
    for endpoint in endpoints {
        // Re-walk the loop each time because split_edge mutates topology.
        let pairs = walk_face_loop(brep, face)?;

        // Check if any existing vertex already matches.
        let already_exists = pairs
            .iter()
            .any(|&(_, vh)| tolerance.points_coincident(&brep.vertices[vh].position, endpoint));
        if already_exists {
            continue;
        }

        // Find which edge the endpoint lies on.
        // Compute t relative to the EDGE's own parameterization
        // (from edge.half_edge.origin to its twin.origin), not the loop
        // half-edge direction which may be reversed.
        let mut best: Option<(EdgeHandle, f64)> = None;
        for &(he_handle, _) in &pairs {
            let he = &brep.half_edges[he_handle];
            let origin_pos = brep.vertices[he.origin].position;
            let next_he = &brep.half_edges[he.next];
            let end_pos = brep.vertices[next_he.origin].position;

            if let Some(_t_loop) =
                point_on_segment(&origin_pos, &end_pos, endpoint, tolerance.linear)
            {
                // Compute t relative to the edge's canonical direction.
                let edge_he = brep.edges[he.edge].half_edge;
                let edge_origin = brep.vertices[brep.half_edges[edge_he].origin].position;
                let edge_twin = brep.half_edges[edge_he].twin;
                let edge_end = brep.vertices[brep.half_edges[edge_twin].origin].position;
                if let Some(t_edge) =
                    point_on_segment(&edge_origin, &edge_end, endpoint, tolerance.linear)
                {
                    best = Some((he.edge, t_edge));
                }
                break;
            }
        }

        if let Some((edge, t)) = best {
            brep.split_edge(edge, t)?;
        } else {
            // Endpoint not on any boundary edge — curve extends beyond this face.
            return Ok(false);
        }
    }
    Ok(true)
}

/// Place vertices on the face boundary using UV crossing edge identity.
pub(super) fn ensure_crossing_vertices(
    brep: &mut BRep,
    face: FaceHandle,
    crossings: &[&UvCrossing],
    curve_3d: &NurbsCurve,
    tolerance: &Tolerance,
) -> Result<bool, popkern_core::error::CoreError> {
    for crossing in crossings {
        let he_handle = crossing.boundary_he;
        let crossing_3d = curve_3d.point_at(crossing.t_curve);
        // Check if vertex already exists.
        let pairs = walk_face_loop(brep, face)?;
        if pairs
            .iter()
            .any(|&(_, vh)| tolerance.points_coincident(&brep.vertices[vh].position, &crossing_3d))
        {
            continue;
        }
        // Check half-edge still exists (previous split may have changed topology).
        if !brep.half_edges.contains_key(he_handle) {
            let pairs2 = walk_face_loop(brep, face)?;
            if pairs2.iter().any(|&(_, vh)| {
                tolerance.points_coincident(&brep.vertices[vh].position, &crossing_3d)
            }) {
                continue;
            }
            return Ok(false);
        }
        let edge_handle = brep.half_edges[he_handle].edge;
        if !brep.edges.contains_key(edge_handle) {
            return Ok(false);
        }
        let canonical_he = brep.edges[edge_handle].half_edge;
        let t_edge = if canonical_he == he_handle {
            crossing.t_boundary
        } else {
            1.0 - crossing.t_boundary
        };
        let eps = 1e-8;
        if t_edge <= eps || t_edge >= 1.0 - eps {
            continue;
        }
        brep.split_edge(edge_handle, t_edge)?;
    }
    Ok(true)
}

/// Returns `true` if the face has planar surface geometry (or no geometry at all).
///
/// Non-planar faces (sphere, cylinder, etc.) are not supported by the
/// Cyrus-Beck polygon clipping used during splitting, so they must be
/// skipped and classified as whole faces instead.
pub(super) fn is_planar_face(brep: &BRep, face: FaceHandle) -> bool {
    match brep.faces.get(face) {
        Some(f) => match &f.surface {
            Some(sof) => matches!(sof.geometry, SurfaceGeometry::Plane { .. }),
            None => true, // No geometry attached — treat as planar (legacy)
        },
        None => true, // Face doesn't exist — caller will skip anyway
    }
}

/// Returns `true` if the face's outer loop half-edges all have pcurves.
pub(super) fn face_has_pcurves(brep: &BRep, face: FaceHandle) -> bool {
    face_uv_boundary(brep, face).is_some()
}

/// Build a [`BoolDiagnostics`] for a split-stage failure.
pub(super) fn split_err(
    face: FaceHandle,
    reason: String,
    records: &[IntersectionRecord],
) -> BoolDiagnostics {
    BoolDiagnostics {
        error: Some(BoolError::SplitFailed { face, reason }),
        intersection_curves: records.iter().map(|r| r.curve.clone()).collect(),
        face_classifications: Vec::new(),
        failed_stage: Some(Stage::Split),
    }
}

/// Try UV-space splitting for a curved face: find where the intersection
/// curve crosses the face boundary in UV space, clip to interior, and split.
///
/// Returns `true` if the face was split, `false` if it should be skipped.
pub(super) fn try_uv_split(
    brep: &mut BRep,
    face_handle: FaceHandle,
    curve_3d: &NurbsCurve,
    curve_uv: &NurbsCurve,
    tolerance: &Tolerance,
    fragments: &mut FragmentMap,
    records: &[IntersectionRecord],
) -> Result<bool, BoolDiagnostics> {
    // Check if the original curve is closed (start ≈ end).
    let (t0_orig, t1_orig) = curve_3d.parameter_range();
    let is_closed =
        tolerance.points_coincident(&curve_3d.point_at(t0_orig), &curve_3d.point_at(t1_orig));

    let crossings = match find_uv_crossings(brep, face_handle, curve_uv, tolerance) {
        Some(c) => c,
        None => return Ok(false),
    };

    // Sub-problem A: curve crosses boundary (>= 2 crossings).
    // Clip to the interior segment and split the face.
    if crossings.len() >= 2 {
        let t_enter = crossings[0].t_curve;
        let t_exit = crossings[crossings.len() - 1].t_curve;
        let clipped_3d = match clip_curve_to_crossings(curve_3d, t_enter, t_exit) {
            Some(c) => c,
            None => return Ok(false),
        };
        let clipped_uv = match clip_curve_to_crossings(curve_uv, t_enter, t_exit) {
            Some(c) => c,
            None => return Ok(false),
        };
        let (ts, te) = clipped_3d.parameter_range();
        let start_pt = clipped_3d.point_at(ts);
        let end_pt = clipped_3d.point_at(te);
        // Use UV crossing edge identity for all face types. With orthonormal
        // plane UV axes (from plane_from_points), the pcurve parameter maps
        // exactly to the edge parameter for both planar and curved faces.
        let endpoints_ok = ensure_crossing_vertices(
            brep,
            face_handle,
            &[&crossings[0], &crossings[crossings.len() - 1]],
            curve_3d,
            tolerance,
        );
        match endpoints_ok {
            Ok(true) => {
                // Rebuild the clipped 3D curve endpoints from the actual vertex
                // positions on the face boundary (which may differ slightly from
                // curve_3d.point_at(t) due to NURBS evaluation precision).
                let pairs = walk_face_loop(brep, face_handle)
                    .map_err(|e| split_err(face_handle, e.to_string(), records))?;
                let snap_to_vertex = |pt: Point3| -> Point3 {
                    pairs
                        .iter()
                        .map(|&(_, vh)| brep.vertices[vh].position)
                        .min_by(|a, b| {
                            a.distance_to(&pt)
                                .partial_cmp(&b.distance_to(&pt))
                                .unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .unwrap_or(pt)
                };
                let snapped_start = snap_to_vertex(start_pt);
                let snapped_end = snap_to_vertex(end_pt);
                let snapped_3d = NurbsCurve {
                    degree: 1,
                    knots: vec![0.0, 0.0, 1.0, 1.0],
                    control_points: vec![snapped_start, snapped_end],
                    weights: vec![1.0, 1.0],
                };
                match brep.split_face(face_handle, &clipped_uv, &snapped_3d, tolerance) {
                    Ok(result) => {
                        fragments.origin.insert(result.face_a, face_handle);
                        fragments.origin.insert(result.face_b, face_handle);
                        for &edge in &result.new_edges {
                            fragments.intersection_edges.insert(edge);
                        }
                        return Ok(true);
                    }
                    Err(e) => return Err(split_err(face_handle, e.to_string(), records)),
                }
            }
            Ok(false) if is_closed => {
                // Endpoints not on boundary but curve is closed — the UV
                // crossings are spurious (near-tangent). Fall through to
                // the inner-loop path below.
            }
            Ok(false) => return Ok(false),
            Err(e) => {
                return Err(split_err(
                    face_handle,
                    format!("pre-split edge failed: {e}"),
                    records,
                ));
            }
        }
    }

    // Sub-problem B: closed curve entirely inside face (0 real crossings).
    // Also reached when crossings were spurious (closed curve that fell through above).
    if is_closed {
        // Use the 3D midpoint projected to the face's UV space, rather than
        // evaluating the UV curve directly. For non-orthogonal UV axes, the
        // UV curve's midpoint may fall outside the face's UV polygon even when
        // the 3D point is inside the face.
        let mid_3d = curve_3d.point_at((t0_orig + t1_orig) / 2.0);
        let mid_uv = if let Some(sof) = brep.faces[face_handle].surface.as_ref() {
            popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &mid_3d)
                .map(|(u, v)| Point3::new(u, v, 0.0))
                .unwrap_or_else(|| curve_uv.point_at((t0_orig + t1_orig) / 2.0))
        } else {
            curve_uv.point_at((t0_orig + t1_orig) / 2.0)
        };
        if point_inside_face_uv(brep, face_handle, mid_uv, tolerance) {
            let n = (curve_3d.control_points.len() * 2).max(12);
            match brep.add_inner_loop(face_handle, curve_3d, curve_uv, n, tolerance) {
                Ok(result) => {
                    for &edge in &result.edges {
                        fragments.intersection_edges.insert(edge);
                    }
                    return Ok(true);
                }
                Err(e) => {
                    let msg = format!("add_inner_loop failed: {e}");
                    return Err(split_err(face_handle, msg, records));
                }
            }
        }
    }

    Ok(false)
}

/// Try the planar Cyrus-Beck clipping path for a face.
///
/// Returns `true` if the face was split, `false` if it should be skipped.
pub(super) fn try_planar_split(
    brep: &mut BRep,
    face_handle: FaceHandle,
    curve_uv: &NurbsCurve,
    full_curve: &CurveOnSurface,
    tolerance: &Tolerance,
    fragments: &mut FragmentMap,
    records: &[IntersectionRecord],
) -> Result<bool, BoolDiagnostics> {
    let clipped = match clip_curve_to_single_face(brep, face_handle, full_curve) {
        Some(c) => c,
        None => return Ok(false),
    };

    let (t_start, t_end) = clipped.curve_3d.parameter_range();
    let start_pt = clipped.curve_3d.point_at(t_start);
    let end_pt = clipped.curve_3d.point_at(t_end);
    match ensure_endpoint_vertices(brep, face_handle, &[start_pt, end_pt], tolerance) {
        Ok(true) => {}
        Ok(false) => return Ok(false),
        Err(e) => {
            return Err(split_err(
                face_handle,
                format!("pre-split edge failed: {e}"),
                records,
            ));
        }
    }

    // Determine which UV curve to use: for face A it's curve_uv1, for face B
    // it's curve_uv2. Compare the input curve_uv with the record's uv1 vs uv2.
    let use_uv1 = curve_uv.control_points.first().map(|p| p.x)
        == full_curve.curve_uv1.control_points.first().map(|p| p.x);
    let trim_uv = if use_uv1 {
        &clipped.curve_uv1
    } else {
        &clipped.curve_uv2
    };

    match brep.split_face(face_handle, trim_uv, &clipped.curve_3d, tolerance) {
        Ok(result) => {
            fragments.origin.insert(result.face_a, face_handle);
            fragments.origin.insert(result.face_b, face_handle);
            for &edge in &result.new_edges {
                fragments.intersection_edges.insert(edge);
            }
            Ok(true)
        }
        Err(e) => Err(split_err(face_handle, e.to_string(), records)),
    }
}

/// Split a single face along an intersection curve using UV-space or planar
/// clipping depending on the face's geometry.
///
/// Tries UV-space splitting first (for faces with pcurves), then falls back
/// to Cyrus-Beck planar clipping for planar faces without pcurves.
#[allow(clippy::too_many_arguments)]
pub(super) fn split_one_face(
    brep: &mut BRep,
    face_handle: FaceHandle,
    curve_3d: &NurbsCurve,
    curve_uv: &NurbsCurve,
    full_curve: &CurveOnSurface,
    tolerance: &Tolerance,
    records: &[IntersectionRecord],
    fragments: &mut FragmentMap,
) -> Result<(), BoolDiagnostics> {
    if face_has_pcurves(brep, face_handle) {
        // UV-space splitting for faces with pcurves (curved or planar).
        if try_uv_split(
            brep,
            face_handle,
            curve_3d,
            curve_uv,
            tolerance,
            fragments,
            records,
        )? {
            return Ok(());
        }
        // UV split didn't produce crossings — try planar fallback if applicable.
    }

    if is_planar_face(brep, face_handle) {
        // Planar Cyrus-Beck clipping path (fallback).
        try_planar_split(
            brep,
            face_handle,
            curve_uv,
            full_curve,
            tolerance,
            fragments,
            records,
        )?;
    }
    // Non-planar face without pcurves: skip (classified as whole face).

    Ok(())
}

/// Find which fragment face (descended from `original_face`) contains the
/// midpoint of a UV curve, using the UV point-in-face test.
///
/// Returns the fragment face handle, or `None` if no fragment contains the
/// point (e.g. the curve lies outside all fragments).
pub(super) fn find_fragment_for_curve(
    brep: &BRep,
    original_face: FaceHandle,
    curve_uv: &NurbsCurve,
    fragments: &FragmentMap,
    tolerance: &Tolerance,
) -> Option<FaceHandle> {
    let (t0, t1) = curve_uv.parameter_range();
    let mid_uv = curve_uv.point_at((t0 + t1) / 2.0);

    // Collect all fragment faces that descended from `original_face`.
    let candidates: Vec<FaceHandle> = fragments
        .origin
        .iter()
        .filter(|&(_, &orig)| orig == original_face)
        .map(|(&frag, _)| frag)
        .collect();

    for candidate in &candidates {
        if !brep.faces.contains_key(*candidate) {
            continue;
        }
        if point_inside_face_uv(brep, *candidate, mid_uv, tolerance) {
            return Some(*candidate);
        }
    }

    // Fallback: if the original face still exists and no fragment matched,
    // try the original face itself.
    if brep.faces.contains_key(original_face) {
        return Some(original_face);
    }

    None
}
