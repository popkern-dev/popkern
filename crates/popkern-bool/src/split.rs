//! Stage 2: Split faces of both solids along intersection curves.
//!
//! Groups intersection records by face, then calls `BRep::split_face` on each
//! affected face to produce face fragments. Tracks provenance (which original
//! face each fragment came from) and which edges were created by splitting.

use std::collections::{HashMap, HashSet};

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
    clip_curve_to_crossings, clip_curve_to_single_face, face_uv_boundary, find_uv_crossings,
};
use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::intersect::IntersectionRecord;
use crate::split_inner_loop::point_inside_face_uv;

/// Provenance map tracking which original face each fragment came from,
/// plus which edges were created by the splitting process.
#[derive(Debug, Clone)]
pub(crate) struct FragmentMap {
    /// Maps each new fragment face handle back to the original face it was split from.
    pub origin: HashMap<FaceHandle, FaceHandle>,
    /// The set of edges created during splitting (trim curve edges).
    pub intersection_edges: HashSet<EdgeHandle>,
}

/// Result of splitting both solids along intersection curves.
#[derive(Debug)]
pub(crate) struct SplitResult {
    /// The modified solid A (with faces split).
    pub a: BRep,
    /// The modified solid B (with faces split).
    pub b: BRep,
    /// Fragment provenance for solid A.
    pub fragments_a: FragmentMap,
    /// Fragment provenance for solid B.
    pub fragments_b: FragmentMap,
}

/// Walk the outer loop of a face, collecting `(half_edge, origin_vertex)` pairs.
fn walk_face_loop(
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
fn point_on_segment(a: &Point3, b: &Point3, point: &Point3, tol: f64) -> Option<f64> {
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

/// Ensure that vertices exist on the face boundary at the given 3D points.
///
/// For each endpoint, if no existing boundary vertex is within tolerance, the
/// boundary edge containing the point is split via `split_edge` to create one.
/// This must be called before `split_face` when intersection curve endpoints
/// land on edges rather than at existing vertices.
///
/// Returns `true` if all endpoints were successfully matched (either to an
/// existing vertex or by splitting an edge). Returns `false` if any endpoint
/// could not be located on the face boundary.
fn ensure_endpoint_vertices(
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

/// Returns `true` if the face has planar surface geometry (or no geometry at all).
///
/// Non-planar faces (sphere, cylinder, etc.) are not supported by the
/// Cyrus-Beck polygon clipping used during splitting, so they must be
/// skipped and classified as whole faces instead.
fn is_planar_face(brep: &BRep, face: FaceHandle) -> bool {
    match brep.faces.get(face) {
        Some(f) => match &f.surface {
            Some(sof) => matches!(sof.geometry, SurfaceGeometry::Plane { .. }),
            None => true, // No geometry attached — treat as planar (legacy)
        },
        None => true, // Face doesn't exist — caller will skip anyway
    }
}

/// Returns `true` if the face's outer loop half-edges all have pcurves.
fn face_has_pcurves(brep: &BRep, face: FaceHandle) -> bool {
    face_uv_boundary(brep, face).is_some()
}

/// Build a [`BoolDiagnostics`] for a split-stage failure.
fn split_err(face: FaceHandle, reason: String, records: &[IntersectionRecord]) -> BoolDiagnostics {
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
fn try_uv_split(
    brep: &mut BRep,
    face_handle: FaceHandle,
    curve_3d: &NurbsCurve,
    curve_uv: &NurbsCurve,
    tolerance: &Tolerance,
    fragments: &mut FragmentMap,
    records: &[IntersectionRecord],
) -> Result<bool, BoolDiagnostics> {
    let crossings = match find_uv_crossings(brep, face_handle, curve_uv, tolerance) {
        Some(c) => c,
        None => return Ok(false),
    };

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
        match brep.split_face(face_handle, &clipped_uv, &clipped_3d, tolerance) {
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

    // 0 crossings — check for closed inner-loop case.
    if crossings.is_empty() {
        let (t_start, t_end) = curve_3d.parameter_range();
        let start_pt = curve_3d.point_at(t_start);
        let end_pt = curve_3d.point_at(t_end);
        if tolerance.points_coincident(&start_pt, &end_pt) {
            let mid_uv = curve_uv.point_at((t_start + t_end) / 2.0);
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
    }

    Ok(false)
}

/// Try the planar Cyrus-Beck clipping path for a face.
///
/// Returns `true` if the face was split, `false` if it should be skipped.
fn try_planar_split(
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
fn split_one_face(
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

/// Split the faces of two BRep solids along intersection curves.
///
/// For each intersection record, splits the corresponding face of A using
/// `curve_uv1` (the UV curve on surface A) and the face of B using
/// `curve_uv2` (the UV curve on surface B).
///
/// Only the first intersection curve per face is applied. Faces with no
/// intersection curves map to themselves in the fragment map.
///
/// # Errors
///
/// Returns [`BoolDiagnostics`] if any `split_face` call fails.
pub(crate) fn split_along_curves(
    mut a: BRep,
    mut b: BRep,
    records: &[IntersectionRecord],
    tolerance: &Tolerance,
) -> Result<SplitResult, BoolDiagnostics> {
    let mut fragments_a = FragmentMap {
        origin: HashMap::new(),
        intersection_edges: HashSet::new(),
    };
    let mut fragments_b = FragmentMap {
        origin: HashMap::new(),
        intersection_edges: HashSet::new(),
    };

    // Group records by face_a and face_b (take only the first curve per face).
    let mut first_curve_a: HashMap<FaceHandle, usize> = HashMap::new();
    let mut first_curve_b: HashMap<FaceHandle, usize> = HashMap::new();

    for (idx, record) in records.iter().enumerate() {
        first_curve_a.entry(record.face_a).or_insert(idx);
        first_curve_b.entry(record.face_b).or_insert(idx);
    }

    // Split faces in solid A.
    for (&face_handle, &record_idx) in &first_curve_a {
        if !a.faces.contains_key(face_handle) {
            continue;
        }
        let record = &records[record_idx];
        split_one_face(
            &mut a,
            face_handle,
            &record.curve.curve_3d,
            &record.curve.curve_uv1,
            &record.curve,
            tolerance,
            records,
            &mut fragments_a,
        )?;
    }

    // Split faces in solid B.
    for (&face_handle, &record_idx) in &first_curve_b {
        if !b.faces.contains_key(face_handle) {
            continue;
        }
        let record = &records[record_idx];
        split_one_face(
            &mut b,
            face_handle,
            &record.curve.curve_3d,
            &record.curve.curve_uv2,
            &record.curve,
            tolerance,
            records,
            &mut fragments_b,
        )?;
    }

    // All unsplit faces in A map to themselves.
    for face_handle in a.faces.keys() {
        fragments_a.origin.entry(face_handle).or_insert(face_handle);
    }

    // All unsplit faces in B map to themselves.
    for face_handle in b.faces.keys() {
        fragments_b.origin.entry(face_handle).or_insert(face_handle);
    }

    Ok(SplitResult {
        a,
        b,
        fragments_a,
        fragments_b,
    })
}

#[cfg(test)]
#[path = "split_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "split_extra_tests.rs"]
mod extra_tests;
