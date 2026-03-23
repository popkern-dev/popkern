//! Stage 2: Split faces of both solids along intersection curves.
//!
//! Groups intersection records by face, then calls `BRep::split_face` on each
//! affected face to produce face fragments. Tracks provenance (which original
//! face each fragment came from) and which edges were created by splitting.

use std::collections::{HashMap, HashSet};

use popkern_core::{
    geometry::curve::CurveEval,
    geometry::surface::SurfaceGeometry,
    geometry::{Point3, Vec3},
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle},
};

use crate::clip::clip_curve_to_single_face;
use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::intersect::IntersectionRecord;

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
        // Skip faces that no longer exist (e.g. stale handles).
        if !a.faces.contains_key(face_handle) {
            continue;
        }

        // Skip splitting for non-planar faces — curved faces (sphere, cylinder,
        // etc.) don't support Cyrus-Beck polygon clipping. The face will be
        // classified as a whole unit instead.
        if !is_planar_face(&a, face_handle) {
            continue;
        }

        let record = &records[record_idx];

        // Re-clip the curve to face A's boundary so endpoints land precisely
        // on A's edges (the shared curve may have endpoints on B's boundary).
        let clipped = match clip_curve_to_single_face(&a, face_handle, &record.curve) {
            Some(c) => c,
            None => continue, // curve doesn't intersect this face
        };

        // Ensure boundary vertices exist at the curve 3D endpoints before splitting.
        let (t_start, t_end) = clipped.curve_3d.parameter_range();
        let start_pt = clipped.curve_3d.point_at(t_start);
        let end_pt = clipped.curve_3d.point_at(t_end);
        match ensure_endpoint_vertices(&mut a, face_handle, &[start_pt, end_pt], tolerance) {
            Ok(true) => { /* all endpoints placed */ }
            Ok(false) => continue, // endpoints outside face boundary — skip
            Err(core_err) => {
                return Err(BoolDiagnostics {
                    error: Some(BoolError::SplitFailed {
                        face: face_handle,
                        reason: format!("pre-split edge failed: {core_err}"),
                    }),
                    intersection_curves: records.iter().map(|r| r.curve.clone()).collect(),
                    face_classifications: Vec::new(),
                    failed_stage: Some(Stage::Split),
                });
            }
        }

        // For face A: curve_uv1 is the UV parameterization on surface A.
        // split_face expects (face, trim_curve_uv, trim_curve_3d, tolerance).
        match a.split_face(
            face_handle,
            &clipped.curve_uv1,
            &clipped.curve_3d,
            tolerance,
        ) {
            Ok(result) => {
                fragments_a.origin.insert(result.face_a, face_handle);
                fragments_a.origin.insert(result.face_b, face_handle);
                for &edge in &result.new_edges {
                    fragments_a.intersection_edges.insert(edge);
                }
            }
            Err(core_err) => {
                return Err(BoolDiagnostics {
                    error: Some(BoolError::SplitFailed {
                        face: face_handle,
                        reason: core_err.to_string(),
                    }),
                    intersection_curves: records.iter().map(|r| r.curve.clone()).collect(),
                    face_classifications: Vec::new(),
                    failed_stage: Some(Stage::Split),
                });
            }
        }
    }

    // Split faces in solid B.
    for (&face_handle, &record_idx) in &first_curve_b {
        // Skip faces that no longer exist (e.g. stale handles).
        if !b.faces.contains_key(face_handle) {
            continue;
        }

        // Skip splitting for non-planar faces — curved faces (sphere, cylinder,
        // etc.) don't support Cyrus-Beck polygon clipping. The face will be
        // classified as a whole unit instead.
        if !is_planar_face(&b, face_handle) {
            continue;
        }

        let record = &records[record_idx];

        // Re-clip the curve to face B's boundary so endpoints land precisely
        // on B's edges (the shared curve may have endpoints on A's boundary).
        let clipped = match clip_curve_to_single_face(&b, face_handle, &record.curve) {
            Some(c) => c,
            None => {
                eprintln!("  B split: clip returned None");
                continue;
            }
        };

        // Ensure boundary vertices exist at the curve 3D endpoints before splitting.
        let (t_start, t_end) = clipped.curve_3d.parameter_range();
        let start_pt = clipped.curve_3d.point_at(t_start);
        let end_pt = clipped.curve_3d.point_at(t_end);
        match ensure_endpoint_vertices(&mut b, face_handle, &[start_pt, end_pt], tolerance) {
            Ok(true) => { /* all endpoints placed */ }
            Ok(false) => continue, // endpoints outside face boundary — skip
            Err(core_err) => {
                return Err(BoolDiagnostics {
                    error: Some(BoolError::SplitFailed {
                        face: face_handle,
                        reason: format!("pre-split edge failed: {core_err}"),
                    }),
                    intersection_curves: records.iter().map(|r| r.curve.clone()).collect(),
                    face_classifications: Vec::new(),
                    failed_stage: Some(Stage::Split),
                });
            }
        }

        // For face B: curve_uv2 is the UV parameterization on surface B.
        // It becomes "uv1" relative to face B's surface in the split_face call.
        match b.split_face(
            face_handle,
            &clipped.curve_uv2,
            &clipped.curve_3d,
            tolerance,
        ) {
            Ok(result) => {
                fragments_b.origin.insert(result.face_a, face_handle);
                fragments_b.origin.insert(result.face_b, face_handle);
                for &edge in &result.new_edges {
                    fragments_b.intersection_edges.insert(edge);
                }
            }
            Err(core_err) => {
                return Err(BoolDiagnostics {
                    error: Some(BoolError::SplitFailed {
                        face: face_handle,
                        reason: core_err.to_string(),
                    }),
                    intersection_curves: records.iter().map(|r| r.curve.clone()).collect(),
                    face_classifications: Vec::new(),
                    failed_stage: Some(Stage::Split),
                });
            }
        }
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
