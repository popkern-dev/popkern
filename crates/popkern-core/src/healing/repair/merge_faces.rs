//! Coplanar face merging: detect and merge adjacent faces on the same surface.
//!
//! Two faces are coplanar candidates if they share an edge and their
//! surface normals are parallel at the shared edge midpoint.  For planar
//! faces an additional distance-from-origin check ensures the planes
//! are truly coincident (not just parallel).

use std::collections::HashSet;

use crate::geometry::surface::SurfaceEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};

/// Result of a face merge operation.
#[derive(Debug, Clone)]
pub struct FaceMergeResult {
    /// Number of coplanar face pairs detected.
    pub pairs_detected: usize,
    /// Number of face pairs actually merged (topology modified).
    /// Currently always 0; full face merging is deferred to a future
    /// iteration that implements the required Euler operator sequences.
    pub pairs_merged: usize,
}

/// Detect pairs of adjacent coplanar faces.
///
/// Two faces are coplanar if they share an edge and their surface normals
/// are parallel (angle < `tolerance.angular()`) at the shared edge midpoint,
/// and their surface types match.  For planes, an additional check verifies
/// they are truly coplanar (same normal direction **and** same signed
/// distance from the origin within `tolerance.linear()`).
///
/// Returns deduplicated pairs sorted by handle order.
pub fn detect_coplanar_pairs(brep: &BRep, tolerance: &Tolerance) -> Vec<(FaceHandle, FaceHandle)> {
    let angular_tol = tolerance.angular();
    let linear_tol = tolerance.linear();
    let mut seen = HashSet::new();

    for (_edge_handle, edge_data) in &brep.edges {
        // Get the two half-edges for this edge.
        let Some(he) = brep.half_edges.get(edge_data.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };

        // Get the faces from the half-edge loops.
        let Some(loop_a) = brep.loops.get(he.loop_) else {
            continue;
        };
        let Some(loop_b) = brep.loops.get(twin.loop_) else {
            continue;
        };

        let face_a = loop_a.face;
        let face_b = loop_b.face;

        // Skip seam edges (same face on both sides).
        if face_a == face_b {
            continue;
        }

        // Skip if already seen (canonical ordering for deduplication).
        let pair = canonical_pair(face_a, face_b);
        if seen.contains(&pair) {
            continue;
        }

        // Both faces must have surface geometry.
        let Some(face_a_data) = brep.faces.get(face_a) else {
            continue;
        };
        let Some(sof_a) = face_a_data.surface.as_ref() else {
            continue;
        };

        let Some(face_b_data) = brep.faces.get(face_b) else {
            continue;
        };
        let Some(sof_b) = face_b_data.surface.as_ref() else {
            continue;
        };

        // Compute a sample point on the shared edge (midpoint of the two
        // endpoint vertex positions).
        let start_pos = brep.vertices.get(he.origin).map(|v| v.position);
        let end_pos = brep.vertices.get(twin.origin).map(|v| v.position);

        let Some(sp) = start_pos else { continue };
        let Some(ep) = end_pos else { continue };

        let mid_3d = crate::geometry::Point3::new(
            (sp.x + ep.x) * 0.5,
            (sp.y + ep.y) * 0.5,
            (sp.z + ep.z) * 0.5,
        );

        // Project the midpoint onto each surface to get (u, v) parameters,
        // then evaluate normals there.
        let Some((ua, va)) = sof_a.geometry.closest_point(&mid_3d) else {
            continue;
        };
        let Some((ub, vb)) = sof_b.geometry.closest_point(&mid_3d) else {
            continue;
        };

        let normal_a = sof_a.geometry.normal_at(ua, va);
        let normal_b = sof_b.geometry.normal_at(ub, vb);

        // Normalise; skip if either is degenerate.
        let Some(unit_a) = normal_a.normalized() else {
            continue;
        };
        let Some(unit_b) = normal_b.normalized() else {
            continue;
        };

        // Check normals are parallel (or anti-parallel).
        let angle = unit_a.angle_to(&unit_b);
        let deviation = angle.min(std::f64::consts::PI - angle);
        if deviation > angular_tol {
            continue;
        }

        // For planes, additionally check they are coplanar (same signed
        // distance from origin).
        if !surfaces_coplanar(&sof_a.geometry, &sof_b.geometry, linear_tol) {
            continue;
        }

        // Suppress cross-type matches: Plane vs Cylinder should not merge.
        if !surface_types_compatible(&sof_a.geometry, &sof_b.geometry) {
            continue;
        }

        seen.insert(pair);
    }

    let mut pairs: Vec<(FaceHandle, FaceHandle)> = seen.into_iter().collect();
    pairs.sort_by(|a, b| {
        let ka = (a.0, a.1);
        let kb = (b.0, b.1);
        ka.partial_cmp(&kb).unwrap_or(std::cmp::Ordering::Equal)
    });
    pairs
}

/// Merge detected coplanar face pairs by removing shared edges.
///
/// For each pair, attempts to remove the shared edge, merging the two
/// faces into one.  This simplifies topology without changing geometry.
/// Only pairs that share exactly one edge are merged; complex multi-edge
/// boundaries are skipped.
pub fn merge_coplanar_faces(brep: &mut BRep, tolerance: &Tolerance) -> FaceMergeResult {
    let pairs = detect_coplanar_pairs(brep, tolerance);
    let pairs_detected = pairs.len();
    let mut pairs_merged = 0;

    for (face_a, face_b) in &pairs {
        // Skip pairs where one face has already been removed by
        // a prior merge in this batch.
        if !brep.faces.contains_key(*face_a) || !brep.faces.contains_key(*face_b) {
            continue;
        }
        if attempt_face_merge(brep, *face_a, *face_b).is_ok() {
            pairs_merged += 1;
        }
    }

    FaceMergeResult {
        pairs_detected,
        pairs_merged,
    }
}

/// Attempt to merge two faces by removing their shared edge.
///
/// This is a simplified merge that works when the two faces share
/// exactly one edge. It reassigns the second face's half-edges to
/// the first face's loop and removes the shared edge and second face.
///
/// Returns `Err` if the faces don't share exactly one edge or if
/// the topology is too complex for this simplified merge.
fn attempt_face_merge(
    brep: &mut BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
) -> Result<(), crate::error::CoreError> {
    use crate::error::CoreError;
    use crate::topology::{EdgeHandle, HalfEdgeHandle};

    // Find shared edges: edges where one half-edge's face == face_a
    // and the twin's face == face_b.
    let mut shared_edges: Vec<(EdgeHandle, HalfEdgeHandle, HalfEdgeHandle)> = Vec::new();

    for (edge_handle, edge_data) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge_data.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };

        let Some(loop_he) = brep.loops.get(he.loop_) else {
            continue;
        };
        let Some(loop_twin) = brep.loops.get(twin.loop_) else {
            continue;
        };

        let (f_he, f_twin) = (loop_he.face, loop_twin.face);

        if (f_he == face_a && f_twin == face_b) || (f_he == face_b && f_twin == face_a) {
            // Ensure he_a is on face_a, he_b is on face_b.
            let (he_a, he_b) = if f_he == face_a {
                (edge_data.half_edge, he.twin)
            } else {
                (he.twin, edge_data.half_edge)
            };
            shared_edges.push((edge_handle, he_a, he_b));
        }
    }

    // Only handle the simple case: exactly one shared edge.
    if shared_edges.len() != 1 {
        return Err(CoreError::PreconditionFailed(format!(
            "face merge requires exactly 1 shared edge, found {}",
            shared_edges.len()
        )));
    }

    let (shared_edge, he_a, he_b) = shared_edges[0];

    // Read link pointers before mutation.
    let he_a_prev = brep
        .half_edges
        .get(he_a)
        .ok_or_else(|| CoreError::InvalidHandle("he_a".into()))?
        .prev;
    let he_a_next = brep
        .half_edges
        .get(he_a)
        .ok_or_else(|| CoreError::InvalidHandle("he_a".into()))?
        .next;
    let he_b_prev = brep
        .half_edges
        .get(he_b)
        .ok_or_else(|| CoreError::InvalidHandle("he_b".into()))?
        .prev;
    let he_b_next = brep
        .half_edges
        .get(he_b)
        .ok_or_else(|| CoreError::InvalidHandle("he_b".into()))?
        .next;

    // Relink: bypass the shared edge in both loops.
    if let Some(prev) = brep.half_edges.get_mut(he_a_prev) {
        prev.next = he_b_next;
    }
    if let Some(next) = brep.half_edges.get_mut(he_b_next) {
        next.prev = he_a_prev;
    }
    if let Some(prev) = brep.half_edges.get_mut(he_b_prev) {
        prev.next = he_a_next;
    }
    if let Some(next) = brep.half_edges.get_mut(he_a_next) {
        next.prev = he_b_prev;
    }

    // Get face_a's outer loop and face_b's outer loop before removal.
    let loop_a = brep
        .faces
        .get(face_a)
        .ok_or_else(|| CoreError::InvalidHandle("face_a".into()))?
        .outer_loop;
    let loop_b = brep
        .faces
        .get(face_b)
        .ok_or_else(|| CoreError::InvalidHandle("face_b".into()))?
        .outer_loop;

    // Reassign all of face_b's loop half-edges to face_a's loop.
    let mut current = he_b_next;
    let mut guard = 0usize;
    let max_steps = brep.half_edges.len();
    while let Some(he) = brep.half_edges.get_mut(current) {
        he.loop_ = loop_a;
        current = he.next;
        guard += 1;
        if current == he_b_next || guard > max_steps {
            break;
        }
    }

    // Update the loop_a half_edge pointer in case it pointed to the
    // removed half-edge.
    if let Some(la) = brep.loops.get_mut(loop_a)
        && la.half_edge == he_a
    {
        la.half_edge = he_b_next;
    }

    // Remove the shared edge and its two half-edges.
    brep.half_edges.remove(he_a);
    brep.half_edges.remove(he_b);
    brep.edges.remove(shared_edge);

    // Remove face_b and its loop.
    brep.loops.remove(loop_b);
    brep.faces.remove(face_b);

    Ok(())
}

// -----------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------

/// Return a canonical (smaller, larger) ordering of two face handles.
fn canonical_pair(a: FaceHandle, b: FaceHandle) -> (FaceHandle, FaceHandle) {
    // SlotMap keys implement Ord through their internal representation.
    // Use debug-format comparison as a stable tie-breaker since
    // FaceHandle doesn't implement Ord directly.
    let sa = format!("{a:?}");
    let sb = format!("{b:?}");
    if sa <= sb { (a, b) } else { (b, a) }
}

/// Check that two surfaces are type-compatible for merging.
///
/// Only same-type analytical surfaces are candidates (Plane+Plane,
/// Cylinder+Cylinder, etc.).  NURBS surfaces are not candidates.
fn surface_types_compatible(a: &SurfaceGeometry, b: &SurfaceGeometry) -> bool {
    matches!(
        (a, b),
        (SurfaceGeometry::Plane { .. }, SurfaceGeometry::Plane { .. })
            | (
                SurfaceGeometry::Cylinder { .. },
                SurfaceGeometry::Cylinder { .. }
            )
            | (SurfaceGeometry::Cone { .. }, SurfaceGeometry::Cone { .. })
            | (
                SurfaceGeometry::Sphere { .. },
                SurfaceGeometry::Sphere { .. }
            )
            | (SurfaceGeometry::Torus { .. }, SurfaceGeometry::Torus { .. })
    )
}

/// Additional coplanarity check for planar surfaces.
///
/// Two planes are coplanar if their normals are parallel (already checked
/// by the caller) **and** they have the same signed distance from the
/// origin.  For non-planar surfaces this always returns `true` (the
/// normal-parallelism check is sufficient).
fn surfaces_coplanar(a: &SurfaceGeometry, b: &SurfaceGeometry, linear_tol: f64) -> bool {
    match (a, b) {
        (
            SurfaceGeometry::Plane {
                origin: origin_a,
                normal: normal_a,
                ..
            },
            SurfaceGeometry::Plane {
                origin: origin_b,
                normal: normal_b,
                ..
            },
        ) => {
            // Signed distance from origin for each plane: d = origin . normal.
            let da = origin_a.x * normal_a.x + origin_a.y * normal_a.y + origin_a.z * normal_a.z;
            let db = origin_b.x * normal_b.x + origin_b.y * normal_b.y + origin_b.z * normal_b.z;
            // Account for anti-parallel normals: the signed distances would
            // have opposite sign when normals point in opposite directions.
            // Compare absolute distances when normals are anti-parallel.
            let dot = normal_a.dot(normal_b);
            if dot >= 0.0 {
                (da - db).abs() < linear_tol
            } else {
                (da + db).abs() < linear_tol
            }
        }
        // For non-planar surfaces the normal-parallelism + type check is
        // sufficient; no extra distance check needed.
        _ => true,
    }
}

#[cfg(test)]
#[path = "merge_faces_tests.rs"]
mod tests;
