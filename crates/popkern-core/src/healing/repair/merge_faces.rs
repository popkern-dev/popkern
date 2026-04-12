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
///
/// Currently a detection-only stub: returns the detected pairs count
/// but does not modify topology (face merging requires careful Euler
/// operator sequences that are deferred to a future iteration).
pub fn merge_coplanar_faces(brep: &mut BRep, tolerance: &Tolerance) -> FaceMergeResult {
    let pairs = detect_coplanar_pairs(brep, tolerance);
    let pairs_detected = pairs.len();

    FaceMergeResult {
        pairs_detected,
        pairs_merged: 0,
    }
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
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn unit_cube_has_no_coplanar_pairs() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let pairs = detect_coplanar_pairs(&brep, &tol);
        assert!(
            pairs.is_empty(),
            "all cube faces are perpendicular, expected 0 coplanar pairs, got {}",
            pairs.len()
        );
    }

    #[test]
    fn merge_coplanar_faces_on_cube_detects_zero() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let result = merge_coplanar_faces(&mut brep, &tol);
        assert_eq!(result.pairs_detected, 0);
        assert_eq!(result.pairs_merged, 0);
    }

    #[test]
    fn empty_brep_has_no_coplanar_pairs() {
        let brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let pairs = detect_coplanar_pairs(&brep, &tol);
        assert!(pairs.is_empty());
    }

    #[test]
    fn canonical_pair_is_stable() {
        // Verify the deduplication helper is commutative.
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let faces: Vec<FaceHandle> = brep.faces.keys().collect();
        if faces.len() >= 2 {
            let a = faces[0];
            let b = faces[1];
            assert_eq!(canonical_pair(a, b), canonical_pair(b, a));
        }
    }

    #[test]
    fn surface_types_compatible_same() {
        use crate::geometry::{Point3, Vec3};
        let plane_a = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let plane_b = SurfaceGeometry::Plane {
            origin: Point3::new(1.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        assert!(surface_types_compatible(&plane_a, &plane_b));
    }

    #[test]
    fn surface_types_incompatible_cross() {
        use crate::geometry::{Point3, Vec3};
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        assert!(!surface_types_compatible(&plane, &cyl));
    }
}
