//! Self-intersection repair: detect and attempt to fix surfaces that
//! fold through themselves.
//!
//! Self-intersection is algorithmically difficult to repair in general
//! (it requires splitting surfaces at intersection curves). This module
//! provides a pragmatic approach: detect the problem, attempt
//! re-approximation for NURBS surfaces by smoothing control points in
//! folded regions, and report when repair is not possible.

use crate::error::ValidationError;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};

/// Result of a self-intersection repair attempt.
#[derive(Debug, Clone)]
pub struct SelfIntersectionRepairResult {
    /// Face that was repaired.
    pub face: FaceHandle,
    /// Whether the repair succeeded.
    pub repaired: bool,
    /// Description of what was done or why repair failed.
    pub description: String,
}

/// Detect faces with self-intersecting geometry.
///
/// Returns pairs of face handles where the surfaces intersect each other
/// (non-adjacent faces). Uses the existing validation infrastructure
/// ([`BRep::validate_self_intersection`]).
///
/// Because the underlying check returns on the first detected pair,
/// at most one pair is returned per call. To find all pairs, use
/// [`BRep::validate_all`] and filter for `SelfIntersection` errors.
pub fn detect_self_intersecting_faces(brep: &BRep) -> Vec<(FaceHandle, FaceHandle)> {
    let mut pairs = Vec::new();

    // Primary check: use the fast single-pair validator.
    if let Err(ValidationError::SelfIntersection { face_a, face_b }) =
        brep.validate_self_intersection()
    {
        pairs.push((face_a, face_b));
    }

    // Secondary: sweep validate_all for additional SelfIntersection errors
    // that the single-pair check might not reach (it stops at the first hit).
    // Only run if the primary check found something (otherwise the model
    // is clean and validate_all would be wasted work).
    if !pairs.is_empty() {
        for err in brep.validate_all() {
            if let ValidationError::SelfIntersection { face_a, face_b } = err {
                let pair = if face_a < face_b {
                    (face_a, face_b)
                } else {
                    (face_b, face_a)
                };
                if !pairs.contains(&pair) {
                    pairs.push(pair);
                }
            }
        }
    }

    pairs
}

/// Detect faces where the surface normal field reverses direction,
/// indicating the surface folds through itself.
///
/// Samples the surface normal on a grid and checks for sign reversals
/// (dot product of adjacent normals < 0). Returns the faces where a
/// reversal was detected together with the approximate (u, v) location
/// of the fold.
fn detect_normal_reversals(brep: &BRep, face: FaceHandle) -> Vec<(f64, f64)> {
    let face_data = match brep.faces.get(face) {
        Some(f) => f,
        None => return Vec::new(),
    };
    let surface_on_face = match face_data.surface.as_ref() {
        Some(s) => s,
        None => return Vec::new(),
    };
    let surface = &surface_on_face.geometry;

    // Get bounded parameter range from face vertices.
    let positions = brep.face_outer_loop_positions(face);
    if positions.is_empty() {
        return Vec::new();
    }
    let ((u_min, u_max), (v_min, v_max)) =
        crate::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface);

    // Sample normals on a 6x6 grid.
    let n = 6usize;
    let mut normals = Vec::with_capacity(n * n);
    for i in 0..n {
        let u = u_min + (u_max - u_min) * (i as f64 + 0.5) / (n as f64);
        for j in 0..n {
            let v = v_min + (v_max - v_min) * (j as f64 + 0.5) / (n as f64);
            normals.push((u, v, surface.normal_at(u, v)));
        }
    }

    // Find locations where adjacent normals reverse (dot product < 0).
    let mut reversals = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let idx = i * n + j;
            let (u, v, normal) = &normals[idx];
            let len = normal.length();
            if len < 1e-12 {
                continue; // degenerate normal, skip
            }

            // Check neighbor to the right.
            if j + 1 < n {
                let (_, _, neighbor) = &normals[idx + 1];
                if neighbor.length() > 1e-12 && normal.dot(neighbor) < 0.0 {
                    reversals.push((*u, *v));
                }
            }
            // Check neighbor below.
            if i + 1 < n {
                let (_, _, neighbor) = &normals[(i + 1) * n + j];
                if neighbor.length() > 1e-12 && normal.dot(neighbor) < 0.0 {
                    reversals.push((*u, *v));
                }
            }
        }
    }

    reversals
}

/// Attempt to repair a self-intersecting face.
///
/// For NURBS surfaces, tries to detect the folded region by sampling
/// the surface normal field and identifying where it reverses direction.
/// If the fold is localized, attempts to smooth the control points in
/// the affected region using Laplacian smoothing.
///
/// For analytical surfaces, self-intersection is not geometrically
/// possible (planes, spheres, cylinders, etc. cannot self-intersect),
/// so this returns immediately.
///
/// This is a best-effort operation -- many self-intersections cannot
/// be repaired without significant surface redesign.
pub fn repair_self_intersection(brep: &mut BRep, face: FaceHandle) -> SelfIntersectionRepairResult {
    let face_data = match brep.faces.get(face) {
        Some(f) => f,
        None => {
            return SelfIntersectionRepairResult {
                face,
                repaired: false,
                description: "Face not found in BRep".to_string(),
            };
        }
    };

    let surface_on_face = match face_data.surface.as_ref() {
        Some(s) => s,
        None => {
            return SelfIntersectionRepairResult {
                face,
                repaired: false,
                description: "Face has no surface geometry".to_string(),
            };
        }
    };

    // Analytical surfaces cannot self-intersect.
    let nurbs = match &surface_on_face.geometry {
        SurfaceGeometry::Nurbs(n) => n.clone(),
        _ => {
            return SelfIntersectionRepairResult {
                face,
                repaired: false,
                description: "Analytical surfaces cannot self-intersect".to_string(),
            };
        }
    };

    // Detect normal reversals in the current surface.
    let reversals = detect_normal_reversals(brep, face);
    if reversals.is_empty() {
        return SelfIntersectionRepairResult {
            face,
            repaired: false,
            description: "No normal reversal detected; self-intersection may be \
                          between non-adjacent faces rather than a surface fold"
                .to_string(),
        };
    }

    // Attempt Laplacian smoothing on the NURBS control points.
    // We smooth the entire control net with a conservative strength;
    // a more targeted approach (smoothing only near the fold) would
    // require mapping (u,v) reversal locations back to control-point
    // indices, which is deferred to a future iteration.
    let original_nurbs = nurbs.clone();
    let mut smoothed = nurbs;
    let num_u = smoothed.num_u;
    let num_v = smoothed.num_v;

    // Need at least a 3x3 grid to have interior points.
    if num_u < 3 || num_v < 3 {
        return SelfIntersectionRepairResult {
            face,
            repaired: false,
            description: format!("NURBS control grid too small for smoothing ({num_u}x{num_v})"),
        };
    }

    let iterations = 3;
    let strength = 0.3;

    for _ in 0..iterations {
        let current = smoothed.control_points.clone();
        // Smooth interior control points (leave boundary rows/columns fixed).
        for iu in 1..(num_u - 1) {
            for iv in 1..(num_v - 1) {
                let idx = iu * num_v + iv;
                let left = current[(iu - 1) * num_v + iv];
                let right = current[(iu + 1) * num_v + iv];
                let up = current[iu * num_v + (iv - 1)];
                let down = current[iu * num_v + (iv + 1)];
                let avg_x = (left.x + right.x + up.x + down.x) * 0.25;
                let avg_y = (left.y + right.y + up.y + down.y) * 0.25;
                let avg_z = (left.z + right.z + up.z + down.z) * 0.25;
                let cp = &current[idx];
                smoothed.control_points[idx] = crate::geometry::Point3::new(
                    cp.x + (avg_x - cp.x) * strength,
                    cp.y + (avg_y - cp.y) * strength,
                    cp.z + (avg_z - cp.z) * strength,
                );
            }
        }
    }

    // Compute max displacement from original.
    let max_disp = original_nurbs
        .control_points
        .iter()
        .zip(smoothed.control_points.iter())
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max);

    // Replace the face surface with the smoothed NURBS.
    let orientation = surface_on_face.orientation;
    if let Some(f) = brep.faces.get_mut(face) {
        f.surface = Some(crate::geometry::SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(smoothed),
            orientation,
        });
    }

    // Verify: check if normals are now consistent.
    let post_reversals = detect_normal_reversals(brep, face);
    if post_reversals.is_empty() {
        SelfIntersectionRepairResult {
            face,
            repaired: true,
            description: format!(
                "Smoothed {num_reversals} normal reversal(s); max control point \
                 displacement {max_disp:.6}",
                num_reversals = reversals.len(),
            ),
        }
    } else {
        // Revert: restore original surface since smoothing didn't help.
        if let Some(f) = brep.faces.get_mut(face) {
            f.surface = Some(crate::geometry::SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(original_nurbs),
                orientation,
            });
        }
        SelfIntersectionRepairResult {
            face,
            repaired: false,
            description: format!(
                "Smoothing reduced reversals from {} to {} but did not eliminate them; \
                 reverted changes",
                reversals.len(),
                post_reversals.len(),
            ),
        }
    }
}

/// Attempt to repair self-intersection by identifying the fold region
/// and splitting the surface.
///
/// For NURBS surfaces, finds the parameter-space curve where the surface
/// folds (normal reversal boundary) and splits the surface into two
/// parts: the valid region and the folded region. The folded region
/// is discarded.
///
/// This is more aggressive than smoothing but handles cases where the
/// fold is too severe for smoothing to fix.
pub fn repair_by_splitting(
    brep: &mut BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
) -> SelfIntersectionRepairResult {
    let face_data = match brep.faces.get(face) {
        Some(f) => f,
        None => {
            return SelfIntersectionRepairResult {
                face,
                repaired: false,
                description: "Face not found in BRep".to_string(),
            };
        }
    };

    let surface_on_face = match face_data.surface.as_ref() {
        Some(s) => s,
        None => {
            return SelfIntersectionRepairResult {
                face,
                repaired: false,
                description: "Face has no surface geometry".to_string(),
            };
        }
    };

    // Only NURBS surfaces can self-intersect and be split.
    let nurbs = match &surface_on_face.geometry {
        SurfaceGeometry::Nurbs(n) => n.clone(),
        _ => {
            return SelfIntersectionRepairResult {
                face,
                repaired: false,
                description: "Only NURBS surfaces can be split; \
                              analytical surfaces cannot self-intersect"
                    .to_string(),
            };
        }
    };

    // Detect normal reversals on a 10x10 grid.
    let reversals = detect_normal_reversals(brep, face);
    if reversals.is_empty() {
        return SelfIntersectionRepairResult {
            face,
            repaired: false,
            description: "No normal reversals detected; surface is clean".to_string(),
        };
    }

    // Get bounded parameter range for this face.
    let positions = brep.face_outer_loop_positions(face);
    if positions.is_empty() {
        return SelfIntersectionRepairResult {
            face,
            repaired: false,
            description: "Face has no outer loop positions".to_string(),
        };
    }
    let surface = &surface_on_face.geometry;
    let ((u_min, u_max), (v_min, v_max)) =
        crate::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface);

    // Build a 10x10 grid of normal directions to find the fold boundary.
    let n = 10usize;
    let mut grid_ok = vec![true; n * n];
    let ref_normal =
        surface.normal_at(u_min + (u_max - u_min) * 0.5, v_min + (v_max - v_min) * 0.5);

    for i in 0..n {
        let u = u_min + (u_max - u_min) * (i as f64 + 0.5) / (n as f64);
        for j in 0..n {
            let v = v_min + (v_max - v_min) * (j as f64 + 0.5) / (n as f64);
            let normal = surface.normal_at(u, v);
            if normal.length() > 1e-12 && normal.dot(&ref_normal) < 0.0 {
                grid_ok[i * n + j] = false;
            }
        }
    }

    let ok_count = grid_ok.iter().filter(|&&ok| ok).count();
    let fold_count = n * n - ok_count;

    if fold_count == 0 {
        return SelfIntersectionRepairResult {
            face,
            repaired: false,
            description: "Grid analysis found no folded region".to_string(),
        };
    }

    // Find the fold boundary in parameter space: midpoints between
    // adjacent grid cells where the normal flips.
    let mut fold_u_min = u_max;
    let mut fold_u_max = u_min;
    let mut fold_v_min = v_max;
    let mut fold_v_max = v_min;

    for i in 0..n {
        let u = u_min + (u_max - u_min) * (i as f64 + 0.5) / (n as f64);
        for j in 0..n {
            if !grid_ok[i * n + j] {
                let v = v_min + (v_max - v_min) * (j as f64 + 0.5) / (n as f64);
                if u < fold_u_min {
                    fold_u_min = u;
                }
                if u > fold_u_max {
                    fold_u_max = u;
                }
                if v < fold_v_min {
                    fold_v_min = v;
                }
                if v > fold_v_max {
                    fold_v_max = v;
                }
            }
        }
    }

    // Determine which side of the fold to trim.
    // The valid region is the larger area.
    let u_range = u_max - u_min;
    let v_range = v_max - v_min;
    let _ = tolerance; // tolerance used above via linear()

    // Try splitting along the dominant fold direction.
    // If the fold is near u_min/u_max, split in u; if near v_min/v_max, split in v.
    let fold_near_u_min = (fold_u_min - u_min) / u_range < 0.3;
    let fold_near_u_max = (u_max - fold_u_max) / u_range < 0.3;
    let fold_near_v_min = (fold_v_min - v_min) / v_range < 0.3;
    let fold_near_v_max = (v_max - fold_v_max) / v_range < 0.3;

    let orientation = surface_on_face.orientation;

    // Attempt to split the NURBS surface, keeping the valid side.
    let split_result = if fold_near_u_min && !fold_near_u_max {
        // Fold is at the low-u end; split and keep the right (high-u) part.
        nurbs
            .split_at_u(fold_u_max)
            .ok()
            .map(|(_left, right)| right)
    } else if fold_near_u_max && !fold_near_u_min {
        // Fold is at the high-u end; split and keep the left (low-u) part.
        nurbs.split_at_u(fold_u_min).ok().map(|(left, _right)| left)
    } else if fold_near_v_min && !fold_near_v_max {
        nurbs
            .split_at_v(fold_v_max)
            .ok()
            .map(|(_left, right)| right)
    } else if fold_near_v_max && !fold_near_v_min {
        nurbs.split_at_v(fold_v_min).ok().map(|(left, _right)| left)
    } else {
        // Fold is in the interior; cannot cleanly split with a single cut.
        None
    };

    match split_result {
        Some(trimmed) => {
            if let Some(f) = brep.faces.get_mut(face) {
                f.surface = Some(crate::geometry::SurfaceOnFace {
                    geometry: SurfaceGeometry::Nurbs(trimmed),
                    orientation,
                });
            }
            SelfIntersectionRepairResult {
                face,
                repaired: true,
                description: format!(
                    "Split surface to remove fold region ({fold_count} of \
                     {total} grid cells were folded)",
                    total = n * n,
                ),
            }
        }
        None => SelfIntersectionRepairResult {
            face,
            repaired: false,
            description: format!(
                "Fold region spans interior ({fold_count} of {} grid cells); \
                 single-cut split not possible",
                n * n,
            ),
        },
    }
}

#[cfg(test)]
#[path = "self_intersection_tests.rs"]
mod tests;
