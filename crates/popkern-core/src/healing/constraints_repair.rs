//! Constrained repair: run healing with constraint preservation.
//!
//! Performs the standard healing pipeline, then checks all constraints.
//! If any constraint is violated by the repair, attempts to project
//! the repair to satisfy constraints. If projection fails, rolls back
//! the specific repair that caused the violation.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::healing::constraints::{
    ConstraintSet, ConstraintViolation, GeometricConstraint, check_constraints,
};
use crate::healing::constraints_helpers::{best_fit_plane, reflect_point};
use crate::healing::options::HealOptions;
use crate::healing::repair::pipeline::HealResult;
use crate::topology::BRep;

// ---------------------------------------------------------------------------
// Constrained repair result
// ---------------------------------------------------------------------------

/// Result of a constrained repair attempt.
///
/// Wraps the underlying [`HealResult`] with constraint-related metadata:
/// how many constraints were checked, satisfied before and after, and
/// whether the repair was modified or rolled back to respect constraints.
#[derive(Debug, Clone)]
pub struct ConstrainedRepairResult {
    /// The underlying repair result.
    pub heal_result: HealResult,
    /// Number of constraints that were checked.
    pub constraints_checked: usize,
    /// Number of constraints satisfied before repair.
    pub satisfied_before: usize,
    /// Number of constraints satisfied after repair.
    pub satisfied_after: usize,
    /// Constraints that would have been violated (repair was modified or rejected).
    pub violations_prevented: Vec<ConstraintViolation>,
    /// Whether the repair was modified to respect constraints.
    pub was_modified: bool,
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Run healing with constraint preservation.
///
/// Performs the standard healing pipeline, then checks all constraints.
/// If any constraint is violated by the repair, attempts to project
/// the repair to satisfy constraints. If projection fails, rolls back
/// the entire repair.
///
/// **Acceptance (from spec):** all specified constraints hold after
/// repair; repairs that would violate constraints are rejected with
/// an explanation.
///
/// # Errors
///
/// Returns `CoreError::Healing` if the underlying pipeline fails and
/// rollback is required.
pub fn heal_with_constraints(
    brep: &mut BRep,
    options: &HealOptions,
    constraints: &ConstraintSet,
) -> Result<ConstrainedRepairResult, CoreError> {
    let total = constraints.len();

    // Fast path: no constraints → just run the pipeline.
    if constraints.is_empty() {
        let heal_result = super::repair::pipeline::run_pipeline(brep, options)?;
        return Ok(ConstrainedRepairResult {
            heal_result,
            constraints_checked: 0,
            satisfied_before: 0,
            satisfied_after: 0,
            violations_prevented: Vec::new(),
            was_modified: false,
        });
    }

    // Step 1: check constraints before repair.
    let violations_before = check_constraints(brep, constraints, &options.tolerance);
    let satisfied_before = total - violations_before.len();

    // Collect indices of already-violated constraints (not our fault).
    let pre_violated: Vec<usize> = violations_before
        .iter()
        .map(|v| v.constraint_index)
        .collect();

    // Step 2: save backup.
    let backup = brep.clone();

    // Step 3: run standard healing.
    let heal_result = super::repair::pipeline::run_pipeline(brep, options)?;

    // Step 4: check constraints after repair.
    let violations_after = check_constraints(brep, constraints, &options.tolerance);
    let satisfied_after_initial = total - violations_after.len();

    // Step 5: find NEW violations (satisfied before but violated after).
    let new_violations: Vec<&ConstraintViolation> = violations_after
        .iter()
        .filter(|v| !pre_violated.contains(&v.constraint_index))
        .collect();

    // Step 6: no new violations → success.
    if new_violations.is_empty() {
        return Ok(ConstrainedRepairResult {
            heal_result,
            constraints_checked: total,
            satisfied_before,
            satisfied_after: satisfied_after_initial,
            violations_prevented: Vec::new(),
            was_modified: false,
        });
    }

    // Step 7: attempt projection for each new violation.
    let violation_indices: Vec<usize> = new_violations.iter().map(|v| v.constraint_index).collect();

    for &idx in &violation_indices {
        if let Some(constraint) = constraints.iter().nth(idx) {
            attempt_projection(brep, constraint, options);
        }
    }

    // Step 8: re-check constraints after projection.
    let violations_post_projection = check_constraints(brep, constraints, &options.tolerance);
    let still_new: Vec<ConstraintViolation> = violations_post_projection
        .into_iter()
        .filter(|v| !pre_violated.contains(&v.constraint_index))
        .collect();

    if still_new.is_empty() {
        // Projection succeeded.
        let satisfied_after =
            total - check_constraints(brep, constraints, &options.tolerance).len();
        return Ok(ConstrainedRepairResult {
            heal_result,
            constraints_checked: total,
            satisfied_before,
            satisfied_after,
            violations_prevented: new_violations.into_iter().cloned().collect(),
            was_modified: true,
        });
    }

    // Step 9: projection failed — restore backup.
    let violations_prevented: Vec<ConstraintViolation> =
        new_violations.into_iter().cloned().collect();
    *brep = backup;

    Ok(ConstrainedRepairResult {
        heal_result,
        constraints_checked: total,
        satisfied_before,
        satisfied_after: satisfied_before, // Same as before — we rolled back.
        violations_prevented,
        was_modified: false,
    })
}

// ---------------------------------------------------------------------------
// Projection helpers
// ---------------------------------------------------------------------------

/// Attempt to project a repaired BRep back to satisfy a single constraint.
///
/// Modifies the BRep in place. For constraints that cannot be projected
/// (e.g., complex tangency or symmetry), this is a best-effort operation.
fn attempt_projection(brep: &mut BRep, constraint: &GeometricConstraint, options: &HealOptions) {
    match constraint {
        GeometricConstraint::FixedPoint {
            vertex,
            position,
            max_displacement: _,
        } => {
            project_fixed_point(brep, *vertex, position);
        }
        GeometricConstraint::Tangency {
            face_a: _,
            face_b: _,
            max_angle: _,
        } => {
            // Tangency projection requires surface adjustment which is
            // beyond simple vertex moves. We rely on the G1 enforcement
            // pass in the pipeline. If it still fails, rollback handles it.
        }
        GeometricConstraint::Symmetry {
            faces_a,
            faces_b,
            plane,
            max_deviation: _,
        } => {
            project_symmetry(brep, faces_a, faces_b, &plane.0, &plane.1);
        }
        GeometricConstraint::Coplanarity {
            faces,
            max_deviation: _,
        } => {
            project_coplanarity(brep, faces);
        }
        GeometricConstraint::Distance {
            entity_a,
            entity_b,
            target_distance,
            tolerance: _,
        } => {
            project_distance(brep, entity_a, entity_b, *target_distance, options);
        }
    }
}

/// Move a vertex back to its fixed position.
fn project_fixed_point(brep: &mut BRep, vertex: crate::topology::VertexHandle, position: &Point3) {
    if let Some(v) = brep.vertices.get_mut(vertex) {
        v.position = *position;
    }
}

/// Restore symmetry by averaging face vertex positions with their reflections.
fn project_symmetry(
    brep: &mut BRep,
    faces_a: &[crate::topology::FaceHandle],
    faces_b: &[crate::topology::FaceHandle],
    plane_point: &Point3,
    plane_normal: &crate::geometry::Vec3,
) {
    let pair_count = faces_a.len().min(faces_b.len());

    for i in 0..pair_count {
        // Collect vertex handles and positions for both sides.
        let verts_a = collect_face_vertices(brep, faces_a[i]);
        let verts_b = collect_face_vertices(brep, faces_b[i]);

        // For each vertex on side A, compute its reflection and average
        // with the corresponding vertex on side B (by index).
        let count = verts_a.len().min(verts_b.len());
        for j in 0..count {
            let (vh_a, pos_a) = verts_a[j];
            let (vh_b, pos_b) = verts_b[j];

            let reflected_a = reflect_point(&pos_a, plane_point, plane_normal);
            let reflected_b = reflect_point(&pos_b, plane_point, plane_normal);

            // Average: move A toward the reflection of B, and B toward
            // the reflection of A.
            let new_a = Point3::new(
                (pos_a.x + reflected_b.x) * 0.5,
                (pos_a.y + reflected_b.y) * 0.5,
                (pos_a.z + reflected_b.z) * 0.5,
            );
            let new_b = Point3::new(
                (pos_b.x + reflected_a.x) * 0.5,
                (pos_b.y + reflected_a.y) * 0.5,
                (pos_b.z + reflected_a.z) * 0.5,
            );

            if let Some(va) = brep.vertices.get_mut(vh_a) {
                va.position = new_a;
            }
            if let Some(vb) = brep.vertices.get_mut(vh_b) {
                vb.position = new_b;
            }
        }
    }
}

/// Project face vertices back onto the best-fit plane.
fn project_coplanarity(brep: &mut BRep, faces: &[crate::topology::FaceHandle]) {
    // Collect all vertex positions.
    let mut all_points = Vec::new();
    let mut all_verts = Vec::new();
    for &fh in faces {
        let verts = collect_face_vertices(brep, fh);
        for (vh, pos) in &verts {
            all_points.push(*pos);
            all_verts.push(*vh);
        }
    }

    if let Some((plane_pt, plane_n)) = best_fit_plane(&all_points) {
        // Project each vertex onto the plane.
        for (i, vh) in all_verts.iter().enumerate() {
            let diff = all_points[i] - plane_pt;
            let dist = diff.dot(&plane_n);
            let projected = all_points[i] - plane_n * dist;
            if let Some(v) = brep.vertices.get_mut(*vh) {
                v.position = projected;
            }
        }
    }
}

/// Attempt to adjust vertex positions to satisfy a distance constraint.
fn project_distance(
    brep: &mut BRep,
    entity_a: &crate::healing::diagnostics::EntityRef,
    entity_b: &crate::healing::diagnostics::EntityRef,
    target_distance: f64,
    _options: &HealOptions,
) {
    // Only handle vertex-vertex distance for projection.
    use crate::healing::diagnostics::EntityRef;
    if let (EntityRef::Vertex(va), EntityRef::Vertex(vb)) = (entity_a, entity_b) {
        let pos_a = match brep.vertices.get(*va) {
            Some(v) => v.position,
            None => return,
        };
        let pos_b = match brep.vertices.get(*vb) {
            Some(v) => v.position,
            None => return,
        };
        let actual = pos_a.distance_to(&pos_b);
        if actual < f64::EPSILON {
            return;
        }
        let scale = target_distance / actual;
        let midpoint = Point3::new(
            (pos_a.x + pos_b.x) * 0.5,
            (pos_a.y + pos_b.y) * 0.5,
            (pos_a.z + pos_b.z) * 0.5,
        );
        let half_dir = crate::geometry::Vec3::new(
            (pos_b.x - pos_a.x) * 0.5 * scale,
            (pos_b.y - pos_a.y) * 0.5 * scale,
            (pos_b.z - pos_a.z) * 0.5 * scale,
        );
        if let Some(va_data) = brep.vertices.get_mut(*va) {
            va_data.position = midpoint - half_dir;
        }
        if let Some(vb_data) = brep.vertices.get_mut(*vb) {
            vb_data.position = midpoint + half_dir;
        }
    }
}

/// Collect vertex handles and positions from a face's outer loop.
fn collect_face_vertices(
    brep: &BRep,
    face: crate::topology::FaceHandle,
) -> Vec<(crate::topology::VertexHandle, Point3)> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let Some(loop_) = brep.loops.get(face_data.outer_loop) else {
        return Vec::new();
    };
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut result = Vec::new();
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let vertex = &brep.vertices[he.origin];
        result.push((he.origin, vertex.position));
        current = he.next;
        if current == start_he {
            break;
        }
    }
    result
}

#[cfg(test)]
#[path = "constraints_repair_tests.rs"]
mod constraints_repair_tests;
