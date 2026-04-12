//! Constraint-preserving healing: respect design-intent constraints
//! during geometry repair.
//!
//! Provides a 3D geometric constraint system for BRep entities.
//! Repairs that would violate specified constraints are either
//! projected to satisfy them or rejected with an explanation.
//!
//! This is a standalone 3D constraint system within the healing module,
//! distinct from the 2D sketch constraints in `popkern-constraint`.

use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};
use crate::healing::constraints_helpers::{clamp_param, compute_centroid};
use crate::healing::diagnostics::EntityRef;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle, VertexHandle};

// Re-export helpers so tests (via `super::*`) and sibling modules can use them.
pub(crate) use crate::healing::constraints_helpers::{
    best_fit_plane, entity_representative_point, reflect_point,
};

// ---------------------------------------------------------------------------
// Constraint types
// ---------------------------------------------------------------------------

/// A 3D geometric constraint between BRep entities.
///
/// Constraints express design intent that must be preserved during
/// healing. Each variant carries the entities involved, the constraint
/// parameters, and a tolerance for checking satisfaction.
#[derive(Debug, Clone)]
pub enum GeometricConstraint {
    /// Two surfaces must be tangent (G1) at their shared boundary.
    Tangency {
        /// First face of the tangency pair.
        face_a: FaceHandle,
        /// Second face of the tangency pair.
        face_b: FaceHandle,
        /// Maximum allowed normal angle deviation (radians).
        max_angle: f64,
    },

    /// Geometry must be symmetric about a plane.
    Symmetry {
        /// Faces on one side of the symmetry plane.
        faces_a: Vec<FaceHandle>,
        /// Corresponding faces on the other side.
        faces_b: Vec<FaceHandle>,
        /// The symmetry plane: (point on plane, unit normal).
        plane: (Point3, Vec3),
        /// Maximum allowed asymmetry (distance).
        max_deviation: f64,
    },

    /// A set of faces must be coplanar.
    Coplanarity {
        /// Faces that must lie in the same plane.
        faces: Vec<FaceHandle>,
        /// Maximum allowed deviation from the best-fit plane.
        max_deviation: f64,
    },

    /// A vertex must remain at a fixed position.
    FixedPoint {
        /// The constrained vertex.
        vertex: VertexHandle,
        /// The required position.
        position: Point3,
        /// Maximum allowed displacement from the fixed position.
        max_displacement: f64,
    },

    /// Two entities must maintain a specific distance.
    Distance {
        /// First entity in the distance pair.
        entity_a: EntityRef,
        /// Second entity in the distance pair.
        entity_b: EntityRef,
        /// Target distance between the entities.
        target_distance: f64,
        /// Maximum allowed deviation from target distance.
        tolerance: f64,
    },
}

// ---------------------------------------------------------------------------
// Constraint set
// ---------------------------------------------------------------------------

/// A set of geometric constraints to be preserved during repair.
///
/// Collects constraints and provides iteration for checking. The set
/// is ordered: constraints are checked in insertion order.
#[derive(Debug, Clone, Default)]
pub struct ConstraintSet {
    constraints: Vec<GeometricConstraint>,
}

impl ConstraintSet {
    /// Create an empty constraint set.
    pub fn new() -> Self {
        Self {
            constraints: Vec::new(),
        }
    }

    /// Add a constraint to the set.
    pub fn add(&mut self, constraint: GeometricConstraint) {
        self.constraints.push(constraint);
    }

    /// Number of constraints in the set.
    pub fn len(&self) -> usize {
        self.constraints.len()
    }

    /// Returns true if the set contains no constraints.
    pub fn is_empty(&self) -> bool {
        self.constraints.is_empty()
    }

    /// Iterate over all constraints.
    pub fn iter(&self) -> impl Iterator<Item = &GeometricConstraint> {
        self.constraints.iter()
    }
}

// ---------------------------------------------------------------------------
// Violation reporting
// ---------------------------------------------------------------------------

/// Result of checking a single constraint.
///
/// Captures which constraint was violated, a human-readable description
/// of the violation, and the magnitude of the violation for ranking.
#[derive(Debug, Clone)]
pub struct ConstraintViolation {
    /// Index of the violated constraint in the [`ConstraintSet`].
    pub constraint_index: usize,
    /// Description of the violation.
    pub description: String,
    /// How much the constraint is violated (distance or angle).
    pub violation_magnitude: f64,
}

// ---------------------------------------------------------------------------
// Constraint checking
// ---------------------------------------------------------------------------

/// Check all constraints against the current BRep state.
///
/// Returns a list of violated constraints. An empty list means all
/// constraints are satisfied.
///
/// For each constraint type:
/// - **Tangency**: samples normals on both faces near the midpoint
///   and checks that the angle between them is within `max_angle`.
/// - **Symmetry**: samples points on each face pair, reflects across
///   the symmetry plane, and measures maximum deviation.
/// - **Coplanarity**: samples points on all faces, fits a best-fit
///   plane, and checks maximum deviation.
/// - **FixedPoint**: checks vertex position against the fixed position.
/// - **Distance**: measures distance between entity representative
///   points and checks against the target.
pub fn check_constraints(
    brep: &BRep,
    constraints: &ConstraintSet,
    tolerance: &Tolerance,
) -> Vec<ConstraintViolation> {
    let mut violations = Vec::new();
    for (idx, constraint) in constraints.iter().enumerate() {
        if let Some(v) = check_one(brep, idx, constraint, tolerance) {
            violations.push(v);
        }
    }
    violations
}

/// Check a single constraint, returning a violation if unsatisfied.
fn check_one(
    brep: &BRep,
    index: usize,
    constraint: &GeometricConstraint,
    _tolerance: &Tolerance,
) -> Option<ConstraintViolation> {
    match constraint {
        GeometricConstraint::Tangency {
            face_a,
            face_b,
            max_angle,
        } => check_tangency(brep, index, *face_a, *face_b, *max_angle),

        GeometricConstraint::Symmetry {
            faces_a,
            faces_b,
            plane,
            max_deviation,
        } => check_symmetry(brep, index, faces_a, faces_b, plane, *max_deviation),

        GeometricConstraint::Coplanarity {
            faces,
            max_deviation,
        } => check_coplanarity(brep, index, faces, *max_deviation),

        GeometricConstraint::FixedPoint {
            vertex,
            position,
            max_displacement,
        } => check_fixed_point(brep, index, *vertex, position, *max_displacement),

        GeometricConstraint::Distance {
            entity_a,
            entity_b,
            target_distance,
            tolerance,
        } => check_distance(
            brep,
            index,
            entity_a,
            entity_b,
            *target_distance,
            *tolerance,
        ),
    }
}

/// Check tangency between two faces by comparing surface normals.
fn check_tangency(
    brep: &BRep,
    index: usize,
    face_a: FaceHandle,
    face_b: FaceHandle,
    max_angle: f64,
) -> Option<ConstraintViolation> {
    let surf_a = brep.faces.get(face_a)?.surface.as_ref()?;
    let surf_b = brep.faces.get(face_b)?.surface.as_ref()?;

    // Sample at the midpoint of each face's parameter range.
    let ((ua0, ua1), (va0, va1)) = surf_a.geometry.parameter_range();
    let ua_mid = clamp_param(ua0, ua1);
    let va_mid = clamp_param(va0, va1);

    let ((ub0, ub1), (vb0, vb1)) = surf_b.geometry.parameter_range();
    let ub_mid = clamp_param(ub0, ub1);
    let vb_mid = clamp_param(vb0, vb1);

    let n_a = surf_a.geometry.normal_at(ua_mid, va_mid);
    let n_b = surf_b.geometry.normal_at(ub_mid, vb_mid);

    // Tangency means normals should be parallel (angle ~ 0 or ~ pi).
    let angle = n_a.angle_to(&n_b);
    let tangency_angle = angle.min(std::f64::consts::PI - angle);

    if tangency_angle > max_angle {
        Some(ConstraintViolation {
            constraint_index: index,
            description: format!(
                "Tangency violated: normal angle {tangency_angle:.6} rad > max {max_angle:.6} rad"
            ),
            violation_magnitude: tangency_angle - max_angle,
        })
    } else {
        None
    }
}

/// Check symmetry by reflecting face sample points across the plane.
fn check_symmetry(
    brep: &BRep,
    index: usize,
    faces_a: &[FaceHandle],
    faces_b: &[FaceHandle],
    plane: &(Point3, Vec3),
    max_deviation: f64,
) -> Option<ConstraintViolation> {
    let mut max_dist: f64 = 0.0;

    let pair_count = faces_a.len().min(faces_b.len());
    for i in 0..pair_count {
        let pts_a = brep.face_outer_loop_positions(faces_a[i]);
        let pts_b = brep.face_outer_loop_positions(faces_b[i]);

        // Reflect points from side A and measure distance to closest on side B.
        for pt in &pts_a {
            let reflected = reflect_point(pt, &plane.0, &plane.1);
            let min_d = pts_b
                .iter()
                .map(|pb| reflected.distance_to(pb))
                .fold(f64::INFINITY, f64::min);
            max_dist = max_dist.max(min_d);
        }
    }

    if max_dist > max_deviation {
        Some(ConstraintViolation {
            constraint_index: index,
            description: format!(
                "Symmetry violated: max deviation {max_dist:.6} > limit {max_deviation:.6}"
            ),
            violation_magnitude: max_dist - max_deviation,
        })
    } else {
        None
    }
}

/// Check coplanarity by fitting a best-fit plane and measuring deviation.
fn check_coplanarity(
    brep: &BRep,
    index: usize,
    faces: &[FaceHandle],
    max_deviation: f64,
) -> Option<ConstraintViolation> {
    let mut all_points = Vec::new();
    for &fh in faces {
        all_points.extend(brep.face_outer_loop_positions(fh));
    }

    if all_points.len() < 3 {
        return None; // Not enough points to define a plane.
    }

    let (_plane_point, plane_normal) = best_fit_plane(&all_points)?;

    // Measure max deviation from the best-fit plane.
    let centroid = compute_centroid(&all_points);
    let mut max_dev: f64 = 0.0;
    for pt in &all_points {
        let diff = *pt - centroid;
        let dist = diff.dot(&plane_normal).abs();
        max_dev = max_dev.max(dist);
    }

    if max_dev > max_deviation {
        Some(ConstraintViolation {
            constraint_index: index,
            description: format!(
                "Coplanarity violated: max deviation {max_dev:.6} > limit {max_deviation:.6}"
            ),
            violation_magnitude: max_dev - max_deviation,
        })
    } else {
        None
    }
}

/// Check that a vertex is at its fixed position.
fn check_fixed_point(
    brep: &BRep,
    index: usize,
    vertex: VertexHandle,
    position: &Point3,
    max_displacement: f64,
) -> Option<ConstraintViolation> {
    let v = brep.vertices.get(vertex)?;
    let dist = v.position.distance_to(position);
    if dist > max_displacement {
        Some(ConstraintViolation {
            constraint_index: index,
            description: format!(
                "FixedPoint violated: displacement {dist:.6} > limit {max_displacement:.6}"
            ),
            violation_magnitude: dist - max_displacement,
        })
    } else {
        None
    }
}

/// Check distance constraint between two entities.
fn check_distance(
    brep: &BRep,
    index: usize,
    entity_a: &EntityRef,
    entity_b: &EntityRef,
    target_distance: f64,
    tolerance: f64,
) -> Option<ConstraintViolation> {
    let pa = entity_representative_point(brep, entity_a)?;
    let pb = entity_representative_point(brep, entity_b)?;
    let actual = pa.distance_to(&pb);
    let error = (actual - target_distance).abs();
    if error > tolerance {
        Some(ConstraintViolation {
            constraint_index: index,
            description: format!(
                "Distance violated: actual {actual:.6}, target {target_distance:.6}, \
                 error {error:.6} > tolerance {tolerance:.6}"
            ),
            violation_magnitude: error - tolerance,
        })
    } else {
        None
    }
}

#[cfg(test)]
#[path = "constraints_tests.rs"]
mod constraints_tests;
