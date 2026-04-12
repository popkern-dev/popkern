//! Topological defect detection (Section 5.1).
//!
//! Wraps the existing BRep validation methods and converts their
//! `ValidationError` results into structured `DefectRecord`s for
//! the healing report.

use crate::error::ValidationError;
use crate::healing::diagnostics::{
    DefectCategory, DefectRecord, EntityRef, HealingReport, Severity,
};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Run all topological defect checks.
pub fn check_all(brep: &BRep, tolerance: &Tolerance, report: &mut HealingReport) {
    check_boundary_closure(brep, tolerance, report);
    check_half_edge_pairing(brep, report);
    check_half_edge_cycles(brep, report);
    check_vertex_consistency(brep, report);
    check_manifold_edges(brep, report);
    check_face_orientation(brep, report);
    check_shell_closure(brep, report);
}

/// Check that every loop is geometrically closed within tolerance.
///
/// Calls `brep.validate_loop_closure()` and converts any
/// `ValidationError::OpenLoop` into a defect record.
pub fn check_boundary_closure(brep: &BRep, _tolerance: &Tolerance, report: &mut HealingReport) {
    if let Err(e) = brep.validate_loop_closure() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that every half-edge has a valid twin and twin-of-twin consistency.
///
/// Calls `brep.validate_half_edge_pairing()` and converts any
/// `ValidationError::BrokenTwin` into a defect record.
pub fn check_half_edge_pairing(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_half_edge_pairing() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that following `next` from any half-edge returns to the start.
///
/// Calls `brep.validate_half_edge_cycles()` and converts any
/// `ValidationError::BrokenCycle` into a defect record.
pub fn check_half_edge_cycles(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_half_edge_cycles() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that each vertex's outgoing half-edge originates from it.
///
/// Calls `brep.validate_vertex_consistency()` and converts any
/// `ValidationError::InconsistentVertex` into a defect record.
pub fn check_vertex_consistency(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_vertex_consistency() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that every edge is manifold (exactly 2 half-edges).
///
/// Calls `brep.validate_manifold_edges()` and converts any
/// `ValidationError::NonManifoldEdge` into a defect record.
pub fn check_manifold_edges(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_manifold_edges() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that face surface normals are not degenerate.
///
/// Calls `brep.validate_face_orientation()` and converts any
/// `ValidationError::InconsistentOrientation` into a defect record.
/// Uses Warning severity since orientation issues are often repairable.
pub fn check_face_orientation(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_face_orientation() {
        report.push(validation_error_to_defect(&e, Severity::Warning));
    }
}

/// Check that the shell is watertight (no gaps at shared edges).
///
/// Calls `brep.validate_watertight()` and converts any
/// `ValidationError::ShellGap` into a defect record.
pub fn check_shell_closure(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_watertight() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Convert a `ValidationError` into a `DefectRecord`.
///
/// Extracts the affected entity handle from the error variant and
/// uses the error's `Display` string as the human-readable description.
fn validation_error_to_defect(error: &ValidationError, severity: Severity) -> DefectRecord {
    let entities = extract_entities(error);
    let (measured_value, threshold) = extract_measurements(error);

    DefectRecord {
        category: DefectCategory::Topological,
        severity,
        description: error.to_string(),
        entities,
        location: None,
        measured_value,
        threshold,
        remediation: None,
    }
}

/// Extract entity references from a `ValidationError` variant.
fn extract_entities(error: &ValidationError) -> Vec<EntityRef> {
    match error {
        ValidationError::EulerPoincare { .. } => Vec::new(),
        ValidationError::BrokenTwin { half_edge } => {
            vec![EntityRef::HalfEdge(*half_edge)]
        }
        ValidationError::BrokenCycle {
            half_edge, loop_, ..
        } => {
            vec![EntityRef::HalfEdge(*half_edge), EntityRef::Loop(*loop_)]
        }
        ValidationError::InconsistentVertex {
            vertex, half_edge, ..
        } => {
            vec![EntityRef::Vertex(*vertex), EntityRef::HalfEdge(*half_edge)]
        }
        ValidationError::NonManifoldEdge { edge } => {
            vec![EntityRef::Edge(*edge)]
        }
        ValidationError::LowValence { vertex, .. } => {
            vec![EntityRef::Vertex(*vertex)]
        }
        ValidationError::InconsistentOrientation { face } => {
            vec![EntityRef::Face(*face)]
        }
        ValidationError::OpenLoop { loop_, .. } => {
            vec![EntityRef::Loop(*loop_)]
        }
        ValidationError::CurveSurfaceDeviation { edge, .. } => {
            vec![EntityRef::Edge(*edge)]
        }
        ValidationError::ShellGap { edge, .. } => {
            vec![EntityRef::Edge(*edge)]
        }
        ValidationError::SelfIntersection { face_a, face_b } => {
            vec![EntityRef::FacePair(*face_a, *face_b)]
        }
        ValidationError::TwistedSurface { face } => {
            vec![EntityRef::Face(*face)]
        }
        ValidationError::Precision { .. } => Vec::new(),
    }
}

/// Extract measured value and threshold from measurement-carrying error variants.
fn extract_measurements(error: &ValidationError) -> (Option<f64>, Option<f64>) {
    match error {
        ValidationError::OpenLoop { gap, .. } => (Some(*gap), None),
        ValidationError::CurveSurfaceDeviation { max_deviation, .. } => {
            (Some(*max_deviation), None)
        }
        ValidationError::ShellGap { gap, .. } => (Some(*gap), None),
        _ => (None, None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::healing::diagnostics::Severity;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn valid_cube_produces_no_topological_defects() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let mut report = HealingReport::new();
        check_all(&brep, &tol, &mut report);
        assert!(
            report.is_empty(),
            "Valid cube should have no topological defects, got {}",
            report.len()
        );
    }

    #[test]
    fn broken_twin_produces_defect() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let he_handle = brep.half_edges.keys().next().unwrap();
        brep.half_edges[he_handle].twin = Default::default();

        let mut report = HealingReport::new();
        check_half_edge_pairing(&brep, &mut report);
        assert!(!report.is_empty(), "Broken twin should produce a defect");
        assert_eq!(report.defects()[0].category, DefectCategory::Topological);
        assert_eq!(report.defects()[0].severity, Severity::Error);
    }

    #[test]
    fn corrupted_vertex_produces_defect() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let v_handle = brep.vertices.keys().next().unwrap();
        brep.vertices[v_handle].position = Point3::new(999.0, 999.0, 999.0);

        let mut report = HealingReport::new();
        check_shell_closure(&brep, &mut report);
        assert!(
            !report.is_empty(),
            "Corrupted vertex should trigger a shell closure defect"
        );
    }

    #[test]
    fn orientation_defect_is_warning() {
        use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face_handle = brep.faces.keys().next().unwrap();
        brep.faces[face_handle].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 0.0,
            },
            orientation: true,
        });

        let mut report = HealingReport::new();
        check_face_orientation(&brep, &mut report);
        assert!(!report.is_empty());
        assert_eq!(report.defects()[0].severity, Severity::Warning);
    }

    #[test]
    fn extract_entities_covers_all_variants() {
        // BrokenTwin
        let e = ValidationError::BrokenTwin {
            half_edge: Default::default(),
        };
        let ents = extract_entities(&e);
        assert_eq!(ents.len(), 1);
        assert!(matches!(ents[0], EntityRef::HalfEdge(_)));

        // NonManifoldEdge
        let e = ValidationError::NonManifoldEdge {
            edge: Default::default(),
        };
        let ents = extract_entities(&e);
        assert!(matches!(ents[0], EntityRef::Edge(_)));

        // InconsistentOrientation
        let e = ValidationError::InconsistentOrientation {
            face: Default::default(),
        };
        let ents = extract_entities(&e);
        assert!(matches!(ents[0], EntityRef::Face(_)));

        // EulerPoincare (no entities)
        let e = ValidationError::EulerPoincare {
            solid: Default::default(),
            expected: 0,
            actual: 2,
        };
        let ents = extract_entities(&e);
        assert!(ents.is_empty());
    }

    #[test]
    fn measurements_extracted_from_gap_error() {
        let e = ValidationError::ShellGap {
            edge: Default::default(),
            gap: 0.05,
        };
        let (measured, _threshold) = extract_measurements(&e);
        assert_eq!(measured, Some(0.05));
    }

    #[test]
    fn empty_brep_no_defects() {
        let brep = BRep::new(Tolerance::default());
        let mut report = HealingReport::new();
        check_all(&brep, &Tolerance::default(), &mut report);
        assert!(report.is_empty());
    }
}
