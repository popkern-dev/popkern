//! Self-overlapping parameterization detection (Section 4.2).
//!
//! A parameterization is self-overlapping when the surface-to-parameter
//! map is not locally injective — two different UV points map to the
//! same 3D point. This is detected by checking the sign of the Jacobian
//! determinant (the area element from the first fundamental form). If
//! the sign changes within a face, the parameterization folds over itself.

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::tight_bounded_parameter_range;
use crate::healing::curvature::fundamental_forms::first_fundamental_form;
use crate::healing::diagnostics::{
    DefectCategory, DefectRecord, EntityRef, HealingReport, Remediation, RemediationAction,
    Severity,
};
use crate::topology::BRep;

/// Grid size for Jacobian sampling.
const GRID_SIZE: usize = 6;

/// Check for self-overlapping parameterizations (Section 4.2).
///
/// A parameterization is self-overlapping when the surface-to-parameter
/// map is not locally injective — two different UV points map to the
/// same 3D point. This is detected by checking the sign of the Jacobian
/// determinant (the area element from the first fundamental form). If
/// the sign changes within a face, the parameterization folds over itself.
///
/// Only NURBS surfaces are checked — analytical surfaces (plane, cylinder,
/// sphere, cone, torus) have globally injective parameterizations by
/// construction.
pub fn check_self_overlapping_parameterization(brep: &BRep, report: &mut HealingReport) {
    for (face_handle, face_data) in &brep.faces {
        let Some(ref surface_on_face) = face_data.surface else {
            continue;
        };

        // Only NURBS surfaces can self-overlap; analytical surfaces are
        // injective by construction.
        if !matches!(surface_on_face.geometry, SurfaceGeometry::Nurbs(_)) {
            continue;
        }

        let positions = brep.face_outer_loop_positions(face_handle);
        if positions.is_empty() {
            continue;
        }

        let ((u_min, u_max), (v_min, v_max)) =
            tight_bounded_parameter_range(&positions, &surface_on_face.geometry);

        let u_span = u_max - u_min;
        let v_span = v_max - v_min;
        if u_span <= 0.0 || v_span <= 0.0 {
            continue;
        }

        let surface: &dyn SurfaceEval = &surface_on_face.geometry;
        let mut has_positive = false;
        let mut has_negative = false;

        for i in 0..GRID_SIZE {
            let u = u_min + u_span * (i as f64 + 0.5) / GRID_SIZE as f64;
            for j in 0..GRID_SIZE {
                let v = v_min + v_span * (j as f64 + 0.5) / GRID_SIZE as f64;
                let det = first_fundamental_form(surface, u, v).determinant();

                if det > 0.0 {
                    has_positive = true;
                } else if det < 0.0 {
                    has_negative = true;
                }

                // Early exit once both signs observed.
                if has_positive && has_negative {
                    break;
                }
            }
            if has_positive && has_negative {
                break;
            }
        }

        if has_positive && has_negative {
            report.push(DefectRecord {
                category: DefectCategory::Geometric,
                severity: Severity::Warning,
                description: "Self-overlapping parameterization on face \
                     — Jacobian sign reversal detected"
                    .into(),
                entities: vec![EntityRef::Face(face_handle)],
                location: None,
                measured_value: None,
                threshold: None,
                remediation: Some(Remediation {
                    action: RemediationAction::ReparameterizeSurface,
                    parameter: None,
                    value: None,
                    rationale: "Reparameterize surface to eliminate \
                                self-overlapping UV mapping"
                        .into(),
                }),
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    use crate::healing::diagnostics::HealingReport;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn valid_cube_no_self_overlapping_parameterization() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut report = HealingReport::new();
        check_self_overlapping_parameterization(&brep, &mut report);
        assert!(
            report.is_empty(),
            "Valid cube should have no self-overlapping parameterization defects, got {}",
            report.len()
        );
    }

    #[test]
    fn empty_brep_no_self_overlapping_parameterization() {
        let brep = BRep::new(Tolerance::default());
        let mut report = HealingReport::new();
        check_self_overlapping_parameterization(&brep, &mut report);
        assert!(report.is_empty());
    }

    #[test]
    fn nurbs_face_without_overlap_no_defect() {
        // Build a BRep with a NURBS face (plane converted to NURBS).
        // Should not trigger self-overlapping since NURBS-plane has
        // positive Jacobian everywhere.
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let fh = brep.faces.keys().next().unwrap();
        let positions = brep.face_outer_loop_positions(fh);
        assert!(!positions.is_empty());
        if let Some(ref surf) = brep.faces[fh].surface {
            let nurbs = surf.geometry.to_nurbs_in_range((0.0, 1.0), (0.0, 1.0));
            brep.faces[fh].surface = Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(nurbs),
                orientation: true,
            });
        }
        let mut report = HealingReport::new();
        check_self_overlapping_parameterization(&brep, &mut report);
        assert!(report.is_empty());
    }

    #[test]
    fn face_without_surface_skipped() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let fh = brep.faces.keys().next().unwrap();
        brep.faces[fh].surface = None;
        let mut report = HealingReport::new();
        check_self_overlapping_parameterization(&brep, &mut report);
        assert!(report.is_empty());
    }

    #[test]
    fn analytical_surface_skipped() {
        // Verify the "not NURBS -> continue" branch.
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut report = HealingReport::new();
        check_self_overlapping_parameterization(&brep, &mut report);
        // All planes -> skipped, no defects.
        assert!(report.is_empty());
    }
}
