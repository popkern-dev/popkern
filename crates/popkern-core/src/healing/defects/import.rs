//! Import-specific defect detection (Section 5.4).
//!
//! Detects defects that commonly occur during STEP/IGES import:
//! lost periodic topology (seam edges missing on cylinders, tori, etc.)
//! and orientation inconsistencies introduced by translator bugs.

use crate::healing::diagnostics::{
    DefectCategory, DefectRecord, EntityRef, HealingReport, Remediation, RemediationAction,
    Severity,
};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Run all import-specific defect checks.
pub fn check_all(brep: &BRep, _tolerance: &Tolerance, report: &mut HealingReport) {
    check_lost_periodic_topology(brep, report);
    check_orientation_inconsistencies(brep, report);
}

/// Check for lost periodic topology on faces with periodic surfaces.
///
/// Periodic surfaces (cylinder, torus, cone, sphere) require seam edges
/// where both half-edges belong to the same face. When import translators
/// lose this topology, the face cannot be properly trimmed and
/// downstream operations may fail.
///
/// For each face with a periodic surface, walks the outer loop and checks
/// whether any edge has both half-edges belonging to the same face (a seam
/// edge). If no seam edge is found on a face whose surface is periodic in
/// at least one direction, a Warning is recorded.
pub fn check_lost_periodic_topology(brep: &BRep, report: &mut HealingReport) {
    for (face_handle, face_data) in &brep.faces {
        let Some(ref surface_on_face) = face_data.surface else {
            continue;
        };

        let (periodic_u, periodic_v) = surface_on_face.geometry.is_periodic();
        if !periodic_u && !periodic_v {
            continue;
        }

        // Walk the outer loop looking for a seam edge.
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };

        let start_he = loop_data.half_edge;
        let mut cursor = start_he;
        let bound = brep.half_edges.len().max(1);
        let mut found_seam = false;

        for _ in 0..=bound {
            let Some(he) = brep.half_edges.get(cursor) else {
                break;
            };

            // A seam edge is one where both half-edges belong to the same face.
            if let Some(twin) = brep.half_edges.get(he.twin) {
                let twin_face = brep.loops.get(twin.loop_).map(|l| l.face);
                if twin_face == Some(face_handle) {
                    found_seam = true;
                    break;
                }
            }

            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }

        if !found_seam {
            let direction = match (periodic_u, periodic_v) {
                (true, true) => "u and v",
                (true, false) => "u",
                (false, true) => "v",
                (false, false) => unreachable!(),
            };

            report.push(DefectRecord {
                category: DefectCategory::Import,
                severity: Severity::Warning,
                description: format!(
                    "Periodic surface (periodic in {direction}) has no seam edge; \
                     topology may have been lost during import"
                ),
                entities: vec![EntityRef::Face(face_handle)],
                location: None,
                measured_value: None,
                threshold: None,
                remediation: Some(Remediation {
                    action: RemediationAction::RestorePeriodicTopology,
                    parameter: None,
                    value: None,
                    rationale: "Reconstruct missing seam edge to restore \
                                periodic surface closure"
                        .into(),
                }),
            });
        }
    }
}

/// Check for orientation inconsistencies that may have been introduced
/// during import.
///
/// Delegates to `brep.validate_face_orientation()` and converts any
/// error into an Import category defect. Import-induced orientation
/// issues are common in STEP files from certain CAD systems.
pub fn check_orientation_inconsistencies(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_face_orientation() {
        report.push(DefectRecord {
            category: DefectCategory::Import,
            severity: Severity::Warning,
            description: format!("Possible import-induced orientation inconsistency: {e}"),
            entities: match &e {
                crate::error::ValidationError::InconsistentOrientation { face } => {
                    vec![EntityRef::Face(*face)]
                }
                _ => Vec::new(),
            },
            location: None,
            measured_value: None,
            threshold: None,
            remediation: Some(Remediation {
                action: RemediationAction::AdjustGeometry,
                parameter: None,
                value: None,
                rationale: "Flip surface orientation or recompute face normals".into(),
            }),
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn valid_cube_no_import_defects() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let mut report = HealingReport::new();
        check_all(&brep, &tol, &mut report);
        // The unit cube has planar faces (non-periodic), so no lost
        // periodic topology, and valid orientation.
        assert!(
            report.is_empty(),
            "Valid cube should have no import defects, got {}",
            report.len()
        );
    }

    #[test]
    fn orientation_inconsistency_detected() {
        use crate::geometry::Point3;
        use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face_handle = brep.faces.keys().next().unwrap();
        // Attach a degenerate surface to trigger orientation failure.
        brep.faces[face_handle].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 0.0,
            },
            orientation: true,
        });

        let mut report = HealingReport::new();
        check_orientation_inconsistencies(&brep, &mut report);
        assert!(
            !report.is_empty(),
            "Degenerate surface should trigger orientation defect"
        );
        assert_eq!(report.defects()[0].category, DefectCategory::Import);
    }

    #[test]
    fn empty_brep_no_import_defects() {
        let brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let mut report = HealingReport::new();
        check_all(&brep, &tol, &mut report);
        assert!(report.is_empty());
    }

    #[test]
    fn non_periodic_face_skipped() {
        // A cube has only planar faces; none should trigger periodic checks.
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut report = HealingReport::new();
        check_lost_periodic_topology(&brep, &mut report);
        assert!(report.is_empty());
    }
}
