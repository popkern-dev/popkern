//! Import-specific defect detection (Section 5.4).
//!
//! Detects defects that commonly occur during STEP/IGES import:
//! lost periodic topology (seam edges missing on cylinders, tori, etc.)
//! and orientation inconsistencies introduced by translator bugs.

use crate::geometry::surface::SurfaceGeometry;
use crate::healing::diagnostics::{
    DefectCategory, DefectRecord, EntityRef, HealingReport, Remediation, RemediationAction,
    Severity,
};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Run all import-specific defect checks.
pub fn check_all(brep: &BRep, tolerance: &Tolerance, report: &mut HealingReport) {
    check_lost_periodic_topology(brep, report);
    check_orientation_inconsistencies(brep, report);
    check_lost_singularity_topology(brep, report);
    check_coordinate_artifacts(brep, tolerance, report);
    check_tolerance_mismatch(brep, tolerance, report);
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

/// Check for lost singularity topology at surface poles (Section 5.4).
///
/// Surfaces with poles (sphere, cone) should have zero-length degenerate
/// edges at those poles. If the pole is within the face's parameter range
/// but no degenerate edge exists, the topology may have been corrupted
/// during import.
///
/// This is a heuristic check: not all spherical/conical faces necessarily
/// cover a pole. The check looks for faces whose outer loop vertex positions
/// are close to the pole location but lack a zero-length edge there.
pub fn check_lost_singularity_topology(brep: &BRep, report: &mut HealingReport) {
    for (face_handle, face_data) in &brep.faces {
        let Some(ref surface_on_face) = face_data.surface else {
            continue;
        };

        // Identify pole locations based on surface type.
        let poles: Vec<crate::geometry::Point3> = match &surface_on_face.geometry {
            SurfaceGeometry::Sphere { center, radius } => {
                // Sphere poles at top (v=PI) and bottom (v=0) along Z axis.
                vec![
                    crate::geometry::Point3::new(center.x, center.y, center.z - radius),
                    crate::geometry::Point3::new(center.x, center.y, center.z + radius),
                ]
            }
            SurfaceGeometry::Cone { apex, .. } => {
                // Cone singularity at the apex.
                vec![*apex]
            }
            _ => continue,
        };

        // Check if any outer loop vertex is near a pole.
        let positions = brep.face_outer_loop_positions(face_handle);
        let has_zero_length_edge = has_zero_length_edge_in_loop(brep, face_data.outer_loop);

        for pole in &poles {
            let near_pole = positions.iter().any(|p| p.distance_to(pole) < 1e-6);
            if near_pole && !has_zero_length_edge {
                report.push(DefectRecord {
                    category: DefectCategory::Import,
                    severity: Severity::Warning,
                    description: format!(
                        "Face near surface pole at ({:.4}, {:.4}, {:.4}) \
                         but has no degenerate edge; singularity topology \
                         may have been lost during import",
                        pole.x, pole.y, pole.z
                    ),
                    entities: vec![EntityRef::Face(face_handle)],
                    location: None,
                    measured_value: None,
                    threshold: None,
                    remediation: Some(Remediation {
                        action: RemediationAction::RestoreSingularityTopology,
                        parameter: None,
                        value: None,
                        rationale: "Reconstruct degenerate edge at surface pole".into(),
                    }),
                });
            }
        }
    }
}

/// Helper: check if a loop contains any zero-length edge (vertex distance < 1e-10).
fn has_zero_length_edge_in_loop(brep: &BRep, loop_handle: crate::topology::LoopHandle) -> bool {
    let Some(loop_data) = brep.loops.get(loop_handle) else {
        return false;
    };

    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return false;
    }

    let mut cursor = start_he;
    let bound = brep.half_edges.len().max(1);

    for _ in 0..=bound {
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            break;
        };

        if let (Some(v1), Some(v2)) = (brep.vertices.get(he.origin), brep.vertices.get(twin.origin))
            && v1.position.distance_to(&v2.position) < 1e-10
        {
            return true;
        }

        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }

    false
}

/// Check for tolerance mismatch between model entities and BRep tolerance (Section 5.4).
///
/// Flags entities whose tracked tolerance is significantly larger or
/// smaller than the BRep's linear tolerance, suggesting the model was
/// imported from a system with different precision conventions.
pub fn check_tolerance_mismatch(brep: &BRep, tolerance: &Tolerance, report: &mut HealingReport) {
    let linear = tolerance.linear();

    for (vertex_handle, vertex) in &brep.vertices {
        if let Some(tol) = vertex.tolerance {
            if tol > linear * 100.0 {
                report.push(DefectRecord {
                    category: DefectCategory::Import,
                    severity: Severity::Warning,
                    description: format!(
                        "Entity tolerance {tol:.2e} is much larger than model tolerance \
                         {linear:.2e} — possible import from a coarser system"
                    ),
                    entities: vec![EntityRef::Vertex(vertex_handle)],
                    location: None,
                    measured_value: Some(tol),
                    threshold: Some(linear),
                    remediation: None,
                });
            } else if tol < linear * 0.01 {
                report.push(DefectRecord {
                    category: DefectCategory::Import,
                    severity: Severity::Info,
                    description: format!(
                        "Entity tolerance {tol:.2e} is much tighter than model tolerance \
                         {linear:.2e} — possible import from a finer system"
                    ),
                    entities: vec![EntityRef::Vertex(vertex_handle)],
                    location: None,
                    measured_value: Some(tol),
                    threshold: Some(linear),
                    remediation: None,
                });
            }
        }
    }

    for (edge_handle, edge) in &brep.edges {
        if let Some(tol) = edge.tolerance {
            if tol > linear * 100.0 {
                report.push(DefectRecord {
                    category: DefectCategory::Import,
                    severity: Severity::Warning,
                    description: format!(
                        "Entity tolerance {tol:.2e} is much larger than model tolerance \
                         {linear:.2e} — possible import from a coarser system"
                    ),
                    entities: vec![EntityRef::Edge(edge_handle)],
                    location: None,
                    measured_value: Some(tol),
                    threshold: Some(linear),
                    remediation: None,
                });
            } else if tol < linear * 0.01 {
                report.push(DefectRecord {
                    category: DefectCategory::Import,
                    severity: Severity::Info,
                    description: format!(
                        "Entity tolerance {tol:.2e} is much tighter than model tolerance \
                         {linear:.2e} — possible import from a finer system"
                    ),
                    entities: vec![EntityRef::Edge(edge_handle)],
                    location: None,
                    measured_value: Some(tol),
                    threshold: Some(linear),
                    remediation: None,
                });
            }
        }
    }

    for (face_handle, face) in &brep.faces {
        if let Some(tol) = face.tolerance {
            if tol > linear * 100.0 {
                report.push(DefectRecord {
                    category: DefectCategory::Import,
                    severity: Severity::Warning,
                    description: format!(
                        "Entity tolerance {tol:.2e} is much larger than model tolerance \
                         {linear:.2e} — possible import from a coarser system"
                    ),
                    entities: vec![EntityRef::Face(face_handle)],
                    location: None,
                    measured_value: Some(tol),
                    threshold: Some(linear),
                    remediation: None,
                });
            } else if tol < linear * 0.01 {
                report.push(DefectRecord {
                    category: DefectCategory::Import,
                    severity: Severity::Info,
                    description: format!(
                        "Entity tolerance {tol:.2e} is much tighter than model tolerance \
                         {linear:.2e} — possible import from a finer system"
                    ),
                    entities: vec![EntityRef::Face(face_handle)],
                    location: None,
                    measured_value: Some(tol),
                    threshold: Some(linear),
                    remediation: None,
                });
            }
        }
    }
}

/// Check for coordinate system artifacts (Section 5.4).
///
/// Detects vertices that are suspiciously close to zero in all coordinates
/// but not exactly zero, suggesting floating-point errors from coordinate
/// transformations during import. These near-origin vertices may indicate
/// that a transform collapsed geometry to the origin.
pub fn check_coordinate_artifacts(brep: &BRep, tolerance: &Tolerance, report: &mut HealingReport) {
    let linear_tol = tolerance.linear();

    for (vertex_handle, vertex) in &brep.vertices {
        let x = vertex.position.x.abs();
        let y = vertex.position.y.abs();
        let z = vertex.position.z.abs();

        // All coordinates are in (1e-15, linear_tolerance) — close to
        // zero but not exactly zero, suggesting a transform artifact.
        if x > 1e-15 && x < linear_tol && y > 1e-15 && y < linear_tol && z > 1e-15 && z < linear_tol
        {
            let dist = (x * x + y * y + z * z).sqrt();
            report.push(DefectRecord {
                category: DefectCategory::Import,
                severity: Severity::Info,
                description: format!(
                    "Vertex at ({:.2e}, {:.2e}, {:.2e}) is suspiciously \
                     close to origin — possible coordinate transform artifact",
                    vertex.position.x, vertex.position.y, vertex.position.z
                ),
                entities: vec![EntityRef::Vertex(vertex_handle)],
                location: None,
                measured_value: Some(dist),
                threshold: Some(linear_tol),
                remediation: Some(Remediation {
                    action: RemediationAction::AdjustGeometry,
                    parameter: None,
                    value: None,
                    rationale: "Snap vertex to exact origin if intended to be \
                                at the coordinate system origin"
                        .into(),
                }),
            });
        }
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

    #[test]
    fn lost_singularity_topology_on_cube_produces_no_defect() {
        // A cube has planar faces; no sphere/cone surfaces.
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut report = HealingReport::new();
        check_lost_singularity_topology(&brep, &mut report);
        assert!(
            report.is_empty(),
            "Cube with planar faces should have no singularity topology defects"
        );
    }

    #[test]
    fn coordinate_artifacts_on_valid_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let mut report = HealingReport::new();
        check_coordinate_artifacts(&brep, &tol, &mut report);
        // Unit cube vertices are at 0.0 and 1.0 — none in the suspicious range.
        assert!(
            report.is_empty(),
            "Valid cube should have no coordinate artifact defects"
        );
    }

    #[test]
    fn empty_brep_no_singularity_or_artifact_defects() {
        let brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let mut report = HealingReport::new();
        check_lost_singularity_topology(&brep, &mut report);
        check_coordinate_artifacts(&brep, &tol, &mut report);
        assert!(report.is_empty());
    }
}
