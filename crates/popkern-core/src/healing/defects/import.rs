//! Import-specific defect detection (Section 5.4).
//!
//! Detects defects that commonly occur during STEP/IGES import:
//! lost periodic topology (seam edges missing on cylinders, tori, etc.),
//! orientation inconsistencies introduced by translator bugs, and
//! representation discontinuities from format conversion artifacts.

use crate::geometry::curve::CurveGeometry;
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
    check_representation_discontinuities(brep, tolerance, report);
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

    // Collect (entity_ref, entity_tolerance) triples.
    let entities: Vec<(EntityRef, f64)> = brep
        .vertices
        .iter()
        .filter_map(|(h, v)| v.tolerance.map(|t| (EntityRef::Vertex(h), t)))
        .chain(
            brep.edges
                .iter()
                .filter_map(|(h, e)| e.tolerance.map(|t| (EntityRef::Edge(h), t))),
        )
        .chain(
            brep.faces
                .iter()
                .filter_map(|(h, f)| f.tolerance.map(|t| (EntityRef::Face(h), t))),
        )
        .collect();

    for (entity, tol) in entities {
        if tol > linear * 100.0 {
            report.push(DefectRecord {
                category: DefectCategory::Import,
                severity: Severity::Warning,
                description: format!(
                    "Entity tolerance {tol:.2e} is much larger than model tolerance \
                     {linear:.2e} — possible import from a coarser system"
                ),
                entities: vec![entity],
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
                entities: vec![entity],
                location: None,
                measured_value: Some(tol),
                threshold: Some(linear),
                remediation: None,
            });
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

/// Check for representation discontinuities from format conversion (Section 5.4).
///
/// Detects NURBS curves and surfaces that have internal continuity
/// breaks at knot boundaries, which may indicate conversion artifacts
/// from rational-to-polynomial approximation or representation changes.
pub fn check_representation_discontinuities(
    brep: &BRep,
    _tolerance: &Tolerance,
    report: &mut HealingReport,
) {
    for (edge_handle, edge) in &brep.edges {
        if let Some(ref coe) = edge.curve
            && let CurveGeometry::Nurbs(nurbs) = &coe.geometry
        {
            emit_knot_defects(
                &nurbs.knots,
                nurbs.degree,
                "NURBS curve",
                "",
                EntityRef::Edge(edge_handle),
                report,
            );
        }
    }
    for (face_handle, face) in &brep.faces {
        if let Some(ref sof) = face.surface
            && let SurfaceGeometry::Nurbs(nurbs) = &sof.geometry
        {
            let entity = EntityRef::Face(face_handle);
            emit_knot_defects(
                &nurbs.knots_u,
                nurbs.degree_u,
                "NURBS surface",
                " in u-direction",
                entity.clone(),
                report,
            );
            emit_knot_defects(
                &nurbs.knots_v,
                nurbs.degree_v,
                "NURBS surface",
                " in v-direction",
                entity,
                report,
            );
        }
    }
}

/// Emit defect records for high-multiplicity interior knots.
fn emit_knot_defects(
    knots: &[f64],
    degree: usize,
    entity_label: &str,
    dir_suffix: &str,
    entity: EntityRef,
    report: &mut HealingReport,
) {
    for (knot_idx, mult) in knot_multiplicities(knots, degree) {
        if mult >= degree {
            report.push(DefectRecord {
                category: DefectCategory::Import,
                severity: Severity::Warning,
                description: format!(
                    "{entity_label} has C0 discontinuity{dir_suffix} at interior \
                     knot index {knot_idx} (multiplicity {mult} = degree {degree}); \
                     possible conversion artifact"
                ),
                entities: vec![entity.clone()],
                location: None,
                measured_value: Some(mult as f64),
                threshold: Some(degree as f64),
                remediation: Some(Remediation {
                    action: RemediationAction::AdjustGeometry,
                    parameter: None,
                    value: None,
                    rationale: format!(
                        "Remove or reduce knot multiplicity to restore \
                         smoothness{dir_suffix}"
                    ),
                }),
            });
        } else if mult + 1 >= degree && degree >= 2 {
            report.push(DefectRecord {
                category: DefectCategory::Import,
                severity: Severity::Info,
                description: format!(
                    "{entity_label} has C1 discontinuity{dir_suffix} at interior \
                     knot index {knot_idx} (multiplicity {mult} = degree {degree} - 1); \
                     possible conversion artifact"
                ),
                entities: vec![entity.clone()],
                location: None,
                measured_value: Some(mult as f64),
                threshold: Some(degree as f64),
                remediation: None,
            });
        }
    }
}

/// Compute the multiplicity of each distinct interior knot.
///
/// Interior knots are those not in the first `degree + 1` or last
/// `degree + 1` positions of the knot vector. Returns `(knot_index,
/// multiplicity)` pairs for each distinct interior knot.
fn knot_multiplicities(knots: &[f64], degree: usize) -> Vec<(usize, usize)> {
    let n = knots.len();
    let boundary = degree + 1;
    if n <= 2 * boundary || boundary >= n - boundary {
        return Vec::new();
    }
    let (start, end) = (boundary, n - boundary);
    let mut result = Vec::new();
    let mut i = start;
    while i < end {
        let val = knots[i];
        let mut m = 1;
        while i + m < end && (knots[i + m] - val).abs() < 1e-14 {
            m += 1;
        }
        result.push((i, m));
        i += m;
    }
    result
}

#[cfg(test)]
#[path = "import_tests.rs"]
mod tests;
