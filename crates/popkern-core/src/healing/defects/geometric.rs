//! Geometric defect detection (Section 5.2).
//!
//! Checks for edge-surface consistency (deviation between an edge's curve
//! and the adjacent face's surface) and tolerance hierarchy violations
//! (vertex tolerance must be <= edge tolerance <= face tolerance).

use crate::error::ValidationError;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::uv_projection::tight_bounded_parameter_range;
use crate::healing::curvature::fundamental_forms::first_fundamental_form;
use crate::healing::deviation::edge_surface_deviation;
use crate::healing::diagnostics::{
    DefectCategory, DefectRecord, EntityRef, HealingReport, Remediation, RemediationAction,
    Severity,
};
use crate::healing::options::SampleDensity;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Run all geometric defect checks.
pub fn check_all(
    brep: &BRep,
    tolerance: &Tolerance,
    sample_density: &SampleDensity,
    report: &mut HealingReport,
) {
    check_edge_surface_consistency(brep, tolerance, sample_density, report);
    check_tolerance_hierarchy(brep, report);
    check_self_intersection(brep, report);
    check_parameterization_quality(brep, tolerance, report);
    super::parameterization::check_self_overlapping_parameterization(brep, report);
    super::containment::check_surface_containment(brep, report);
    super::containment::check_boundary_containment(brep, tolerance, sample_density, report);
}

/// Check that each edge's curve is consistent with adjacent face surfaces.
///
/// For each edge that has curve geometry, iterates over its half-edges.
/// For each half-edge whose parent face has surface geometry, measures
/// the deviation between the edge curve and the face surface. If the
/// maximum deviation exceeds the linear tolerance, a Geometric defect
/// is recorded with a `RecomputeEdgeCurve` remediation.
pub fn check_edge_surface_consistency(
    brep: &BRep,
    tolerance: &Tolerance,
    sample_density: &SampleDensity,
    report: &mut HealingReport,
) {
    let num_samples = sample_density.samples_per_edge();

    for (edge_handle, edge_data) in &brep.edges {
        let Some(ref curve_on_edge) = edge_data.curve else {
            continue;
        };

        let t_range = (curve_on_edge.t_start, curve_on_edge.t_end);

        // Check against both adjacent faces (via the two half-edges).
        let twin_handle = brep.half_edges.get(edge_data.half_edge).map(|he| he.twin);
        let half_edges_to_check = [Some(edge_data.half_edge), twin_handle];

        for maybe_he_handle in half_edges_to_check.iter().flatten() {
            let Some(he) = brep.half_edges.get(*maybe_he_handle) else {
                continue;
            };
            let Some(loop_data) = brep.loops.get(he.loop_) else {
                continue;
            };
            let Some(face_data) = brep.faces.get(loop_data.face) else {
                continue;
            };
            let Some(ref surface_on_face) = face_data.surface else {
                continue;
            };

            let result = edge_surface_deviation(
                &curve_on_edge.geometry,
                t_range,
                &surface_on_face.geometry,
                num_samples,
            );

            if result.max_deviation > tolerance.linear() {
                report.push(DefectRecord {
                    category: DefectCategory::Geometric,
                    severity: Severity::Error,
                    description: format!(
                        "Edge curve deviates from adjacent face surface by {:.6} \
                         (tolerance: {:.6})",
                        result.max_deviation,
                        tolerance.linear()
                    ),
                    entities: vec![EntityRef::Edge(edge_handle)],
                    location: None,
                    measured_value: Some(result.max_deviation),
                    threshold: Some(tolerance.linear()),
                    remediation: Some(Remediation {
                        action: RemediationAction::RecomputeEdgeCurve,
                        parameter: None,
                        value: None,
                        rationale: "Re-derive edge curve from adjacent surfaces \
                                    to eliminate deviation"
                            .into(),
                    }),
                });
            }
        }
    }
}

/// Check the tolerance hierarchy: vertex <= edge <= face.
///
/// For each edge, if both the vertex and edge have explicit tolerances
/// and the vertex tolerance exceeds the edge tolerance, a Warning is
/// recorded. Similarly, for each face, if an adjacent edge's tolerance
/// exceeds the face tolerance, a Warning is recorded.
///
/// The correct hierarchy is: vertex tolerance <= edge tolerance <= face
/// tolerance. Violations indicate the model was built with inconsistent
/// precision metadata.
pub fn check_tolerance_hierarchy(brep: &BRep, report: &mut HealingReport) {
    // Check vertex <= edge for each edge's origin vertices.
    for (edge_handle, edge_data) in &brep.edges {
        let Some(edge_tol) = edge_data.tolerance else {
            continue;
        };

        let Some(he) = brep.half_edges.get(edge_data.half_edge) else {
            continue;
        };

        // Check origin vertex of canonical half-edge.
        check_vertex_edge_tolerance(brep, he.origin, edge_handle, edge_tol, report);

        // Check origin vertex of twin half-edge.
        if let Some(twin) = brep.half_edges.get(he.twin) {
            check_vertex_edge_tolerance(brep, twin.origin, edge_handle, edge_tol, report);
        }
    }

    // Check edge <= face for each face's boundary edges.
    for (face_handle, face_data) in &brep.faces {
        let Some(face_tol) = face_data.tolerance else {
            continue;
        };

        // Walk the outer loop to find all edges on the face boundary.
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };

        let start_he = loop_data.half_edge;
        let mut cursor = start_he;
        let bound = brep.half_edges.len().max(1);

        for _ in 0..=bound {
            let Some(he) = brep.half_edges.get(cursor) else {
                break;
            };

            if let Some(edge_data) = brep.edges.get(he.edge)
                && let Some(edge_tol) = edge_data.tolerance
                && edge_tol > face_tol
            {
                report.push(DefectRecord {
                    category: DefectCategory::Tolerance,
                    severity: Severity::Warning,
                    description: format!(
                        "Edge tolerance ({:.6}) exceeds face tolerance ({:.6})",
                        edge_tol, face_tol
                    ),
                    entities: vec![EntityRef::Edge(he.edge), EntityRef::Face(face_handle)],
                    location: None,
                    measured_value: Some(edge_tol),
                    threshold: Some(face_tol),
                    remediation: Some(Remediation {
                        action: RemediationAction::EnforceToleranceHierarchy,
                        parameter: None,
                        value: None,
                        rationale: "Edge tolerance should not exceed \
                                    the parent face tolerance"
                            .into(),
                    }),
                });
            }

            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }
    }
}

/// Helper: check that a vertex tolerance does not exceed an edge tolerance.
fn check_vertex_edge_tolerance(
    brep: &BRep,
    vertex_handle: crate::topology::VertexHandle,
    edge_handle: crate::topology::EdgeHandle,
    edge_tol: f64,
    report: &mut HealingReport,
) {
    let Some(vertex) = brep.vertices.get(vertex_handle) else {
        return;
    };
    let Some(vertex_tol) = vertex.tolerance else {
        return;
    };

    if vertex_tol > edge_tol {
        report.push(DefectRecord {
            category: DefectCategory::Tolerance,
            severity: Severity::Warning,
            description: format!(
                "Vertex tolerance ({:.6}) exceeds edge tolerance ({:.6})",
                vertex_tol, edge_tol
            ),
            entities: vec![
                EntityRef::Vertex(vertex_handle),
                EntityRef::Edge(edge_handle),
            ],
            location: None,
            measured_value: Some(vertex_tol),
            threshold: Some(edge_tol),
            remediation: Some(Remediation {
                action: RemediationAction::EnforceToleranceHierarchy,
                parameter: None,
                value: None,
                rationale: "Vertex tolerance should not exceed the parent edge tolerance".into(),
            }),
        });
    }
}

/// Check for self-intersecting faces (Section 5.2.4).
///
/// Wraps the existing `validate_self_intersection()` method, converting
/// the result to a [`DefectRecord`]. Only runs on BReps with at least
/// two faces (same precondition as `validate_all`).
pub fn check_self_intersection(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_self_intersection() {
        let entities = match &e {
            ValidationError::SelfIntersection { face_a, face_b } => {
                vec![EntityRef::FacePair(*face_a, *face_b)]
            }
            _ => Vec::new(),
        };
        report.push(DefectRecord {
            category: DefectCategory::Geometric,
            severity: Severity::Error,
            description: e.to_string(),
            entities,
            location: None,
            measured_value: None,
            threshold: None,
            remediation: Some(Remediation {
                action: RemediationAction::SplitEntity,
                parameter: None,
                value: None,
                rationale: "Split surface at self-intersection curve".into(),
            }),
        });
    }
}

/// Check for parameterization quality issues on face surfaces (Section 4.2).
///
/// Detects severely distorted parameterizations where the area element
/// ratio (max/min) exceeds a threshold, indicating that the UV mapping
/// is highly non-uniform.
pub fn check_parameterization_quality(
    brep: &BRep,
    tolerance: &Tolerance,
    report: &mut HealingReport,
) {
    const GRID_SIZE: usize = 4;
    const DISTORTION_THRESHOLD: f64 = 1000.0;

    let degenerate_threshold = tolerance.linear() * tolerance.linear();

    for (face_handle, face_data) in &brep.faces {
        let Some(ref surface_on_face) = face_data.surface else {
            continue;
        };

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
        let mut min_area = f64::INFINITY;
        let mut max_area = 0.0_f64;
        let mut has_near_degenerate = false;

        for i in 0..GRID_SIZE {
            let u = u_min + u_span * (i as f64 + 0.5) / GRID_SIZE as f64;
            for j in 0..GRID_SIZE {
                let v = v_min + v_span * (j as f64 + 0.5) / GRID_SIZE as f64;
                let area = first_fundamental_form(surface, u, v).area_element();

                if area < degenerate_threshold {
                    has_near_degenerate = true;
                    continue; // skip near-zero for ratio computation
                }

                min_area = min_area.min(area);
                max_area = max_area.max(area);
            }
        }

        if has_near_degenerate {
            report.push(DefectRecord {
                category: DefectCategory::Geometric,
                severity: Severity::Warning,
                description: "Near-degenerate parameterization point detected".into(),
                entities: vec![EntityRef::Face(face_handle)],
                location: None,
                measured_value: None,
                threshold: Some(degenerate_threshold),
                remediation: Some(Remediation {
                    action: RemediationAction::ReparameterizeSurface,
                    parameter: None,
                    value: None,
                    rationale: "Reparameterize surface to remove \
                                near-degenerate UV mapping"
                        .into(),
                }),
            });
        }

        if min_area > 0.0 && max_area / min_area > DISTORTION_THRESHOLD {
            report.push(DefectRecord {
                category: DefectCategory::Geometric,
                severity: Severity::Warning,
                description: format!(
                    "Severely distorted parameterization: area element ratio {:.1}",
                    max_area / min_area
                ),
                entities: vec![EntityRef::Face(face_handle)],
                location: None,
                measured_value: Some(max_area / min_area),
                threshold: Some(DISTORTION_THRESHOLD),
                remediation: Some(Remediation {
                    action: RemediationAction::ReparameterizeSurface,
                    parameter: None,
                    value: None,
                    rationale: "Reparameterize surface to reduce \
                                UV distortion ratio"
                        .into(),
                }),
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::healing::diagnostics::DefectCategory;
    use crate::healing::options::SampleDensity;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn valid_cube_no_geometric_defects() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let density = SampleDensity::Standard;
        let mut report = HealingReport::new();
        check_all(&brep, &tol, &density, &mut report);
        assert!(
            report.is_empty(),
            "Valid cube should have no geometric defects, got {}",
            report.len()
        );
    }

    #[test]
    fn edge_surface_consistency_on_cube_with_geometry() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let density = SampleDensity::Fine;
        let mut report = HealingReport::new();
        check_edge_surface_consistency(&brep, &tol, &density, &mut report);
        assert!(
            report.is_empty(),
            "Unit cube edges should be consistent with face surfaces"
        );
    }

    #[test]
    fn tolerance_hierarchy_no_violations_on_default_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut report = HealingReport::new();
        check_tolerance_hierarchy(&brep, &mut report);
        // Default cube has no per-entity tolerances, so no violations.
        assert!(report.is_empty());
    }

    #[test]
    fn vertex_tolerance_exceeding_edge_produces_defect() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();

        // Set an edge tolerance, then set its vertex tolerance higher.
        let edge_handle = brep.edges.keys().next().unwrap();
        brep.edges[edge_handle].tolerance = Some(1e-6);

        let he_handle = brep.edges[edge_handle].half_edge;
        let vertex_handle = brep.half_edges[he_handle].origin;
        brep.vertices[vertex_handle].tolerance = Some(1e-3);

        let mut report = HealingReport::new();
        check_tolerance_hierarchy(&brep, &mut report);
        assert!(
            !report.is_empty(),
            "Vertex tolerance > edge tolerance should produce a defect"
        );
        assert_eq!(report.defects()[0].category, DefectCategory::Tolerance);
    }

    #[test]
    fn edge_tolerance_exceeding_face_produces_defect() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();

        // Set a face tolerance, then set a boundary edge's tolerance higher.
        let face_handle = brep.faces.keys().next().unwrap();
        brep.faces[face_handle].tolerance = Some(1e-6);

        let outer_loop = brep.faces[face_handle].outer_loop;
        let he_handle = brep.loops[outer_loop].half_edge;
        let edge_handle = brep.half_edges[he_handle].edge;
        brep.edges[edge_handle].tolerance = Some(1e-3);

        let mut report = HealingReport::new();
        check_tolerance_hierarchy(&brep, &mut report);
        assert!(
            !report.is_empty(),
            "Edge tolerance > face tolerance should produce a defect"
        );
    }

    #[test]
    fn self_intersection_on_valid_cube_produces_no_defect() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut report = HealingReport::new();
        check_self_intersection(&brep, &mut report);
        assert!(
            report.is_empty(),
            "Valid cube should have no self-intersection defects"
        );
    }

    #[test]
    fn empty_brep_no_geometric_defects() {
        let brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let density = SampleDensity::Standard;
        let mut report = HealingReport::new();
        check_all(&brep, &tol, &density, &mut report);
        assert!(report.is_empty());
    }

    #[test]
    fn parameterization_quality_on_valid_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let mut report = HealingReport::new();
        check_parameterization_quality(&brep, &tol, &mut report);
        // Plane surfaces have uniform parameterization (area element = 1.0 everywhere).
        assert!(
            report.is_empty(),
            "Valid cube should have no parameterization defects, got {}",
            report.len()
        );
    }

    #[test]
    fn parameterization_quality_on_empty_brep() {
        let brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let mut report = HealingReport::new();
        check_parameterization_quality(&brep, &tol, &mut report);
        assert!(report.is_empty());
    }
}
