//! Boundary containment defect detection (Sections 5.2.2, 5.2.5).
//!
//! Checks that face trimming boundaries lie within the surface parameter
//! domain, and that trim curves (pcurves) agree with the corresponding
//! 3D edge curves when evaluated through the surface.

use crate::geometry::curve::{CurveEval, closest_point_on_curve};
use crate::geometry::surface::SurfaceEval;
use crate::healing::diagnostics::{
    DefectCategory, DefectRecord, EntityRef, HealingReport, Remediation, RemediationAction,
    Severity,
};
use crate::healing::options::SampleDensity;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Check that face trimming boundaries lie within the surface parameter domain
/// (Section 5.2.2).
///
/// For each face with a surface, evaluates the surface parameter range and
/// checks that pcurve sample points lie within bounds. Surfaces with
/// infinite parameter domains (planes, cylinders) are skipped because
/// their trim curves can legitimately occupy any region.
pub fn check_surface_containment(brep: &BRep, report: &mut HealingReport) {
    const SAMPLES_PER_PCURVE: usize = 4;
    // Relative padding fraction: pcurve points up to 5% outside the
    // native parameter range are tolerated.
    const PADDING_FRACTION: f64 = 0.05;

    for (face_handle, face_data) in &brep.faces {
        let Some(ref surface_on_face) = face_data.surface else {
            continue;
        };

        let ((u_min, u_max), (v_min, v_max)) = surface_on_face.geometry.parameter_range();

        // Skip surfaces with any infinite parameter dimension.
        if !u_min.is_finite() || !u_max.is_finite() || !v_min.is_finite() || !v_max.is_finite() {
            continue;
        }

        let u_pad = (u_max - u_min) * PADDING_FRACTION;
        let v_pad = (v_max - v_min) * PADDING_FRACTION;

        // Walk the outer loop half-edges.
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

            if let Some(ref pcurve) = he.pcurve {
                let (t0, t1) = pcurve.parameter_range();
                for k in 0..SAMPLES_PER_PCURVE {
                    let frac = (k as f64 + 0.5) / SAMPLES_PER_PCURVE as f64;
                    let t = t0 + (t1 - t0) * frac;
                    let pt = pcurve.point_at(t);
                    let u = pt.x;
                    let v = pt.y;

                    if u < u_min - u_pad
                        || u > u_max + u_pad
                        || v < v_min - v_pad
                        || v > v_max + v_pad
                    {
                        report.push(DefectRecord {
                            category: DefectCategory::Geometric,
                            severity: Severity::Warning,
                            description: format!(
                                "Pcurve point ({:.4}, {:.4}) outside surface \
                                 domain [{:.4}..{:.4}] x [{:.4}..{:.4}]",
                                u, v, u_min, u_max, v_min, v_max
                            ),
                            entities: vec![
                                EntityRef::Face(face_handle),
                                EntityRef::HalfEdge(cursor),
                            ],
                            location: None,
                            measured_value: None,
                            threshold: None,
                            remediation: Some(Remediation {
                                action: RemediationAction::RecomputePCurve,
                                parameter: None,
                                value: None,
                                rationale: "Recompute pcurve to lie within \
                                            the surface parameter domain"
                                    .into(),
                            }),
                        });
                        // One defect per half-edge is sufficient.
                        break;
                    }
                }
            }

            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }
    }
}

/// Check that all trim curves lie on their respective surfaces (Section 5.2.5).
///
/// For each half-edge with a pcurve, evaluates the pcurve at several points,
/// maps those UV coordinates to 3D via the surface, and checks that each
/// resulting 3D point lies on the edge's 3D curve within tolerance.
///
/// Uses Newton-based closest-point projection onto the edge curve rather
/// than assuming parameter correspondence between pcurve and edge curve.
pub fn check_boundary_containment(
    brep: &BRep,
    tolerance: &Tolerance,
    _sample_density: &SampleDensity,
    report: &mut HealingReport,
) {
    const PC_SAMPLES: usize = 4;

    for (face_handle, face_data) in &brep.faces {
        let Some(ref surface_on_face) = face_data.surface else {
            continue;
        };

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

            if let (Some(pcurve), Some(edge_data)) = (&he.pcurve, brep.edges.get(he.edge))
                && let Some(curve_on_edge) = &edge_data.curve
            {
                let surface: &dyn SurfaceEval = &surface_on_face.geometry;
                let (pc_t0, pc_t1) = pcurve.parameter_range();
                let mut max_dev = 0.0_f64;

                for k in 0..PC_SAMPLES {
                    let frac = (k as f64 + 0.5) / PC_SAMPLES as f64;
                    let pc_t = pc_t0 + (pc_t1 - pc_t0) * frac;
                    let uv = pcurve.point_at(pc_t);
                    let pt_surface = surface.point_at(uv.x, uv.y);

                    // Project onto edge curve to find closest distance.
                    if let Some(t_closest) =
                        closest_point_on_curve(&curve_on_edge.geometry, &pt_surface)
                    {
                        let pt_curve = curve_on_edge.geometry.point_at(t_closest);
                        let dist = pt_surface.distance_to(&pt_curve);
                        max_dev = max_dev.max(dist);
                    }
                }

                if max_dev > tolerance.linear() {
                    report.push(DefectRecord {
                        category: DefectCategory::Geometric,
                        severity: Severity::Warning,
                        description: format!(
                            "Trim curve deviates from edge curve by \
                             {:.6} (tolerance: {:.6})",
                            max_dev,
                            tolerance.linear()
                        ),
                        entities: vec![
                            EntityRef::Face(face_handle),
                            EntityRef::HalfEdge(cursor),
                            EntityRef::Edge(he.edge),
                        ],
                        location: None,
                        measured_value: Some(max_dev),
                        threshold: Some(tolerance.linear()),
                        remediation: Some(Remediation {
                            action: RemediationAction::RecomputePCurve,
                            parameter: None,
                            value: None,
                            rationale: "Recompute pcurve from 3D edge \
                                        curve projection onto face surface"
                                .into(),
                        }),
                    });
                }
            }

            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::healing::options::SampleDensity;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn surface_containment_on_valid_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut report = HealingReport::new();
        check_surface_containment(&brep, &mut report);
        // Cube faces use planes with infinite domains, so check is skipped.
        assert!(
            report.is_empty(),
            "Valid cube should have no surface containment defects"
        );
    }

    #[test]
    fn surface_containment_on_empty_brep() {
        let brep = BRep::new(Tolerance::default());
        let mut report = HealingReport::new();
        check_surface_containment(&brep, &mut report);
        assert!(report.is_empty());
    }

    #[test]
    fn boundary_containment_on_valid_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let density = SampleDensity::Standard;
        let mut report = HealingReport::new();
        check_boundary_containment(&brep, &tol, &density, &mut report);
        assert!(
            report.is_empty(),
            "Valid cube should have no boundary containment defects, got {}",
            report.len()
        );
    }

    #[test]
    fn boundary_containment_on_empty_brep() {
        let brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let density = SampleDensity::Standard;
        let mut report = HealingReport::new();
        check_boundary_containment(&brep, &tol, &density, &mut report);
        assert!(report.is_empty());
    }

    #[test]
    fn surface_containment_on_sphere_face_with_pcurves() {
        use crate::geometry::Point3;
        use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let fh = brep.faces.keys().next().unwrap();
        // Sphere has a finite parameter range, so the check won't skip it.
        brep.faces[fh].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.5, 0.5, 0.5),
                radius: 1.0,
            },
            orientation: true,
        });
        let mut report = HealingReport::new();
        check_surface_containment(&brep, &mut report);
        // Cube pcurves are designed for planes; on a sphere they may lie
        // outside the domain. Either way, the check exercises the main path.
        // We just verify it doesn't panic.
        let _ = report.len();
    }

    #[test]
    fn boundary_containment_detects_deviation_on_corrupted_pcurve() {
        use crate::geometry::Point3;
        use crate::geometry::nurbs_curve::NurbsCurve;

        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // Corrupt a pcurve so it maps to a wildly wrong 3D point.
        let he_h = brep.half_edges.keys().next().unwrap();
        // Replace pcurve with one that maps to (100, 100, 0) in UV space.
        brep.half_edges[he_h].pcurve = Some(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(100.0, 100.0, 0.0),
                Point3::new(200.0, 200.0, 0.0),
            ],
            weights: vec![1.0, 1.0],
        });
        let tol = Tolerance::default();
        let density = SampleDensity::Standard;
        let mut report = HealingReport::new();
        check_boundary_containment(&brep, &tol, &density, &mut report);
        // The badly-mapped pcurve will produce 3D points far from the edge curve.
        assert!(
            !report.is_empty(),
            "Should detect trim curve deviation from corrupted pcurve"
        );
    }

    #[test]
    fn surface_containment_pcurve_outside_domain() {
        use crate::geometry::Point3;
        use crate::geometry::nurbs_curve::NurbsCurve;
        use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let fh = brep.faces.keys().next().unwrap();
        // Sphere has finite parameter domain [0, 2PI] x [0, PI].
        brep.faces[fh].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
            orientation: true,
        });
        // Set a pcurve on the first half-edge of this face's outer loop
        // that lies way outside the sphere's parameter domain.
        let face_data = &brep.faces[fh];
        let loop_data = &brep.loops[face_data.outer_loop];
        let he_h = loop_data.half_edge;
        brep.half_edges[he_h].pcurve = Some(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(50.0, 50.0, 0.0), Point3::new(60.0, 60.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let mut report = HealingReport::new();
        check_surface_containment(&brep, &mut report);
        assert!(
            !report.is_empty(),
            "Should detect pcurve outside sphere parameter domain"
        );
    }
}
