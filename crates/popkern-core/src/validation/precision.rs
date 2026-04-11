//! Precision validation checks for BRep models.
//!
//! Two entry points:
//! - [`BRep::validate_precision_fast`] — returns on first error (used by `validate()`).
//! - [`BRep::validate_precision_thorough`] — collects all errors (used by `validate_all()`).
//!
//! Current checks:
//! 1. **Bounds check:** entity tolerance must be in `[numerical_floor, max_imprecision]`.
//! 2. **Containment hierarchy:** vertex tolerance must not exceed the tolerance
//!    of any edge that contains it.

use crate::error::ValidationError;
use crate::precision::{
    DiagnosticSeverity, EntityHandle, PrecisionDiagnostic, PrecisionDiagnosticKind, Remediation,
};
use crate::topology::BRep;

impl BRep {
    /// Fast precision validation — returns on first error.
    ///
    /// Checks all entities with explicit tolerances against the BRep's
    /// `max_imprecision` upper bound and `tolerance.numerical_floor()` lower
    /// bound, then verifies the containment hierarchy (vertex tolerance must
    /// not exceed edge tolerance for any incident edge).
    ///
    /// This is a no-op on models where all entity tolerances are `None`.
    pub fn validate_precision_fast(&self) -> Result<(), ValidationError> {
        let floor = self.tolerance.numerical_floor();
        let max = self.max_imprecision;

        // --- Bounds checks on vertices ---
        for (vh, v) in &self.vertices {
            if let Some(tol) = v.tolerance {
                self.check_bounds(tol, max, floor, EntityHandle::Vertex(vh))?;
            }
        }

        // --- Bounds checks on edges ---
        for (eh, e) in &self.edges {
            if let Some(tol) = e.tolerance {
                self.check_bounds(tol, max, floor, EntityHandle::Edge(eh))?;
            }
        }

        // --- Bounds checks on faces ---
        for (fh, f) in &self.faces {
            if let Some(tol) = f.tolerance {
                self.check_bounds(tol, max, floor, EntityHandle::Face(fh))?;
            }
        }

        // --- Containment hierarchy: vertex tolerance <= edge tolerance ---
        for (eh, edge) in &self.edges {
            let he = match self.half_edges.get(edge.half_edge) {
                Some(he) => he,
                None => continue,
            };
            let twin_he = match self.half_edges.get(he.twin) {
                Some(t) => t,
                None => continue,
            };

            let edge_tol = self.edge_tolerance(eh);
            let v1 = he.origin;
            let v2 = twin_he.origin;

            if self.vertex_tolerance(v1) > edge_tol {
                return Err(ValidationError::Precision {
                    diagnostic: PrecisionDiagnostic {
                        severity: DiagnosticSeverity::Warning,
                        kind: PrecisionDiagnosticKind::ContainmentViolation,
                        entities: vec![EntityHandle::Vertex(v1), EntityHandle::Edge(eh)],
                        measured: self.vertex_tolerance(v1),
                        threshold: edge_tol,
                        remediation: Remediation::RecoverVertexPrecision(v1),
                    },
                });
            }

            if self.vertex_tolerance(v2) > edge_tol {
                return Err(ValidationError::Precision {
                    diagnostic: PrecisionDiagnostic {
                        severity: DiagnosticSeverity::Warning,
                        kind: PrecisionDiagnosticKind::ContainmentViolation,
                        entities: vec![EntityHandle::Vertex(v2), EntityHandle::Edge(eh)],
                        measured: self.vertex_tolerance(v2),
                        threshold: edge_tol,
                        remediation: Remediation::RecoverVertexPrecision(v2),
                    },
                });
            }
        }

        Ok(())
    }

    /// Thorough precision validation — collects all errors without returning early.
    ///
    /// Runs the same checks as [`validate_precision_fast`] but accumulates
    /// every violation into the returned `Vec` instead of short-circuiting.
    pub fn validate_precision_thorough(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        let floor = self.tolerance.numerical_floor();
        let max = self.max_imprecision;

        // --- Bounds checks on vertices ---
        for (vh, v) in &self.vertices {
            if let Some(tol) = v.tolerance
                && let Err(e) = self.check_bounds(tol, max, floor, EntityHandle::Vertex(vh))
            {
                errors.push(e);
            }
        }

        // --- Bounds checks on edges ---
        for (eh, e) in &self.edges {
            if let Some(tol) = e.tolerance
                && let Err(e) = self.check_bounds(tol, max, floor, EntityHandle::Edge(eh))
            {
                errors.push(e);
            }
        }

        // --- Bounds checks on faces ---
        for (fh, f) in &self.faces {
            if let Some(tol) = f.tolerance
                && let Err(e) = self.check_bounds(tol, max, floor, EntityHandle::Face(fh))
            {
                errors.push(e);
            }
        }

        // --- Containment hierarchy: vertex tolerance <= edge tolerance ---
        for (eh, edge) in &self.edges {
            let he = match self.half_edges.get(edge.half_edge) {
                Some(he) => he,
                None => continue,
            };
            let twin_he = match self.half_edges.get(he.twin) {
                Some(t) => t,
                None => continue,
            };

            let edge_tol = self.edge_tolerance(eh);
            let v1 = he.origin;
            let v2 = twin_he.origin;

            if self.vertex_tolerance(v1) > edge_tol {
                errors.push(ValidationError::Precision {
                    diagnostic: PrecisionDiagnostic {
                        severity: DiagnosticSeverity::Warning,
                        kind: PrecisionDiagnosticKind::ContainmentViolation,
                        entities: vec![EntityHandle::Vertex(v1), EntityHandle::Edge(eh)],
                        measured: self.vertex_tolerance(v1),
                        threshold: edge_tol,
                        remediation: Remediation::RecoverVertexPrecision(v1),
                    },
                });
            }

            if self.vertex_tolerance(v2) > edge_tol {
                errors.push(ValidationError::Precision {
                    diagnostic: PrecisionDiagnostic {
                        severity: DiagnosticSeverity::Warning,
                        kind: PrecisionDiagnosticKind::ContainmentViolation,
                        entities: vec![EntityHandle::Vertex(v2), EntityHandle::Edge(eh)],
                        measured: self.vertex_tolerance(v2),
                        threshold: edge_tol,
                        remediation: Remediation::RecoverVertexPrecision(v2),
                    },
                });
            }
        }

        errors
    }

    /// Check that `tol` is within `[floor, max]`, returning the appropriate
    /// `ValidationError::Precision` on violation.
    fn check_bounds(
        &self,
        tol: f64,
        max: f64,
        floor: f64,
        entity: EntityHandle,
    ) -> Result<(), ValidationError> {
        if tol > max {
            return Err(ValidationError::Precision {
                diagnostic: PrecisionDiagnostic {
                    severity: DiagnosticSeverity::Error,
                    kind: PrecisionDiagnosticKind::ExceedsMaxBound,
                    entities: vec![entity],
                    measured: tol,
                    threshold: max,
                    remediation: match entity {
                        EntityHandle::Vertex(vh) => Remediation::RecoverVertexPrecision(vh),
                        EntityHandle::Edge(eh) => Remediation::RecoverEdgePrecision(eh),
                        EntityHandle::Face(_) => Remediation::Manual("re-fit face surface".into()),
                    },
                },
            });
        }
        if tol < floor {
            return Err(ValidationError::Precision {
                diagnostic: PrecisionDiagnostic {
                    severity: DiagnosticSeverity::Error,
                    kind: PrecisionDiagnosticKind::BelowNumericalFloor,
                    entities: vec![entity],
                    measured: tol,
                    threshold: floor,
                    remediation: Remediation::Manual(format!(
                        "tolerance {tol} is below numerical floor {floor}; \
                         set to at least {floor}"
                    )),
                },
            });
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::Point3;
    use crate::precision::{DiagnosticSeverity, EntityHandle, PrecisionDiagnosticKind};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// Helper: build a minimal BRep with one edge (two vertices, two half-edges).
    fn brep_with_one_edge() -> BRep {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let _mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        brep
    }

    // -----------------------------------------------------------------------
    // No-op on clean models
    // -----------------------------------------------------------------------

    #[test]
    fn clean_model_passes_fast() {
        let brep = brep_with_one_edge();
        assert!(brep.validate_precision_fast().is_ok());
    }

    #[test]
    fn clean_model_passes_thorough() {
        let brep = brep_with_one_edge();
        assert!(brep.validate_precision_thorough().is_empty());
    }

    #[test]
    fn unit_cube_passes_precision() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        assert!(cube.validate_precision_fast().is_ok());
        assert!(cube.validate_precision_thorough().is_empty());
    }

    // -----------------------------------------------------------------------
    // ExceedsMaxBound
    // -----------------------------------------------------------------------

    #[test]
    fn vertex_exceeds_max_bound_fast() {
        let mut brep = brep_with_one_edge();
        let vh = brep.vertices.keys().next().unwrap();
        brep.vertices[vh].tolerance = Some(brep.max_imprecision + 1.0);

        let err = brep.validate_precision_fast().unwrap_err();
        if let crate::error::ValidationError::Precision { diagnostic } = &err {
            assert_eq!(diagnostic.severity, DiagnosticSeverity::Error);
            assert!(matches!(
                diagnostic.kind,
                PrecisionDiagnosticKind::ExceedsMaxBound
            ));
        } else {
            panic!("expected Precision error, got: {err:?}");
        }
    }

    #[test]
    fn edge_exceeds_max_bound_thorough() {
        let mut brep = brep_with_one_edge();
        let eh = brep.edges.keys().next().unwrap();
        brep.edges[eh].tolerance = Some(brep.max_imprecision + 1.0);

        let errors = brep.validate_precision_thorough();
        assert!(!errors.is_empty());
        assert!(errors.iter().any(|e| matches!(
            e,
            crate::error::ValidationError::Precision {
                diagnostic
            } if matches!(diagnostic.kind, PrecisionDiagnosticKind::ExceedsMaxBound)
        )));
    }

    #[test]
    fn face_exceeds_max_bound() {
        let mut brep = brep_with_one_edge();
        let fh = brep.faces.keys().next().unwrap();
        brep.faces[fh].tolerance = Some(brep.max_imprecision + 1.0);

        let err = brep.validate_precision_fast().unwrap_err();
        if let crate::error::ValidationError::Precision { diagnostic } = &err {
            assert!(matches!(
                diagnostic.kind,
                PrecisionDiagnosticKind::ExceedsMaxBound
            ));
        } else {
            panic!("expected Precision error");
        }
    }

    // -----------------------------------------------------------------------
    // BelowNumericalFloor
    // -----------------------------------------------------------------------

    #[test]
    fn vertex_below_floor_fast() {
        let mut brep = brep_with_one_edge();
        let vh = brep.vertices.keys().next().unwrap();
        let floor = brep.tolerance.numerical_floor();
        brep.vertices[vh].tolerance = Some(floor * 0.1);

        let err = brep.validate_precision_fast().unwrap_err();
        if let crate::error::ValidationError::Precision { diagnostic } = &err {
            assert!(matches!(
                diagnostic.kind,
                PrecisionDiagnosticKind::BelowNumericalFloor
            ));
        } else {
            panic!("expected Precision error");
        }
    }

    #[test]
    fn edge_below_floor_thorough() {
        let mut brep = brep_with_one_edge();
        let eh = brep.edges.keys().next().unwrap();
        let floor = brep.tolerance.numerical_floor();
        brep.edges[eh].tolerance = Some(floor * 0.01);

        let errors = brep.validate_precision_thorough();
        assert!(errors.iter().any(|e| matches!(
            e,
            crate::error::ValidationError::Precision {
                diagnostic
            } if matches!(diagnostic.kind, PrecisionDiagnosticKind::BelowNumericalFloor)
        )));
    }

    // -----------------------------------------------------------------------
    // ContainmentViolation
    // -----------------------------------------------------------------------

    #[test]
    fn containment_violation_fast() {
        let mut brep = brep_with_one_edge();
        // Set a vertex tolerance higher than the edge tolerance.
        let vh = brep.vertices.keys().next().unwrap();
        let eh = brep.edges.keys().next().unwrap();
        let edge_tol = brep.tolerance.linear() * 2.0;
        let vertex_tol = brep.tolerance.linear() * 3.0;
        brep.edges[eh].tolerance = Some(edge_tol);
        brep.vertices[vh].tolerance = Some(vertex_tol);

        let err = brep.validate_precision_fast().unwrap_err();
        if let crate::error::ValidationError::Precision { diagnostic } = &err {
            assert!(matches!(
                diagnostic.kind,
                PrecisionDiagnosticKind::ContainmentViolation
            ));
            assert_eq!(diagnostic.entities.len(), 2);
            assert!(
                diagnostic
                    .entities
                    .iter()
                    .any(|e| matches!(e, EntityHandle::Vertex(_)))
            );
            assert!(
                diagnostic
                    .entities
                    .iter()
                    .any(|e| matches!(e, EntityHandle::Edge(_)))
            );
        } else {
            panic!("expected Precision error");
        }
    }

    #[test]
    fn containment_violation_thorough_collects_all() {
        let mut brep = brep_with_one_edge();
        // Set both vertices above the edge tolerance.
        let eh = brep.edges.keys().next().unwrap();
        let edge_tol = brep.tolerance.linear() * 2.0;
        brep.edges[eh].tolerance = Some(edge_tol);

        for (_, v) in brep.vertices.iter_mut() {
            v.tolerance = Some(edge_tol + 1.0);
        }

        let errors = brep.validate_precision_thorough();
        let containment_errors: Vec<_> = errors
            .iter()
            .filter(|e| {
                matches!(
                    e,
                    crate::error::ValidationError::Precision {
                        diagnostic
                    } if matches!(diagnostic.kind, PrecisionDiagnosticKind::ContainmentViolation)
                )
            })
            .collect();
        // Should detect at least one containment violation per vertex.
        assert!(
            containment_errors.len() >= 2,
            "expected >= 2 containment violations, got {}",
            containment_errors.len()
        );
    }

    // -----------------------------------------------------------------------
    // No false positive: valid explicit tolerances
    // -----------------------------------------------------------------------

    #[test]
    fn valid_explicit_tolerances_pass() {
        let mut brep = brep_with_one_edge();
        let floor = brep.tolerance.numerical_floor();
        let max = brep.max_imprecision;
        let safe = (floor + max) / 2.0;

        // Set all entity tolerances to a safe mid-range value.
        for (_, v) in brep.vertices.iter_mut() {
            v.tolerance = Some(safe);
        }
        for (_, e) in brep.edges.iter_mut() {
            // Edge tolerance must be >= vertex tolerance for containment.
            e.tolerance = Some(safe);
        }
        for (_, f) in brep.faces.iter_mut() {
            f.tolerance = Some(safe);
        }

        assert!(brep.validate_precision_fast().is_ok());
        assert!(brep.validate_precision_thorough().is_empty());
    }

    // -----------------------------------------------------------------------
    // Multiple errors collected by thorough
    // -----------------------------------------------------------------------

    #[test]
    fn thorough_collects_multiple_error_types() {
        let mut brep = brep_with_one_edge();
        let floor = brep.tolerance.numerical_floor();

        // One vertex below floor, one face above max.
        let vh = brep.vertices.keys().next().unwrap();
        brep.vertices[vh].tolerance = Some(floor * 0.01);

        let fh = brep.faces.keys().next().unwrap();
        brep.faces[fh].tolerance = Some(brep.max_imprecision + 10.0);

        let errors = brep.validate_precision_thorough();
        assert!(
            errors.len() >= 2,
            "expected >= 2 errors, got {}",
            errors.len()
        );
    }
}
