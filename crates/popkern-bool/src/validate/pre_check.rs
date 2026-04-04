//! Pre-boolean validation of input solids.
//!
//! Validates both input bodies before the boolean operation begins.
//! Catches problems early that would cause failures deep in the pipeline.

use popkern_core::error::ValidationError;
use popkern_core::topology::BRep;

use crate::error::Severity;

/// A single finding from pre-check validation.
#[derive(Debug)]
pub struct PreCheckFinding {
    /// Severity of the finding.
    pub severity: Severity,
    /// Which operand has the issue (0 = A, 1 = B).
    pub operand: usize,
    /// Human-readable description.
    pub message: String,
    /// The underlying validation error, if any.
    pub error: Option<ValidationError>,
}

/// Result of pre-boolean validation.
#[derive(Debug, Default)]
pub struct PreCheckReport {
    /// All findings.
    pub findings: Vec<PreCheckFinding>,
}

impl PreCheckReport {
    /// Returns true if no errors were found (warnings are ok).
    pub fn is_ok(&self) -> bool {
        !self
            .findings
            .iter()
            .any(|f| matches!(f.severity, Severity::Error | Severity::Critical))
    }

    /// Number of errors (not warnings).
    pub fn error_count(&self) -> usize {
        self.findings
            .iter()
            .filter(|f| matches!(f.severity, Severity::Error | Severity::Critical))
            .count()
    }

    /// Number of warnings.
    pub fn warning_count(&self) -> usize {
        self.findings
            .iter()
            .filter(|f| matches!(f.severity, Severity::Warning))
            .count()
    }
}

/// Validate both input bodies before a boolean operation.
///
/// Checks:
/// - Euler-Poincaré formula (closed shell)
/// - Half-edge pairing and cycles
/// - Manifold edges (each shared by exactly 2 faces)
/// - Face orientation consistency
/// - Loop closure (no gaps between edges)
/// - Missing surface/curve geometry
///
/// Optionally checks self-intersection (expensive).
pub fn pre_check(a: &BRep, b: &BRep, check_self_intersection: bool) -> PreCheckReport {
    let mut report = PreCheckReport::default();

    let tol_linear = a.tolerance.linear;
    validate_one_body(a, 0, check_self_intersection, tol_linear, &mut report);
    validate_one_body(b, 1, check_self_intersection, tol_linear, &mut report);

    // Check that both bodies have faces
    if a.faces.is_empty() {
        report.findings.push(PreCheckFinding {
            severity: Severity::Error,
            operand: 0,
            message: "Body A has no faces".to_string(),
            error: None,
        });
    }
    if b.faces.is_empty() {
        report.findings.push(PreCheckFinding {
            severity: Severity::Error,
            operand: 1,
            message: "Body B has no faces".to_string(),
            error: None,
        });
    }

    // Check that faces have surface geometry
    let a_missing_geo = a.faces.values().filter(|f| f.surface.is_none()).count();
    if a_missing_geo > 0 {
        report.findings.push(PreCheckFinding {
            severity: Severity::Warning,
            operand: 0,
            message: format!("Body A: {a_missing_geo} face(s) without surface geometry"),
            error: None,
        });
    }
    let b_missing_geo = b.faces.values().filter(|f| f.surface.is_none()).count();
    if b_missing_geo > 0 {
        report.findings.push(PreCheckFinding {
            severity: Severity::Warning,
            operand: 1,
            message: format!("Body B: {b_missing_geo} face(s) without surface geometry"),
            error: None,
        });
    }

    report
}

/// Run core validation checks on a single body.
fn validate_one_body(
    brep: &BRep,
    operand: usize,
    check_self_intersection: bool,
    tol_linear: f64,
    report: &mut PreCheckReport,
) {
    let label = if operand == 0 { "A" } else { "B" };

    // Structural checks — errors prevent boolean from producing valid output
    let structural_checks: Vec<(&str, Result<(), ValidationError>)> = vec![
        ("Euler-Poincaré", brep.validate_euler_poincare()),
        ("half-edge pairing", brep.validate_half_edge_pairing()),
        ("half-edge cycles", brep.validate_half_edge_cycles()),
        ("vertex consistency", brep.validate_vertex_consistency()),
        ("manifold edges", brep.validate_manifold_edges()),
    ];

    for (name, result) in structural_checks {
        if let Err(e) = result {
            report.findings.push(PreCheckFinding {
                severity: Severity::Error,
                operand,
                message: format!("Body {label}: {name} failed — {e}"),
                error: Some(e),
            });
        }
    }

    // Geometric checks — warnings (may still produce usable output)
    let geometric_checks: Vec<(&str, Result<(), ValidationError>)> = vec![
        ("face orientation", brep.validate_face_orientation()),
        ("loop closure", brep.validate_loop_closure()),
    ];

    for (name, result) in geometric_checks {
        if let Err(e) = result {
            report.findings.push(PreCheckFinding {
                severity: Severity::Warning,
                operand,
                message: format!("Body {label}: {name} — {e}"),
                error: Some(e),
            });
        }
    }

    // Edge validity check
    let edge_report =
        crate::imprint::edge_validity::validate_edges(brep, tol_linear * 0.5, tol_linear * 10.0);
    if !edge_report.is_valid() {
        let msg = format!(
            "Body {label}: {} edge issues ({} collapsed, {} short, {} endpoint mismatch)",
            edge_report.collapsed_range.len()
                + edge_report.too_short.len()
                + edge_report.endpoint_mismatch.len(),
            edge_report.collapsed_range.len(),
            edge_report.too_short.len(),
            edge_report.endpoint_mismatch.len()
        );
        report.findings.push(PreCheckFinding {
            severity: Severity::Warning,
            operand,
            message: msg,
            error: None,
        });
    }

    // Optional expensive check
    if check_self_intersection && let Err(e) = brep.validate_self_intersection() {
        report.findings.push(PreCheckFinding {
            severity: Severity::Error,
            operand,
            message: format!("Body {label}: self-intersection detected — {e}"),
            error: Some(e),
        });
    }
}
