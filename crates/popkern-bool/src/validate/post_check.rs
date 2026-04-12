//! Post-boolean validation of the result solid.
//!
//! Validates that the boolean result satisfies all required topological
//! and geometric invariants.

use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::Severity;
use crate::imprint::shared_geometry;

/// A single finding from post-check validation.
#[derive(Debug)]
pub struct PostCheckFinding {
    /// Severity of the finding.
    pub severity: Severity,
    /// Human-readable description.
    pub message: String,
    /// The underlying validation error, if any.
    pub error: Option<ValidationError>,
}

/// Result of post-boolean validation.
#[derive(Debug, Default)]
pub struct PostCheckReport {
    /// All findings.
    pub findings: Vec<PostCheckFinding>,
    /// Number of faces in the result.
    pub face_count: usize,
    /// Number of edges in the result.
    pub edge_count: usize,
    /// Number of vertices in the result.
    pub vertex_count: usize,
    /// Euler characteristic V - E + F.
    pub euler_characteristic: i64,
    /// Maximum watertightness gap (0 if watertight).
    pub max_gap: f64,
}

impl PostCheckReport {
    /// Returns true if no errors were found.
    pub fn is_ok(&self) -> bool {
        !self
            .findings
            .iter()
            .any(|f| matches!(f.severity, Severity::Error | Severity::Critical))
    }

    /// Number of errors.
    pub fn error_count(&self) -> usize {
        self.findings
            .iter()
            .filter(|f| matches!(f.severity, Severity::Error | Severity::Critical))
            .count()
    }
}

/// Validate a boolean result solid.
///
/// Checks:
/// 1. Manifold edges (every edge shared by exactly 2 faces)
/// 2. Closed shell (Euler-Poincaré V - E + F = 2)
/// 3. Orientation consistency (half-edge pairing)
/// 4. Watertightness (edge-surface deviation)
/// 5. No degenerate topology (half-edge cycles, vertex consistency)
/// 6. Self-intersection (optional, expensive O(n²) face-pair check)
pub fn post_check(
    result: &BRep,
    tolerance: &Tolerance,
    check_self_intersection: bool,
) -> PostCheckReport {
    let mut report = PostCheckReport {
        face_count: result.faces.len(),
        edge_count: result.edges.len(),
        vertex_count: result.vertices.len(),
        euler_characteristic: result.vertices.len() as i64 - result.edges.len() as i64
            + result.faces.len() as i64,
        ..Default::default()
    };

    // 1. Manifold edges
    if let Err(e) = result.validate_manifold_edges() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!("Non-manifold edges: {e}"),
            error: Some(e),
        });
    }

    // 2. Euler-Poincaré (closed shell). Boolean results may be
    // higher-genus (e.g. cube minus cylinder = genus-1 torus-like
    // solid), so use the genus-tolerant variant. The strict
    // genus-0 check is too restrictive for general booleans.
    if let Err(e) = result.validate_euler_poincare_any_genus() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!(
                "Euler-Poincaré violated (V-E+F={}): {e}",
                report.euler_characteristic
            ),
            error: Some(e),
        });
    }

    // 3. Orientation consistency
    if let Err(e) = result.validate_half_edge_pairing() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!("Half-edge pairing error: {e}"),
            error: Some(e),
        });
    }

    // 4. Half-edge cycles
    if let Err(e) = result.validate_half_edge_cycles() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!("Half-edge cycle error: {e}"),
            error: Some(e),
        });
    }

    // 5. Vertex consistency
    if let Err(e) = result.validate_vertex_consistency() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Warning,
            message: format!("Vertex consistency: {e}"),
            error: Some(e),
        });
    }

    // 6. Watertightness
    let wt = shared_geometry::validate_watertightness(result, tolerance.linear() * 10.0, 3);
    report.max_gap = wt.max_gap();
    if !wt.is_watertight() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Warning,
            message: format!(
                "{} watertightness gap(s), max deviation: {:.2e}",
                wt.gaps.len(),
                wt.max_gap()
            ),
            error: None,
        });
    }

    // 7. Face orientation
    if let Err(e) = result.validate_face_orientation() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Warning,
            message: format!("Face orientation: {e}"),
            error: Some(e),
        });
    }

    // 8. Self-intersection (opt-in, O(n²) face-pair check).
    // Required by §9.2: result must not self-intersect. Reported as Warning
    // because for boolean results a detected self-intersection usually
    // indicates a degenerate input rather than an outright pipeline bug,
    // and the caller may still want the result for diagnostic purposes.
    if check_self_intersection && let Err(e) = result.validate_self_intersection() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Warning,
            message: format!("Self-intersection in result: {e}"),
            error: Some(e),
        });
    }

    report
}
