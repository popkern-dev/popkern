//! Pre-boolean validation of input solids.
//!
//! Validates both input bodies before the boolean operation begins.
//! Catches problems early that would cause failures deep in the pipeline.

use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

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

/// A face flagged as geometrically degenerate by pre-check.
#[derive(Debug, Clone)]
pub struct DegenerateFace {
    /// Which operand the face is on (0 = A, 1 = B).
    pub operand: usize,
    /// The face handle.
    pub face: FaceHandle,
    /// Approximate area as computed by `analyze_face_polygon`.
    pub area: f64,
    /// Why it was flagged (`"zero_area"`, `"too_few_vertices"`, etc.).
    pub reason: String,
}

/// A potentially problematic face pair flagged by pre-check (§12.4 [BEYOND]).
///
/// Identified by an AABB-overlap test followed by either a near-coincidence
/// or near-tangent classification. Helps callers preemptively repair or
/// adjust geometry before running the full boolean.
#[derive(Debug, Clone)]
pub struct ProblematicPair {
    /// Face from operand A.
    pub face_a: FaceHandle,
    /// Face from operand B.
    pub face_b: FaceHandle,
    /// Why this pair is potentially problematic.
    pub kind: ProblematicKind,
    /// Quantitative measure: distance for `NearCoincident`, sin(angle) for `NearTangent`.
    pub measure: f64,
}

/// Classification of a problematic face pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProblematicKind {
    /// Faces are within 10× tolerance and same orientation: likely
    /// to trip the coincidence detector at a tolerance cliff edge.
    NearCoincident,
    /// Surface normals are nearly parallel (within angular tolerance):
    /// the SSI solver's Jacobian is rank-deficient at this configuration.
    NearTangent,
}

/// Result of pre-boolean validation.
#[derive(Debug, Default)]
pub struct PreCheckReport {
    /// All findings.
    pub findings: Vec<PreCheckFinding>,
    /// Faces flagged as geometrically degenerate (zero area, too few vertices).
    pub degenerate_faces: Vec<DegenerateFace>,
    /// Face pairs flagged as potentially problematic (§12.4 [BEYOND]).
    pub problematic_pairs: Vec<ProblematicPair>,
    /// Heuristic estimate of how likely the boolean is to succeed (§12.4 [BEYOND]).
    /// `1.0` = no issues detected; `0.0` = many issues. Computed by
    /// [`PreCheckReport::compute_success_likelihood`] from the other fields.
    pub success_likelihood: f64,
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
/// - Degenerate faces (§12.4)
/// - Problematic face pairs and success likelihood (§12.4 [BEYOND])
///
/// Optionally checks self-intersection (expensive).
pub fn pre_check(a: &BRep, b: &BRep, check_self_intersection: bool) -> PreCheckReport {
    let mut report = PreCheckReport::default();

    let tol_linear = a.tolerance.linear();
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

    // Detect degenerate faces (§12.4) — area below tolerance² is a
    // strong indicator that the face will collapse during subdivision.
    let tol = Tolerance::with_linear(Tolerance::default(), tol_linear);
    detect_degenerate_faces(a, 0, &tol, &mut report);
    detect_degenerate_faces(b, 1, &tol, &mut report);

    // §12.4: G0/G1 continuity check at shared edges. G0 = positions
    // agree across the edge; G1 = surface tangents agree.
    super::pre_check_geometry::detect_g0_g1_violations(a, 0, &tol, &mut report);
    super::pre_check_geometry::detect_g0_g1_violations(b, 1, &tol, &mut report);

    // §12.4: detect surface types that the boolean engine handles
    // less robustly (NURBS-NURBS, torus). Reported as warnings — the
    // boolean MAY succeed, but the caller should be aware.
    detect_unsupported_geometry(a, b, &mut report);

    // Identify problematic face pairs (§12.4 [BEYOND]).
    super::pre_check_geometry::detect_problematic_pairs(a, b, &tol, &mut report);

    // Compute success likelihood from accumulated findings.
    report.success_likelihood = report.compute_success_likelihood();

    report
}

/// Warn when either body uses surface types known to be less robust
/// in the boolean engine (NURBS-NURBS pairs, torus surfaces).
fn detect_unsupported_geometry(a: &BRep, b: &BRep, report: &mut PreCheckReport) {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let count = |brep: &BRep, predicate: fn(&SurfaceGeometry) -> bool| -> usize {
        brep.faces
            .values()
            .filter_map(|f| f.surface.as_ref())
            .filter(|s| predicate(&s.geometry))
            .count()
    };

    let a_nurbs = count(a, |s| matches!(s, SurfaceGeometry::Nurbs(_)));
    let b_nurbs = count(b, |s| matches!(s, SurfaceGeometry::Nurbs(_)));
    let a_torus = count(a, |s| matches!(s, SurfaceGeometry::Torus { .. }));
    let b_torus = count(b, |s| matches!(s, SurfaceGeometry::Torus { .. }));

    if a_nurbs > 0 && b_nurbs > 0 {
        report.findings.push(PreCheckFinding {
            severity: Severity::Warning,
            operand: 0,
            message: format!(
                "Both bodies contain NURBS faces (A: {a_nurbs}, B: {b_nurbs}). NURBS-NURBS intersection is less robust than analytic-pair intersection."
            ),
            error: None,
        });
    }
    if a_torus > 0 || b_torus > 0 {
        report.findings.push(PreCheckFinding {
            severity: Severity::Warning,
            operand: 0,
            message: format!(
                "Torus surfaces present (A: {a_torus}, B: {b_torus}). Torus intersections are tested but less stable than other analytic surfaces."
            ),
            error: None,
        });
    }
}

/// Walk every face of `brep` and flag those whose polygon area falls
/// below `tolerance.linear()²` or that have fewer than 3 vertices.
fn detect_degenerate_faces(
    brep: &BRep,
    operand: usize,
    tolerance: &Tolerance,
    report: &mut PreCheckReport,
) {
    use crate::clip::face_vertex_positions;
    use crate::imprint::sliver::{SliverThresholds, analyze_face_polygon};

    let thresholds = SliverThresholds::from_tolerance(tolerance.linear());
    for fh in brep.faces.keys() {
        let positions = face_vertex_positions(brep, fh);
        if positions.len() < 3 {
            report.degenerate_faces.push(DegenerateFace {
                operand,
                face: fh,
                area: 0.0,
                reason: "too_few_vertices".into(),
            });
            continue;
        }
        let analysis = analyze_face_polygon(&positions, &thresholds);
        if analysis.area < thresholds.min_area {
            report.degenerate_faces.push(DegenerateFace {
                operand,
                face: fh,
                area: analysis.area,
                reason: format!(
                    "area {:.2e} < tolerance² {:.2e}",
                    analysis.area, thresholds.min_area
                ),
            });
        }
    }
}

impl PreCheckReport {
    /// Compute the success likelihood (0.0–1.0) from the accumulated
    /// findings, degenerate faces, and problematic pairs.
    ///
    /// Heuristic, not a proof: high values mean the boolean is more
    /// likely to succeed; low values mean significant manual repair
    /// is probably needed before running the operation.
    pub fn compute_success_likelihood(&self) -> f64 {
        let mut score = 1.0_f64;
        for f in &self.findings {
            score -= match f.severity {
                Severity::Critical => 0.50,
                Severity::Error => 0.30,
                Severity::Warning => 0.05,
                Severity::Info => 0.0,
            };
        }
        score -= self.degenerate_faces.len() as f64 * 0.05;
        for p in &self.problematic_pairs {
            score -= match p.kind {
                ProblematicKind::NearCoincident => 0.10,
                ProblematicKind::NearTangent => 0.05,
            };
        }
        score.clamp(0.0, 1.0)
    }
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
