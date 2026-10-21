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
    /// Measured geometric value (§13.1) — e.g. the maximum
    /// watertightness gap in linear units, or `None` if the finding
    /// is intrinsically topological.
    pub measured_value: Option<f64>,
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
    post_check_scoped(
        result,
        tolerance,
        check_self_intersection,
        true,
        &Default::default(),
    )
}

/// [`post_check`]'s GATES only: every Error/Critical check, none of the
/// Warning-only diagnostics (vertex consistency, the watertightness
/// gap measurement, face orientation, self-intersection). The serve
/// path reads only the gates; the watertightness diagnostic alone
/// sampled every edge of the result per boolean, so a body carrying
/// earlier tapped holes paid for re-measuring all of them on every
/// later boolean (PK-731: +43 ms per earlier tap, per boolean).
/// The gating findings are exactly [`post_check`]'s, except that
/// `carried` edges — built from unsplit operand edges, so carrying the
/// operand's own curve and face surfaces — skip the curve-surface
/// check (check 9): it re-measured every far edge of the body on every
/// boolean (PK-802). Serving boundaries (`validate_all`) still measure
/// every edge. A caller that knows the body passes `validate()`
/// (`structure_validated`) skips the gates `validate()` contains —
/// manifold edges, Euler, loop closure, watertightness, half-edge
/// pairing and cycles — and gets the curve-surface check's findings
/// alone.
pub(crate) fn post_check_gates(
    result: &BRep,
    tolerance: &Tolerance,
    carried: &crate::serve_region::CarriedEdges,
    structure_validated: bool,
) -> PostCheckReport {
    if structure_validated {
        let mut report = counted(result);
        push_curve_surface(&mut report, result, carried);
        return report;
    }
    post_check_scoped(result, tolerance, false, false, carried)
}

fn post_check_scoped(
    result: &BRep,
    tolerance: &Tolerance,
    check_self_intersection: bool,
    diagnostics: bool,
    carried: &crate::serve_region::CarriedEdges,
) -> PostCheckReport {
    let mut report = counted(result);

    // 1. Manifold edges
    if let Err(e) = result.validate_manifold_edges() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!("Non-manifold edges: {e}"),
            error: Some(e),
            measured_value: None,
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
            measured_value: Some(report.euler_characteristic as f64),
        });
    }

    // 2b. Loop closure + shell watertightness + curve/surface
    // agreement (PK-342): the coarse-tolerance resonance served
    // bodies whose loops gaped and whose edge curves sat a full cap
    // height off their faces — validate_all catches all three, but
    // the serve gate only consulted the structural set, so the
    // boolean LOGGED "4 issues" and served garbage. These are the
    // geometric closure checks a boolean result must pass to be a
    // boundary representation at all; anything failing them declines
    // as ValidationFailed instead of serving.
    if let Err(e) = result.validate_loop_closure() {
        if popkern_core::trace::trace_enabled("boolean") {
            // Gap-site dump: walk every loop and print curve-end vs
            // next-vertex distances above tolerance scale.
            use popkern_core::geometry::curve::CurveEval as _;
            for (lh, l) in &result.loops {
                let mut cursor = l.half_edge;
                for _ in 0..result.half_edges.len().max(1) {
                    let Some(he) = result.half_edges.get(cursor) else {
                        break;
                    };
                    let Some(edge) = result.edges.get(he.edge) else {
                        break;
                    };
                    if let (Some(c), Some(tw)) = (&edge.curve, result.half_edges.get(he.twin)) {
                        let v0 = result.vertices[he.origin].position;
                        let v1 = result.vertices[tw.origin].position;
                        let p0 = c.geometry.point_at(c.t_start);
                        let p1 = c.geometry.point_at(c.t_end);
                        let d = (p0.distance_to(&v0) + p1.distance_to(&v1))
                            .min(p0.distance_to(&v1) + p1.distance_to(&v0));
                        if d > 1.0e-5 {
                            eprintln!(
                                "  [gap-site] loop {lh:?} he {cursor:?} edge {:?} d={d:.3e} curve ({:.6},{:.6},{:.6})->({:.6},{:.6},{:.6}) v ({:.6},{:.6},{:.6})/({:.6},{:.6},{:.6})",
                                he.edge,
                                p0.x,
                                p0.y,
                                p0.z,
                                p1.x,
                                p1.y,
                                p1.z,
                                v0.x,
                                v0.y,
                                v0.z,
                                v1.x,
                                v1.y,
                                v1.z
                            );
                        }
                    }
                    cursor = he.next;
                    if cursor == l.half_edge {
                        break;
                    }
                }
            }
        }
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!("Loop not geometrically closed: {e}"),
            error: Some(e),
            measured_value: None,
        });
    }
    if let Err(e) = result.validate_watertight() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!("Shell gap: {e}"),
            error: Some(e),
            measured_value: None,
        });
    }
    push_curve_surface(&mut report, result, carried);

    // 3. Orientation consistency
    if let Err(e) = result.validate_half_edge_pairing() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!("Half-edge pairing error: {e}"),
            error: Some(e),
            measured_value: None,
        });
    }

    // 4. Half-edge cycles
    if let Err(e) = result.validate_half_edge_cycles() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!("Half-edge cycle error: {e}"),
            error: Some(e),
            measured_value: None,
        });
    }

    if !diagnostics {
        return report;
    }

    // 5. Vertex consistency
    if let Err(e) = result.validate_vertex_consistency() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Warning,
            message: format!("Vertex consistency: {e}"),
            error: Some(e),
            measured_value: None,
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
            // §13.1: surface the maximum measured gap so callers
            // can compare against their own tolerance budget.
            measured_value: Some(wt.max_gap()),
        });
    }

    // 7. Face orientation
    if let Err(e) = result.validate_face_orientation() {
        report.findings.push(PostCheckFinding {
            severity: Severity::Warning,
            message: format!("Face orientation: {e}"),
            error: Some(e),
            measured_value: None,
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
            measured_value: None,
        });
    }

    report
}

/// An empty report carrying `result`'s entity counts.
fn counted(result: &BRep) -> PostCheckReport {
    PostCheckReport {
        face_count: result.faces.len(),
        edge_count: result.edges.len(),
        vertex_count: result.vertices.len(),
        euler_characteristic: result.vertices.len() as i64 - result.edges.len() as i64
            + result.faces.len() as i64,
        ..Default::default()
    }
}

/// The curve-surface gate (check 9), on every edge but `carried`.
fn push_curve_surface(
    report: &mut PostCheckReport,
    result: &BRep,
    carried: &crate::serve_region::CarriedEdges,
) {
    let curve_surface = if carried.is_empty() {
        result.validate_curve_surface()
    } else {
        let scoped: Vec<_> = result
            .edges
            .keys()
            .filter(|&e| !carried.contains_key(e))
            .collect();
        match result
            .validate_curve_surface_edges(&scoped)
            .into_iter()
            .next()
        {
            Some(e) => Err(e),
            None => Ok(()),
        }
    };
    if let Err(e) = curve_surface {
        report.findings.push(PostCheckFinding {
            severity: Severity::Error,
            message: format!("Edge curve deviates from face surface: {e}"),
            error: Some(e),
            measured_value: None,
        });
    }
}

#[cfg(test)]
#[path = "post_check_tests.rs"]
mod tests;
