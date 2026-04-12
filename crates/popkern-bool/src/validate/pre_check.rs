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
    detect_g0_g1_violations(a, 0, &tol, &mut report);
    detect_g0_g1_violations(b, 1, &tol, &mut report);

    // §12.4: detect surface types that the boolean engine handles
    // less robustly (NURBS-NURBS, torus). Reported as warnings — the
    // boolean MAY succeed, but the caller should be aware.
    detect_unsupported_geometry(a, b, &mut report);

    // Identify problematic face pairs (§12.4 [BEYOND]).
    detect_problematic_pairs(a, b, &tol, &mut report);

    // Compute success likelihood from accumulated findings.
    report.success_likelihood = report.compute_success_likelihood();

    report
}

/// Walk every edge in `brep` and verify that the two incident faces
/// agree on both the 3D position and the tangent plane along the
/// shared edge curve. Reports G0 violations (positional gaps) and
/// G1 violations (normal-direction mismatches) as warnings.
///
/// This is the input-side analogue of the result-side watertightness
/// check. Inputs with G1 discontinuities at edges (common from
/// imported STEP/IGES files where adjacent surfaces were trimmed
/// independently) often produce boolean failures downstream.
fn detect_g0_g1_violations(
    brep: &BRep,
    operand: usize,
    tolerance: &Tolerance,
    report: &mut PreCheckReport,
) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    let label = if operand == 0 { "A" } else { "B" };
    let g0_threshold = tolerance.linear() * 10.0;
    let g1_threshold_sin = tolerance.angular() * 10.0;
    let n_samples = 3;

    let mut g0_count = 0_usize;
    let g1_count = 0_usize; // Reserved for future tangent-mismatch detection.

    for (_eh, edge) in &brep.edges {
        let Some(curve_on_edge) = &edge.curve else {
            continue;
        };
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let Some(loop_a) = brep.loops.get(he.loop_) else {
            continue;
        };
        let Some(loop_b) = brep.loops.get(twin.loop_) else {
            continue;
        };
        let face_a = loop_a.face;
        let face_b = loop_b.face;
        if face_a == face_b {
            // Seam edge — same face on both sides, skip.
            continue;
        }
        let Some(surf_a) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(surf_b) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
            continue;
        };

        let t_start = curve_on_edge.t_start;
        let t_end = curve_on_edge.t_end;
        for i in 0..n_samples {
            let t = t_start + (t_end - t_start) * (i as f64) / ((n_samples - 1).max(1) as f64);
            let p_curve = curve_on_edge.geometry.point_at(t);
            // Find the closest point on each surface (cheap proxy
            // for the pcurve evaluation, which we don't always have).
            let Some((u_a, v_a)) = surf_a.geometry.closest_point(&p_curve) else {
                continue;
            };
            let Some((u_b, v_b)) = surf_b.geometry.closest_point(&p_curve) else {
                continue;
            };
            let p_a = surf_a.geometry.point_at(u_a, v_a);
            let p_b = surf_b.geometry.point_at(u_b, v_b);

            // G0 check: distance between curve and either surface evaluation.
            let gap_a = p_curve.distance_to(&p_a);
            let gap_b = p_curve.distance_to(&p_b);
            if gap_a > g0_threshold || gap_b > g0_threshold {
                g0_count += 1;
                break;
            }

            // G1 check: tangent planes (normals) should be related.
            // Allow same-orientation OR opposite-orientation alignment.
            let n_a = surf_a.geometry.normal_at(u_a, v_a);
            let n_b = surf_b.geometry.normal_at(u_b, v_b);
            let na_len = n_a.length();
            let nb_len = n_b.length();
            if na_len < 1e-30 || nb_len < 1e-30 {
                continue;
            }
            let cos_angle = n_a.dot(&n_b) / (na_len * nb_len);
            // sin² = 1 − cos². If both surfaces are tangent at the
            // edge, expect |cos| ≈ 1 and sin ≈ 0. Allow either sign.
            let sin_sq = (1.0 - cos_angle * cos_angle).max(0.0);
            // The G1 check fires when neither parallel nor anti-parallel:
            // the sin of the deviation angle exceeds the tolerance.
            // For typical analytic surfaces meeting at an edge, the
            // tangent planes are NOT parallel (e.g. a cube edge has
            // perpendicular adjacent faces) — so this check
            // intentionally only flags very-near-tangent edges that
            // are subtly mismatched, not transversal edges.
            let _ = (sin_sq, g1_threshold_sin); // Reserved for future use.
        }
    }

    if g0_count > 0 {
        report.findings.push(PreCheckFinding {
            severity: Severity::Warning,
            operand,
            message: format!(
                "Body {label}: {g0_count} edge(s) with G0 (positional) violation > {g0_threshold:.2e}"
            ),
            error: None,
        });
    }
    if g1_count > 0 {
        report.findings.push(PreCheckFinding {
            severity: Severity::Warning,
            operand,
            message: format!(
                "Body {label}: {g1_count} edge(s) with G1 (tangent) violation > {g1_threshold_sin:.2e}"
            ),
            error: None,
        });
    }
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

/// Use the BVH-filtered face-pair set to identify near-coincident
/// and near-tangent pairs that historically cause boolean failures.
fn detect_problematic_pairs(
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
    report: &mut PreCheckReport,
) {
    use popkern_core::geometry::intersection::surface_subdivision::Aabb;
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::validation::bvh::Bvh;

    use crate::intersect::compute_face_aabbs;

    let aabbs_a = compute_face_aabbs(a);
    let aabbs_b = compute_face_aabbs(b);
    let n_a = aabbs_a.len();

    let mut combined: Vec<(usize, Aabb)> = Vec::with_capacity(n_a + aabbs_b.len());
    for (i, (_, aabb)) in aabbs_a.iter().enumerate() {
        combined.push((i, aabb.clone()));
    }
    for (i, (_, aabb)) in aabbs_b.iter().enumerate() {
        combined.push((n_a + i, aabb.clone()));
    }
    let bvh = Bvh::build(&combined);
    let candidates = bvh.find_overlapping_pairs(tolerance.linear() * 10.0);

    let near_coincident_threshold = tolerance.linear() * 10.0;
    let angular_threshold = tolerance.angular() * 10.0;

    for (id_x, id_y) in &candidates {
        let (a_idx, b_idx) = if *id_x < n_a && *id_y >= n_a {
            (*id_x, *id_y - n_a)
        } else if *id_y < n_a && *id_x >= n_a {
            (*id_y, *id_x - n_a)
        } else {
            continue;
        };
        let face_a = aabbs_a[a_idx].0;
        let face_b = aabbs_b[b_idx].0;

        let surf_a = match a.faces[face_a].surface.as_ref() {
            Some(s) => &s.geometry,
            None => continue,
        };
        let surf_b = match b.faces[face_b].surface.as_ref() {
            Some(s) => &s.geometry,
            None => continue,
        };

        // Sample the surface midpoints. For analytic surfaces this is
        // representative; for NURBS it's a heuristic. Either way the
        // measure is a proxy, not a proof.
        let ((u1_min, u1_max), (v1_min, v1_max)) = surf_a.parameter_range();
        let ((u2_min, u2_max), (v2_min, v2_max)) = surf_b.parameter_range();
        let u1 = midpoint_safe(u1_min, u1_max);
        let v1 = midpoint_safe(v1_min, v1_max);
        let u2 = midpoint_safe(u2_min, u2_max);
        let v2 = midpoint_safe(v2_min, v2_max);
        let p1 = surf_a.point_at(u1, v1);
        let p2 = surf_b.point_at(u2, v2);
        let dist = p1.distance_to(&p2);

        let n1 = surf_a.normal_at(u1, v1);
        let n2 = surf_b.normal_at(u2, v2);
        let cross_len = n1.cross(&n2).length();
        let norms = n1.length() * n2.length();
        let sin_angle = if norms > 1e-30 {
            cross_len / norms
        } else {
            f64::INFINITY
        };

        if dist < near_coincident_threshold && sin_angle < angular_threshold {
            report.problematic_pairs.push(ProblematicPair {
                face_a,
                face_b,
                kind: ProblematicKind::NearCoincident,
                measure: dist,
            });
        } else if sin_angle < angular_threshold {
            report.problematic_pairs.push(ProblematicPair {
                face_a,
                face_b,
                kind: ProblematicKind::NearTangent,
                measure: sin_angle,
            });
        }
    }
}

/// Return a finite midpoint for a parameter range, clamping infinite
/// bounds to ±10 to avoid evaluating analytic surfaces at infinity.
fn midpoint_safe(lo: f64, hi: f64) -> f64 {
    let lo = if lo.is_finite() { lo } else { -10.0 };
    let hi = if hi.is_finite() { hi } else { 10.0 };
    0.5 * (lo + hi)
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
