//! Geometric pre-check passes (G0/G1 continuity, problematic-pair
//! identification). Split out of `pre_check.rs` to keep the parent
//! under the 500-line structural cap.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::Severity;

use super::pre_check::{PreCheckFinding, PreCheckReport, ProblematicKind, ProblematicPair};

/// Walk every edge in `brep` and verify that the two incident faces
/// agree on both the 3D position and the tangent plane along the
/// shared edge curve. Reports G0 violations (positional gaps) and
/// G1 violations (normal-direction mismatches) as warnings.
///
/// This is the input-side analogue of the result-side watertightness
/// check. Inputs with G1 discontinuities at edges (common from
/// imported STEP/IGES files where adjacent surfaces were trimmed
/// independently) often produce boolean failures downstream.
pub(super) fn detect_g0_g1_violations(
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
    let mut g1_count = 0_usize;

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

            // G1 check: tangent-plane mismatch on near-tangent edges.
            //
            // Two surfaces meeting at an edge are G1-continuous when
            // their tangent planes (normals) align. For a SHARP edge
            // (e.g. cube corner), the surfaces meet at a large angle
            // and discontinuity is expected — no violation.
            //
            // We only flag the SUBTLE case: surfaces that are
            // nearly-parallel (designed to be tangent) but slightly
            // off. The detection: |cos(angle)| > 0.99 (so the edge
            // SHOULD be smooth) AND sin(angle) > tolerance.
            let n_a = surf_a.geometry.normal_at(u_a, v_a);
            let n_b = surf_b.geometry.normal_at(u_b, v_b);
            let na_len = n_a.length();
            let nb_len = n_b.length();
            if na_len < 1e-30 || nb_len < 1e-30 {
                continue;
            }
            let cos_angle = n_a.dot(&n_b) / (na_len * nb_len);
            let abs_cos = cos_angle.abs();
            // Only check edges that are designed to be smooth (within
            // ~8° of tangent: |cos| > 0.99). Sharp transitions
            // (cube edges, etc.) are excluded by this gate.
            if abs_cos > 0.99 {
                let sin_sq = (1.0 - cos_angle * cos_angle).max(0.0);
                if sin_sq > g1_threshold_sin * g1_threshold_sin {
                    g1_count += 1;
                    break;
                }
            }
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

/// Use the BVH-filtered face-pair set to identify near-coincident
/// and near-tangent pairs that historically cause boolean failures.
pub(super) fn detect_problematic_pairs(
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
