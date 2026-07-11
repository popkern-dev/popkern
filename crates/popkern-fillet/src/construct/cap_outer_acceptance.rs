//! Cap/support postcondition checks shared by apex cap candidates.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

#[path = "cap_outer_acceptance_measure.rs"]
mod measure;
use measure::{CapSupportSeamEvidence, measure_cap_support_continuity};

/// Validate that a cap candidate preserves the cap/support boundary
/// postconditions expected by commercial-grade apex construction.
pub(crate) fn cap_support_acceptance_error(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    max_g0_gap: f64,
    max_g1_coplanarity: f64,
) -> Option<String> {
    let cap_set = cap_faces
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let mut worst_g0 = 0.0f64;
    let mut worst_g1 = 0.0f64;
    let mut worst_g0_edge = None;
    let mut worst_g1_edge = None;
    let mut worst_g0_seam = None;
    let mut worst_g1_seam = None;
    let mut expected_edges = 0usize;
    let mut measured_edges = 0usize;
    let mut missing_pcurve_edges = 0usize;
    let mut missing_cap_pcurve_edges = 0usize;
    let mut missing_support_pcurve_edges = 0usize;
    let mut missing_surface_edges = 0usize;
    let mut missing_curve_edges = 0usize;

    for (edge_handle, edge) in brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let Some(face_a) = brep.loops.get(ha.loop_).map(|loop_data| loop_data.face) else {
            continue;
        };
        let Some(face_b) = brep.loops.get(hb.loop_).map(|loop_data| loop_data.face) else {
            continue;
        };
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        expected_edges += 1;
        // Per-edge forensics for the aggregate counts in the decline
        // message: which seam edge, which side lacks the pcurve, and
        // what surfaces flank it.
        let trace_missing = |missing_side: &str, side_is_cap: bool| {
            if !crate::direct_apex_trace::enabled() {
                return;
            }
            let kind = |fh: FaceHandle| {
                brep.faces
                    .get(fh)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| match &s.geometry {
                        popkern_core::geometry::SurfaceGeometry::Plane { .. } => "Plane",
                        popkern_core::geometry::SurfaceGeometry::Cylinder { .. } => "Cylinder",
                        popkern_core::geometry::SurfaceGeometry::Sphere { .. } => "Sphere",
                        popkern_core::geometry::SurfaceGeometry::Cone { .. } => "Cone",
                        popkern_core::geometry::SurfaceGeometry::Torus { .. } => "Torus",
                        popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "Nurbs",
                        popkern_core::geometry::SurfaceGeometry::Gregory(_) => "Gregory",
                    })
                    .unwrap_or("None")
            };
            let p = brep
                .half_edges
                .get(he_a)
                .and_then(|h| brep.vertices.get(h.origin))
                .map(|v| {
                    format!(
                        "({:.4},{:.4},{:.4})",
                        v.position.x, v.position.y, v.position.z
                    )
                })
                .unwrap_or_default();
            eprintln!(
                "[direct-apex-trace] seam {edge_handle:?} missing {missing_side} pcurve \
                 (missing side is_cap={side_is_cap}) faces {face_a:?}:{} / {face_b:?}:{} at {p}",
                kind(face_a),
                kind(face_b),
            );
        };
        let Some(pa) = ha.pcurve.as_ref() else {
            missing_pcurve_edges += 1;
            if a_is_cap {
                missing_cap_pcurve_edges += 1;
            } else {
                missing_support_pcurve_edges += 1;
            }
            trace_missing("half-edge-a", a_is_cap);
            continue;
        };
        let Some(pb) = hb.pcurve.as_ref() else {
            missing_pcurve_edges += 1;
            if b_is_cap {
                missing_cap_pcurve_edges += 1;
            } else {
                missing_support_pcurve_edges += 1;
            }
            trace_missing("half-edge-b", b_is_cap);
            continue;
        };
        let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            missing_surface_edges += 1;
            continue;
        };
        let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
            missing_surface_edges += 1;
            continue;
        };
        let Some(coe) = edge.curve.as_ref() else {
            missing_curve_edges += 1;
            continue;
        };

        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let pa_at = |t: f64| {
            let uv = pa.point_at(t);
            sa.geometry.point_at(uv.x, uv.y)
        };
        let pb_at = |t: f64| {
            let uv = pb.point_at(t);
            sb.geometry.point_at(uv.x, uv.y)
        };
        let p3_lo = coe.geometry.point_at(coe.t_start);
        let p3_hi = coe.geometry.point_at(coe.t_end);
        let pa_forward = p3_lo.distance_to(&pa_at(pa_lo)) + p3_hi.distance_to(&pa_at(pa_hi))
            <= p3_lo.distance_to(&pa_at(pa_hi)) + p3_hi.distance_to(&pa_at(pa_lo));
        let pb_forward = p3_lo.distance_to(&pb_at(pb_lo)) + p3_hi.distance_to(&pb_at(pb_hi))
            <= p3_lo.distance_to(&pb_at(pb_hi)) + p3_hi.distance_to(&pb_at(pb_lo));

        measured_edges += 1;
        for i in 0..100 {
            let s = 0.05 + 0.90 * (i as f64) / 99.0;
            let s_a = if pa_forward { s } else { 1.0 - s };
            let s_b = if pb_forward { s } else { 1.0 - s };
            let uv_a = pa.point_at(pa_lo + s_a * (pa_hi - pa_lo));
            let uv_b = pb.point_at(pb_lo + s_b * (pb_hi - pb_lo));
            let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
            let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
            let curve_point = coe
                .geometry
                .point_at(coe.t_start + s * (coe.t_end - coe.t_start));
            let (cap_curve_gap, support_curve_gap) = if a_is_cap {
                (p_a.distance_to(&curve_point), p_b.distance_to(&curve_point))
            } else {
                (p_b.distance_to(&curve_point), p_a.distance_to(&curve_point))
            };
            let g0_gap = p_a.distance_to(&p_b);
            if g0_gap > worst_g0 {
                worst_g0 = g0_gap;
                worst_g0_edge = Some(edge_handle);
                let (cap_face, support_face) = if a_is_cap {
                    (face_a, face_b)
                } else {
                    (face_b, face_a)
                };
                let (cap_uv, support_uv) = if a_is_cap {
                    ((uv_a.x, uv_a.y), (uv_b.x, uv_b.y))
                } else {
                    ((uv_b.x, uv_b.y), (uv_a.x, uv_a.y))
                };
                worst_g0_seam = Some(CapSupportSeamEvidence {
                    cap_face,
                    support_face,
                    edge: edge_handle,
                    sample_parameter: s,
                    cap_uv,
                    support_uv,
                    cap_curve_gap,
                    support_curve_gap,
                });
            }

            let (normal_u_a, normal_v_a) =
                interior_normal_sample_parameters(&sa.geometry, uv_a.x, uv_a.y);
            let (normal_u_b, normal_v_b) =
                interior_normal_sample_parameters(&sb.geometry, uv_b.x, uv_b.y);
            let n_a = sa
                .geometry
                .partial_u(normal_u_a, normal_v_a)
                .cross(&sa.geometry.partial_v(normal_u_a, normal_v_a));
            let n_b = sb
                .geometry
                .partial_u(normal_u_b, normal_v_b)
                .cross(&sb.geometry.partial_v(normal_u_b, normal_v_b));
            let len_a = n_a.length();
            let len_b = n_b.length();
            if len_a > 1e-300 && len_b > 1e-300 {
                let dot = n_a.x * n_b.x + n_a.y * n_b.y + n_a.z * n_b.z;
                let g1_gap = 1.0 - (dot / (len_a * len_b)).abs();
                if g1_gap > worst_g1 {
                    worst_g1 = g1_gap;
                    worst_g1_edge = Some(edge_handle);
                    let (cap_face, support_face) = if a_is_cap {
                        (face_a, face_b)
                    } else {
                        (face_b, face_a)
                    };
                    let (cap_uv, support_uv) = if a_is_cap {
                        ((uv_a.x, uv_a.y), (uv_b.x, uv_b.y))
                    } else {
                        ((uv_b.x, uv_b.y), (uv_a.x, uv_a.y))
                    };
                    worst_g1_seam = Some(CapSupportSeamEvidence {
                        cap_face,
                        support_face,
                        edge: edge_handle,
                        sample_parameter: s,
                        cap_uv,
                        support_uv,
                        cap_curve_gap,
                        support_curve_gap,
                    });
                }
            }
        }
    }

    if expected_edges == 0 {
        return Some("measured 0 cap/support seams, expected at least 1".to_string());
    }
    if measured_edges < expected_edges {
        return Some(format!(
            "measured {measured_edges} cap/support seams, expected {expected_edges} topological cap/support seams; missing_pcurve_edges={missing_pcurve_edges}; missing_cap_pcurve_edges={missing_cap_pcurve_edges}; missing_support_pcurve_edges={missing_support_pcurve_edges}; missing_surface_edges={missing_surface_edges}; missing_curve_edges={missing_curve_edges}",
        ));
    }
    if worst_g0 > max_g0_gap {
        let Some(seam) = worst_g0_seam else {
            return Some(format!(
                "cap_outer_g0 {worst_g0:.3e} on {:?} exceeds {max_g0_gap:.3e}",
                worst_g0_edge.unwrap_or_default()
            ));
        };
        return Some(format!(
            "cap_outer_g0 {worst_g0:.3e} on {:?} cap_face={:?} support_face={:?} sample={:.3} cap_uv=({:.3},{:.3}) support_uv=({:.3},{:.3}) cap_curve_gap={:.3e} support_curve_gap={:.3e} exceeds {max_g0_gap:.3e}",
            seam.edge,
            seam.cap_face,
            seam.support_face,
            seam.sample_parameter,
            seam.cap_uv.0,
            seam.cap_uv.1,
            seam.support_uv.0,
            seam.support_uv.1,
            seam.cap_curve_gap,
            seam.support_curve_gap
        ));
    }
    if worst_g1 > max_g1_coplanarity {
        let Some(seam) = worst_g1_seam else {
            return Some(format!(
                "cap_outer_g1 {worst_g1:.3e} on {:?} exceeds {max_g1_coplanarity:.3e}",
                worst_g1_edge.unwrap_or_default()
            ));
        };
        return Some(format!(
            "cap_outer_g1 {worst_g1:.3e} on {:?} cap_face={:?} support_face={:?} sample={:.3} cap_uv=({:.3},{:.3}) support_uv=({:.3},{:.3}) cap_curve_gap={:.3e} support_curve_gap={:.3e} exceeds {max_g1_coplanarity:.3e}",
            seam.edge,
            seam.cap_face,
            seam.support_face,
            seam.sample_parameter,
            seam.cap_uv.0,
            seam.cap_uv.1,
            seam.support_uv.0,
            seam.support_uv.1,
            seam.cap_curve_gap,
            seam.support_curve_gap
        ));
    }
    None
}

/// Promote cap/support seam acceptance failure into the shared diagnostics
/// surface used by direct-apex and exact-boundary routes.
pub(crate) fn cap_support_continuity_diagnostic(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    max_g0_gap: f64,
    max_g1_coplanarity: f64,
) -> Option<crate::types::ApexConstructionDiagnostic> {
    let evidence = measure_cap_support_continuity(brep, cap_faces)?;
    if evidence.worst_g0_gap <= max_g0_gap && evidence.worst_g1_coplanarity <= max_g1_coplanarity {
        return None;
    }
    let target = crate::types::CapContinuityTarget::native_cap_g1(max_g0_gap, max_g1_coplanarity);
    let continuity = if evidence.worst_g0_gap > max_g0_gap {
        crate::types::CapContinuityReport::declined_g0(
            target,
            evidence.measured_edges,
            evidence.worst_g0_gap,
            evidence.worst_g1_coplanarity,
        )
    } else {
        crate::types::CapContinuityReport::declined_g1(
            target,
            evidence.measured_edges,
            evidence.worst_g0_gap,
            evidence.worst_g1_coplanarity,
        )
    };
    let failed_seam = if evidence.worst_g0_gap > max_g0_gap {
        evidence.worst_g0_seam
    } else {
        evidence.worst_g1_seam
    };
    Some(
        crate::types::ApexConstructionDiagnostic::CapSupportContinuity(
            crate::types::CapSupportContinuityDiagnostic {
                continuity,
                cap_face: failed_seam.cap_face,
                support_face: failed_seam.support_face,
                edge: failed_seam.edge,
                measurement_scope: "live_committed_topology".to_string(),
                sample_parameter: Some(failed_seam.sample_parameter),
                cap_uv: Some(failed_seam.cap_uv),
                support_uv: Some(failed_seam.support_uv),
                cap_curve_gap: Some(failed_seam.cap_curve_gap),
                support_curve_gap: Some(failed_seam.support_curve_gap),
                worst_g1_coplanarity: evidence.worst_g1_coplanarity,
                gate: max_g1_coplanarity,
                edges_measured: evidence.measured_edges,
            },
        ),
    )
}

fn interior_normal_sample_parameters(surface: &impl SurfaceEval, u: f64, v: f64) -> (f64, f64) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    (
        interior_parameter(u, u_min, u_max),
        interior_parameter(v, v_min, v_max),
    )
}

fn interior_parameter(value: f64, min: f64, max: f64) -> f64 {
    if !value.is_finite() || !min.is_finite() || !max.is_finite() || max <= min {
        return value;
    }
    let margin = ((max - min).abs() * 1.0e-5).max(1.0e-5);
    value.clamp(min + margin, max - margin)
}
