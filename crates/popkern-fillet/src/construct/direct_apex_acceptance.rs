//! Native Gregory acceptance checks for DirectApex dispatch.

#[path = "direct_apex_acceptance_helpers.rs"]
mod direct_apex_acceptance_helpers;

use super::*;
pub(super) use direct_apex_acceptance_helpers::native_gregory_has_expected_cylinder_seams;
use direct_apex_acceptance_helpers::{
    measure_exact_boundary_cap_outer_g0, native_gregory_cap_faces,
    trace_loop_for_acceptance_failure, trace_native_gregory_cap_seam_g0,
};

pub(super) fn native_gregory_path_a_error(brep: &BRep) -> Option<String> {
    use popkern_core::error::ValidationError;

    brep.validate_all().into_iter().find_map(|err| match err {
        ValidationError::BrokenCycle { .. } | ValidationError::ShellGap { .. } => {
            Some(format!("{err:?}"))
        }
        ValidationError::OpenLoop { loop_, .. } => {
            if crate::direct_apex_trace::enabled() {
                trace_loop_for_acceptance_failure(brep, loop_);
            }
            Some(format!("{err:?}"))
        }
        _ => None,
    })
}

pub(super) fn native_gregory_cap_outer_g0_error(
    brep: &BRep,
    result: &crate::apex_construction::construct::ApexBuildResult,
    max_gap: f64,
) -> Option<String> {
    use popkern_core::geometry::surface::SurfaceEval;

    let cap_faces = native_gregory_cap_faces(result);
    if cap_faces.is_empty() {
        return Some("no native cap patch faces reported".to_string());
    }

    let cap_face_list = cap_faces.iter().copied().collect::<Vec<_>>();
    if let Some(error) = crate::construct::cap_outer_acceptance::cap_support_acceptance_error(
        brep,
        &cap_face_list,
        max_gap,
        f64::INFINITY,
    ) {
        return Some(error);
    }

    let mut worst_gap = 0.0f64;
    let mut worst_edge = None;
    let mut measured_edges = 0usize;
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
        let a_is_cap = cap_faces.contains(&face_a);
        let b_is_cap = cap_faces.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let (Some(pa), Some(pb)) = (ha.pcurve.as_ref(), hb.pcurve.as_ref()) else {
            continue;
        };
        let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
            continue;
        };

        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let uv_a = pa.point_at(0.5 * (pa_lo + pa_hi));
        let uv_b = pb.point_at(0.5 * (pb_lo + pb_hi));
        let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
        let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
        let gap = p_a.distance_to(&p_b);
        measured_edges += 1;
        if crate::direct_apex_trace::cap_g0_enabled() {
            trace_native_gregory_cap_seam_g0(
                brep,
                edge_handle,
                he_a,
                ha.twin,
                face_a,
                face_b,
                a_is_cap,
                gap,
            );
        }
        if gap > worst_gap {
            worst_gap = gap;
            worst_edge = Some(edge_handle);
        }
        if gap > max_gap {
            if crate::direct_apex_trace::cap_g0_enabled() {
                continue;
            }
            return Some(format!(
                "cap_outer_g0 {gap:.3e} on {edge_handle:?} exceeds {max_gap:.3e}"
            ));
        }
    }

    measure_exact_boundary_cap_outer_g0(
        brep,
        &cap_faces,
        max_gap,
        &mut worst_gap,
        &mut worst_edge,
        &mut measured_edges,
    );

    if measured_edges == 0 {
        return Some("no cap-vs-noncap seams measured".to_string());
    }
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] native Gregory cap_outer_g0 postcondition: worst={worst_gap:.3e} over {measured_edges} edge(s)"
        );
    }
    if worst_gap > max_gap {
        let Some(worst_edge) = worst_edge else {
            return Some(format!(
                "cap_outer_g0 {worst_gap:.3e} exceeds {max_gap:.3e}, but no worst edge was recorded"
            ));
        };
        return Some(format!(
            "cap_outer_g0 {worst_gap:.3e} on {:?} exceeds {max_gap:.3e}",
            worst_edge
        ));
    }
    None
}

pub(super) fn native_gregory_cap_outer_g1_diagnostic(
    brep: &BRep,
    result: &crate::apex_construction::construct::ApexBuildResult,
    max_g0_gap: f64,
    gate: f64,
) -> Option<crate::types::ApexConstructionDiagnostic> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    let cap_faces = native_gregory_cap_faces(result);
    if cap_faces.is_empty() {
        return None;
    }

    let mut worst = 0.0f64;
    let mut worst_g0_gap = 0.0f64;
    let mut worst_edge = None;
    let mut worst_faces = None;
    let mut measured_edges = 0usize;
    for (edge_handle, edge) in brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_faces.contains(&face_a);
        let b_is_cap = cap_faces.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let (Some(pa), Some(pb)) = (ha.pcurve.as_ref(), hb.pcurve.as_ref()) else {
            continue;
        };
        let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
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
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let (t_lo, t_hi) = (coe.t_start, coe.t_end);
        let p3_lo = coe.geometry.point_at(t_lo);
        let p3_hi = coe.geometry.point_at(t_hi);
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
            worst_g0_gap = worst_g0_gap.max(p_a.distance_to(&p_b));
            let n_a = sa
                .geometry
                .partial_u(uv_a.x, uv_a.y)
                .cross(&sa.geometry.partial_v(uv_a.x, uv_a.y));
            let n_b = sb
                .geometry
                .partial_u(uv_b.x, uv_b.y)
                .cross(&sb.geometry.partial_v(uv_b.x, uv_b.y));
            let nn_a = (n_a.x * n_a.x + n_a.y * n_a.y + n_a.z * n_a.z).sqrt();
            let nn_b = (n_b.x * n_b.x + n_b.y * n_b.y + n_b.z * n_b.z).sqrt();
            if nn_a <= 1e-300 || nn_b <= 1e-300 {
                continue;
            }
            let deficit =
                1.0 - ((n_a.x * n_b.x + n_a.y * n_b.y + n_a.z * n_b.z) / (nn_a * nn_b)).abs();
            if deficit > worst {
                worst = deficit;
                worst_edge = Some(edge_handle);
                let (cap_face, support_face) = if a_is_cap {
                    (face_a, face_b)
                } else {
                    (face_b, face_a)
                };
                worst_faces = Some((cap_face, support_face));
            }
        }
    }

    if worst <= gate {
        return None;
    }
    let (cap_face, support_face) = worst_faces?;
    let target = crate::types::CapContinuityTarget::native_cap_g1(max_g0_gap, gate);
    Some(
        crate::types::ApexConstructionDiagnostic::CapSupportContinuity(
            crate::types::CapSupportContinuityDiagnostic {
                continuity: crate::types::CapContinuityReport::declined_g1(
                    target,
                    measured_edges,
                    worst_g0_gap,
                    worst,
                ),
                cap_face,
                support_face,
                edge: worst_edge?,
                worst_g1_coplanarity: worst,
                gate,
                edges_measured: measured_edges,
            },
        ),
    )
}

#[cfg(test)]
mod tests {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::{BRep, EdgeHandle};

    use crate::types::{
        ApexConstructionDiagnostic, CapContinuityReport, CapContinuityTarget, ContinuityLevel,
        FilletOptions, FilletRadius, FilletSpec,
    };

    use super::native_gregory_cap_outer_g1_diagnostic;

    #[test]
    fn native_gregory_g1_decline_populates_cap_continuity_report() {
        let result = mixed_radius_direct_apex_cube();
        let apex_result = crate::apex_construction::construct::ApexBuildResult {
            created_faces: result
                .corner_patches
                .iter()
                .map(|patch| patch.face)
                .collect(),
            apex_vertices: Vec::new(),
            apex_edges: Vec::new(),
            corner_patches: result.corner_patches.clone(),
            exact_boundary_sectors: Vec::new(),
            exact_boundary_prepared_curve_count: 0,
            exact_boundary_inner_spoke_count: 0,
            exact_boundary_ribbon_count: 0,
            exact_boundary_solved_surface_count: 0,
            exact_boundary_pcurve_pair_count: 0,
            exact_boundary_performance: Default::default(),
            exact_boundary_report: result.repair_report.apex_routes[0]
                .exact_boundary_report
                .clone(),
            diagnostics: Vec::new(),
            telemetry: Default::default(),
        };

        let diagnostic =
            native_gregory_cap_outer_g1_diagnostic(&result.brep, &apex_result, 5.0e-3, 0.0)
                .expect("zero G1 gate should force a native cap/support continuity decline");
        let ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) = diagnostic else {
            panic!("expected cap/support continuity diagnostic");
        };

        assert_eq!(
            diagnostic.continuity.target,
            CapContinuityTarget::native_cap_g1(5.0e-3, 0.0)
        );
        assert_eq!(
            diagnostic.continuity,
            CapContinuityReport::declined_g1(
                diagnostic.continuity.target,
                diagnostic.edges_measured,
                diagnostic.continuity.worst_g0_gap.unwrap(),
                diagnostic.worst_g1_coplanarity,
            )
        );
        assert_eq!(diagnostic.continuity.target.requested, ContinuityLevel::G1);
        assert_eq!(
            diagnostic.continuity.target.accepted_floor,
            ContinuityLevel::G1
        );
        assert_eq!(diagnostic.continuity.achieved, Some(ContinuityLevel::G0));
        assert_eq!(
            diagnostic.continuity.decline_code.as_deref(),
            Some("cap_outer_g1_exceeds_gate")
        );
        assert_eq!(
            diagnostic.continuity.measured_edges,
            diagnostic.edges_measured
        );
        assert!(diagnostic.edges_measured > 0);
        assert_eq!(
            diagnostic.continuity.worst_g1_normal_deficit,
            Some(diagnostic.worst_g1_coplanarity)
        );
    }

    fn mixed_radius_direct_apex_cube() -> crate::types::BlendResult {
        let tol = Tolerance::default();
        let cube = BRep::unit_cube(tol).unwrap();
        let apex = cube.vertices.keys().next().unwrap();
        let incident_edges = incident_edges(&cube, apex);
        let specs = incident_edges
            .iter()
            .zip([0.08, 0.10, 0.12])
            .map(|(&edge, radius)| FilletSpec {
                edge,
                radius: FilletRadius::Constant(radius),
            })
            .collect::<Vec<_>>();
        let options = FilletOptions {
            use_direct_apex_construction: true,
            ..FilletOptions::default()
        };
        crate::fillet_edges(&cube, &specs, &options, tol, None).unwrap()
    }

    fn incident_edges(
        brep: &BRep,
        vertex: popkern_core::topology::VertexHandle,
    ) -> Vec<EdgeHandle> {
        brep.half_edges
            .values()
            .filter(|half_edge| half_edge.origin == vertex)
            .map(|half_edge| half_edge.edge)
            .collect()
    }
}
