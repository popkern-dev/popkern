use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

pub(super) struct CapSupportContinuityEvidence {
    pub(super) worst_g0_gap: f64,
    pub(super) worst_g1_coplanarity: f64,
    pub(super) measured_edges: usize,
    pub(super) worst_g0_seam: CapSupportSeamEvidence,
    pub(super) worst_g1_seam: CapSupportSeamEvidence,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct CapSupportSeamEvidence {
    pub(super) cap_face: FaceHandle,
    pub(super) support_face: FaceHandle,
    pub(super) edge: EdgeHandle,
    pub(super) sample_parameter: f64,
    pub(super) cap_uv: (f64, f64),
    pub(super) support_uv: (f64, f64),
    pub(super) cap_curve_gap: f64,
    pub(super) support_curve_gap: f64,
}

pub(super) fn measure_cap_support_continuity(
    brep: &BRep,
    cap_faces: &[FaceHandle],
) -> Option<CapSupportContinuityEvidence> {
    let cap_set = cap_faces
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let mut evidence: Option<CapSupportContinuityEvidence> = None;
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
        let Some(sample) = sample_cap_support_edge(brep, edge_handle, face_a, face_b, a_is_cap)
        else {
            continue;
        };
        match &mut evidence {
            Some(existing) => {
                existing.measured_edges += 1;
                if sample.worst_g0_gap > existing.worst_g0_gap {
                    existing.worst_g0_gap = sample.worst_g0_gap;
                    existing.worst_g0_seam = sample.worst_g0_seam;
                }
                if sample.worst_g1_coplanarity > existing.worst_g1_coplanarity {
                    existing.worst_g1_coplanarity = sample.worst_g1_coplanarity;
                    existing.worst_g1_seam = sample.worst_g1_seam;
                }
            }
            None => evidence = Some(sample),
        }
    }
    evidence
}

fn sample_cap_support_edge(
    brep: &BRep,
    edge_handle: EdgeHandle,
    face_a: FaceHandle,
    face_b: FaceHandle,
    a_is_cap: bool,
) -> Option<CapSupportContinuityEvidence> {
    let edge = brep.edges.get(edge_handle)?;
    let ha = brep.half_edges.get(edge.half_edge)?;
    let hb = brep.half_edges.get(ha.twin)?;
    let (Some(pa), Some(pb)) = (ha.pcurve.as_ref(), hb.pcurve.as_ref()) else {
        return None;
    };
    let sa = brep.faces.get(face_a).and_then(|f| f.surface.as_ref())?;
    let sb = brep.faces.get(face_b).and_then(|f| f.surface.as_ref())?;
    let coe = edge.curve.as_ref()?;
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
    let mut worst_g0_gap = 0.0f64;
    let mut worst_g1_coplanarity = 0.0f64;
    let mut worst_g0_parameter = 0.0f64;
    let mut worst_g1_parameter = 0.0f64;
    let mut worst_g0_cap_uv = (0.0f64, 0.0f64);
    let mut worst_g0_support_uv = (0.0f64, 0.0f64);
    let mut worst_g1_cap_uv = (0.0f64, 0.0f64);
    let mut worst_g1_support_uv = (0.0f64, 0.0f64);
    let mut worst_g0_cap_curve_gap = 0.0f64;
    let mut worst_g0_support_curve_gap = 0.0f64;
    let mut worst_g1_cap_curve_gap = 0.0f64;
    let mut worst_g1_support_curve_gap = 0.0f64;
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
        if g0_gap > worst_g0_gap {
            worst_g0_gap = g0_gap;
            worst_g0_parameter = s;
            let (cap_uv, support_uv) = cap_and_support_uv(a_is_cap, uv_a, uv_b);
            worst_g0_cap_uv = cap_uv;
            worst_g0_support_uv = support_uv;
            worst_g0_cap_curve_gap = cap_curve_gap;
            worst_g0_support_curve_gap = support_curve_gap;
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
            let g1_coplanarity = 1.0 - (dot / (len_a * len_b)).abs();
            if g1_coplanarity > worst_g1_coplanarity {
                worst_g1_coplanarity = g1_coplanarity;
                worst_g1_parameter = s;
                let (cap_uv, support_uv) = cap_and_support_uv(a_is_cap, uv_a, uv_b);
                worst_g1_cap_uv = cap_uv;
                worst_g1_support_uv = support_uv;
                worst_g1_cap_curve_gap = cap_curve_gap;
                worst_g1_support_curve_gap = support_curve_gap;
            }
        }
    }
    let (cap_face, support_face) = if a_is_cap {
        (face_a, face_b)
    } else {
        (face_b, face_a)
    };
    Some(CapSupportContinuityEvidence {
        worst_g0_gap,
        worst_g1_coplanarity,
        measured_edges: 1,
        worst_g0_seam: CapSupportSeamEvidence {
            cap_face,
            support_face,
            edge: edge_handle,
            sample_parameter: worst_g0_parameter,
            cap_uv: worst_g0_cap_uv,
            support_uv: worst_g0_support_uv,
            cap_curve_gap: worst_g0_cap_curve_gap,
            support_curve_gap: worst_g0_support_curve_gap,
        },
        worst_g1_seam: CapSupportSeamEvidence {
            cap_face,
            support_face,
            edge: edge_handle,
            sample_parameter: worst_g1_parameter,
            cap_uv: worst_g1_cap_uv,
            support_uv: worst_g1_support_uv,
            cap_curve_gap: worst_g1_cap_curve_gap,
            support_curve_gap: worst_g1_support_curve_gap,
        },
    })
}

fn cap_and_support_uv(
    a_is_cap: bool,
    uv_a: popkern_core::geometry::Point3,
    uv_b: popkern_core::geometry::Point3,
) -> ((f64, f64), (f64, f64)) {
    if a_is_cap {
        ((uv_a.x, uv_a.y), (uv_b.x, uv_b.y))
    } else {
        ((uv_b.x, uv_b.y), (uv_a.x, uv_a.y))
    }
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
