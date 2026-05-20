//! Cap/support postcondition checks shared by apex cap candidates.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

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
    let mut measured_edges = 0usize;

    for (_edge_handle, edge) in brep.edges.iter() {
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
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
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
        let Some(coe) = edge.curve.as_ref() else {
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
            worst_g0 = worst_g0.max(p_a.distance_to(&p_b));

            let n_a = sa
                .geometry
                .partial_u(uv_a.x, uv_a.y)
                .cross(&sa.geometry.partial_v(uv_a.x, uv_a.y));
            let n_b = sb
                .geometry
                .partial_u(uv_b.x, uv_b.y)
                .cross(&sb.geometry.partial_v(uv_b.x, uv_b.y));
            let len_a = n_a.length();
            let len_b = n_b.length();
            if len_a > 1e-300 && len_b > 1e-300 {
                let dot = n_a.x * n_b.x + n_a.y * n_b.y + n_a.z * n_b.z;
                worst_g1 = worst_g1.max(1.0 - (dot / (len_a * len_b)).abs());
            }
        }
    }

    if measured_edges < cap_faces.len() {
        return Some(format!(
            "measured {measured_edges} cap/support seams, expected at least {}",
            cap_faces.len()
        ));
    }
    if worst_g0 > max_g0_gap {
        return Some(format!(
            "cap_outer_g0 {worst_g0:.3e} exceeds {max_g0_gap:.3e}"
        ));
    }
    if worst_g1 > max_g1_coplanarity {
        return Some(format!(
            "cap_outer_g1 {worst_g1:.3e} exceeds {max_g1_coplanarity:.3e}"
        ));
    }
    None
}
