//! Cap-outer boundary metric helpers for K-P commercial gates.
//! Included only in test mode; .unwrap() is allowed in tests.
// #[cfg(test)] — marker for the no_unwrap linter.

use popkern_core::topology::{BRep, FaceHandle};

use super::{dot, vec_norm};

/// Worst-case continuity at the *outer* boundary of the corner-cap
/// ring, where each cap sub-face meets the cylinder fillet roll
/// or other non-cap neighbour it sits on top of.
pub struct CapOuterMetrics {
    /// Max 3-D gap `||cyl(cyl_pc(t)) - greg(greg_pc(t))||`.
    pub worst_g0_gap: f64,
    /// Max unit-normal parallelism deficit at edge samples.
    pub worst_g1_coplanarity: f64,
    /// Number of cap-vs-non-cap edges measured.
    pub edges_measured: usize,
}

/// Measure G0/G1 at every edge between a Gregory cap face and a
/// non-cap neighbour.
///
/// Edges without a 3D curve are skipped: the curve-less
/// cap/support seam class this metric once had a permissive mode
/// for (Stage 4.7e's `include_pcurve_only_edges`, retired with
/// PK-181) is now a hard failure in
/// `direct_apex_native_cap_support_seam_ownership_inventory` and
/// L1 of `validate_coherence`.
pub fn measure_cap_outer_boundary(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    n_samples: usize,
) -> CapOuterMetrics {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let mut worst_g0 = 0.0_f64;
    let mut worst_g1 = 0.0_f64;
    let mut edges_measured = 0usize;
    let mut worst_g0_trace = String::new();
    let mut worst_g1_trace = String::new();
    let n_samples = n_samples.max(2);

    for (_eh, edge) in brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let Some(loop_a) = brep.loops.get(ha.loop_) else {
            continue;
        };
        let Some(loop_b) = brep.loops.get(hb.loop_) else {
            continue;
        };
        let face_a = loop_a.face;
        let face_b = loop_b.face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let cap_half_edge = if a_is_cap { he_a } else { ha.twin };
        let cap_face = if a_is_cap { face_a } else { face_b };
        if !half_edge_is_on_face_outer_loop(brep, cap_face, cap_half_edge) {
            continue;
        }
        let Some(pa) = ha.pcurve.as_ref() else {
            continue;
        };
        let Some(pb) = hb.pcurve.as_ref() else {
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
        let (pa_forward, pb_forward) = if let Some(coe) = edge.curve.as_ref() {
            let (t_lo, t_hi) = (coe.t_start, coe.t_end);
            let p3_lo = coe.geometry.point_at(t_lo);
            let p3_hi = coe.geometry.point_at(t_hi);
            (
                p3_lo.distance_to(&pa_at(pa_lo)) + p3_hi.distance_to(&pa_at(pa_hi))
                    <= p3_lo.distance_to(&pa_at(pa_hi)) + p3_hi.distance_to(&pa_at(pa_lo)),
                p3_lo.distance_to(&pb_at(pb_lo)) + p3_hi.distance_to(&pb_at(pb_hi))
                    <= p3_lo.distance_to(&pb_at(pb_hi)) + p3_hi.distance_to(&pb_at(pb_lo)),
            )
        } else {
            continue;
        };

        edges_measured += 1;
        for i in 0..n_samples {
            let s = 0.05 + 0.90 * (i as f64) / ((n_samples - 1) as f64);
            let s_a = if pa_forward { s } else { 1.0 - s };
            let s_b = if pb_forward { s } else { 1.0 - s };
            let t_a = pa_lo + s_a * (pa_hi - pa_lo);
            let t_b = pb_lo + s_b * (pb_hi - pb_lo);
            let uv_a = pa.point_at(t_a);
            let uv_b = pb.point_at(t_b);
            let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
            let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
            let g0 = p_a.distance_to(&p_b);
            if g0 > worst_g0 {
                worst_g0 = g0;
                if popkern_core::trace::trace_enabled("cap-outer-measure") {
                    let edge_trace = edge
                        .curve
                        .as_ref()
                        .map(|coe| {
                            let p0 = coe.geometry.point_at(coe.t_start);
                            let p1 = coe.geometry.point_at(coe.t_end);
                            format!(
                                " curve=({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                                p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                            )
                        })
                        .unwrap_or_default();
                    worst_g0_trace = format!(
                        "edge={_eh:?} faces=({face_a:?}:{},{face_b:?}:{}) cap={} sample={i}{edge_trace} uv_a=({:.3},{:.3}) uv_b=({:.3},{:.3}) p_a=({:.4},{:.4},{:.4}) p_b=({:.4},{:.4},{:.4})",
                        surface_kind_label(&sa.geometry),
                        surface_kind_label(&sb.geometry),
                        if a_is_cap { "a" } else { "b" },
                        uv_a.x,
                        uv_a.y,
                        uv_b.x,
                        uv_b.y,
                        p_a.x,
                        p_a.y,
                        p_a.z,
                        p_b.x,
                        p_b.y,
                        p_b.z
                    );
                }
            }

            let (normal_u_a, normal_v_a) =
                interior_normal_sample_parameters(&sa.geometry, uv_a.x, uv_a.y);
            let (normal_u_b, normal_v_b) =
                interior_normal_sample_parameters(&sb.geometry, uv_b.x, uv_b.y);
            let su_a = sa.geometry.partial_u(normal_u_a, normal_v_a);
            let sv_a = sa.geometry.partial_v(normal_u_a, normal_v_a);
            let su_b = sb.geometry.partial_u(normal_u_b, normal_v_b);
            let sv_b = sb.geometry.partial_v(normal_u_b, normal_v_b);
            let n_a = su_a.cross(&sv_a);
            let n_b = su_b.cross(&sv_b);
            let nn_a = vec_norm(n_a);
            let nn_b = vec_norm(n_b);
            if nn_a > 1e-300 && nn_b > 1e-300 {
                let deficit = 1.0 - (dot(n_a, n_b) / (nn_a * nn_b)).abs();
                if deficit > worst_g1 {
                    worst_g1 = deficit;
                    if popkern_core::trace::trace_enabled("cap-outer-measure") {
                        worst_g1_trace = format!(
                            "edge={_eh:?} faces=({face_a:?},{face_b:?}) sample={i} uv_a=({:.3},{:.3}) uv_b=({:.3},{:.3}) n_a=({:.3},{:.3},{:.3}) n_b=({:.3},{:.3},{:.3}) su_a=({:.3},{:.3},{:.3}) sv_a=({:.3},{:.3},{:.3}) su_b=({:.3},{:.3},{:.3}) sv_b=({:.3},{:.3},{:.3})",
                            uv_a.x,
                            uv_a.y,
                            uv_b.x,
                            uv_b.y,
                            n_a.x / nn_a,
                            n_a.y / nn_a,
                            n_a.z / nn_a,
                            n_b.x / nn_b,
                            n_b.y / nn_b,
                            n_b.z / nn_b,
                            su_a.x,
                            su_a.y,
                            su_a.z,
                            sv_a.x,
                            sv_a.y,
                            sv_a.z,
                            su_b.x,
                            su_b.y,
                            su_b.z,
                            sv_b.x,
                            sv_b.y,
                            sv_b.z
                        );
                    }
                }
            }
        }
    }
    if !worst_g0_trace.is_empty() {
        eprintln!("[cap_outer measure worst_g0] {worst_g0:.6e} {worst_g0_trace}");
    }
    if !worst_g1_trace.is_empty() {
        eprintln!("[cap_outer measure worst_g1] {worst_g1:.6e} {worst_g1_trace}");
    }

    CapOuterMetrics {
        worst_g0_gap: worst_g0,
        worst_g1_coplanarity: worst_g1,
        edges_measured,
    }
}

fn surface_kind_label(surface: &popkern_core::geometry::surface::SurfaceGeometry) -> &'static str {
    use popkern_core::geometry::surface::SurfaceGeometry;
    match surface {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
        SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
    }
}

fn interior_normal_sample_parameters(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    u: f64,
    v: f64,
) -> (f64, f64) {
    use popkern_core::geometry::surface::SurfaceEval;

    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    // Periodic directions WRAP instead of clamping: a loop legally
    // chained onto an adjacent period branch (u - 2π on a cylinder
    // chart) evaluates exactly, but clamping such a parameter to the
    // domain edge samples the normal up to a full period away
    // (PK-125: pyramid4_mixed read a phantom cap_outer_g1 of 0.99
    // through a clamped support sample).
    let (period_u, period_v) = surface.periodicity();
    (
        interior_parameter(wrap_periodic(u, u_min, u_max, period_u), u_min, u_max),
        interior_parameter(wrap_periodic(v, v_min, v_max, period_v), v_min, v_max),
    )
}

fn wrap_periodic(value: f64, min: f64, max: f64, period: Option<f64>) -> f64 {
    let Some(period) = period else {
        return value;
    };
    if !(value.is_finite() && min.is_finite() && max.is_finite()) || period <= 0.0 {
        return value;
    }
    let mut wrapped = value;
    while wrapped < min {
        wrapped += period;
    }
    while wrapped > max {
        wrapped -= period;
    }
    wrapped
}

fn interior_parameter(value: f64, min: f64, max: f64) -> f64 {
    if !value.is_finite() || !min.is_finite() || !max.is_finite() || max <= min {
        return value;
    }
    let margin = ((max - min).abs() * 1.0e-5).max(1.0e-5);
    value.clamp(min + margin, max - margin)
}

fn half_edge_is_on_face_outer_loop(
    brep: &BRep,
    face: FaceHandle,
    half_edge: popkern_core::topology::HalfEdgeHandle,
) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return false;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        if cursor == half_edge {
            return true;
        }
        let Some(half_edge_data) = brep.half_edges.get(cursor) else {
            return false;
        };
        cursor = half_edge_data.next;
        if cursor == start {
            return false;
        }
    }
    false
}
