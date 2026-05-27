use super::*;

#[ignore]
#[test]
fn direct_apex_cube_mixed_cap_outer_per_edge_breakdown() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let result = build_under_direct_apex(false, &[0.08, 0.10, 0.12]);
    let cap_set = native_cap_faces(&result);
    eprintln!(
        "\n=== direct-apex cube_mixed cap-outer per-edge breakdown: {} cap faces ===",
        cap_set.len()
    );
    for route in &result.repair_report.apex_routes {
        eprintln!(
            "route status={:?} pcurve_pairs={} diagnostics={:?}",
            route.exact_boundary_report.candidate_status,
            route.exact_boundary_pcurve_pair_count,
            route.exact_boundary_report.diagnostics
        );
        for sector in &route.exact_boundary_report.sectors {
            eprintln!(
                "sector={} cap={:?} fillet_edge={:?} face_edge={:?} support=({:?},{:?}) cap_hes={:?}",
                sector.sector_name,
                sector.cap_face,
                sector.fillet_end_edge,
                sector.face_edge,
                sector.fillet_support_face,
                sector.bridge_support_face,
                sector.cap_half_edges
            );
        }
    }
    let kind = |g: &SurfaceGeometry| -> &'static str {
        match g {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cyl",
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Gregory(_) => "Gregory",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
        }
    };
    let pcurve_kind = |p: &popkern_core::geometry::pcurve::Pcurve| -> &'static str {
        match &p.geometry {
            popkern_core::geometry::pcurve::PcurveGeometry::Linear { .. } => "Linear",
            popkern_core::geometry::pcurve::PcurveGeometry::Elliptical { .. } => "Elliptical",
            popkern_core::geometry::pcurve::PcurveGeometry::Sinusoidal { .. } => "Sinusoidal",
            popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(_) => "Nurbs",
            popkern_core::geometry::pcurve::PcurveGeometry::CylinderArcViaRationalQuad {
                ..
            } => "CylinderArcViaRationalQuad",
        }
    };
    for (eh, edge) in result.brep.edges.iter() {
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let he_a = edge.half_edge;
        let Some(ha) = result.brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = result.brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = result.brep.loops[ha.loop_].face;
        let face_b = result.brep.loops[hb.loop_].face;
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if face_a == face_b || a_is_cap == b_is_cap {
            continue;
        }
        let (Some(pa), Some(pb)) = (ha.pcurve.as_ref(), hb.pcurve.as_ref()) else {
            continue;
        };
        let Some(sa) = result
            .brep
            .faces
            .get(face_a)
            .and_then(|f| f.surface.as_ref())
        else {
            continue;
        };
        let Some(sb) = result
            .brep
            .faces
            .get(face_b)
            .and_then(|f| f.surface.as_ref())
        else {
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
        let (pa_forward, pb_forward) = {
            let (t_lo, t_hi) = (coe.t_start, coe.t_end);
            let p3_lo = coe.geometry.point_at(t_lo);
            let p3_hi = coe.geometry.point_at(t_hi);
            (
                p3_lo.distance_to(&pa_at(pa_lo)) + p3_hi.distance_to(&pa_at(pa_hi))
                    <= p3_lo.distance_to(&pa_at(pa_hi)) + p3_hi.distance_to(&pa_at(pa_lo)),
                p3_lo.distance_to(&pb_at(pb_lo)) + p3_hi.distance_to(&pb_at(pb_hi))
                    <= p3_lo.distance_to(&pb_at(pb_hi)) + p3_hi.distance_to(&pb_at(pb_lo)),
            )
        };
        let mut worst = 0.0_f64;
        let mut worst_sample = 0usize;
        for i in 0..100 {
            let s = 0.05 + 0.90 * (i as f64) / 99.0;
            let s_a = if pa_forward { s } else { 1.0 - s };
            let s_b = if pb_forward { s } else { 1.0 - s };
            let t_a = pa_lo + s_a * (pa_hi - pa_lo);
            let t_b = pb_lo + s_b * (pb_hi - pb_lo);
            let gap = pa_at(t_a).distance_to(&pb_at(t_b));
            if gap > worst {
                worst = gap;
                worst_sample = i;
            }
        }
        let t_a = 0.5 * (pa_lo + pa_hi);
        let t_b = 0.5 * (pb_lo + pb_hi);
        let uv_a = pa.point_at(t_a);
        let uv_b = pb.point_at(t_b);
        let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
        let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
        let (t_lo, t_hi) = (coe.t_start, coe.t_end);
        let p_mid = coe.geometry.point_at(0.5 * (t_lo + t_hi));
        eprintln!(
            "edge={eh:?} faces=({face_a:?}:{}, {face_b:?}:{}) cap_is_a={} g0={:.3e} dev_a={:.3e} dev_b={:.3e} curve_mid=({:.3},{:.3},{:.3}) a=({:.3},{:.3},{:.3}) b=({:.3},{:.3},{:.3})",
            kind(&sa.geometry),
            kind(&sb.geometry),
            a_is_cap,
            p_a.distance_to(&p_b),
            p_mid.distance_to(&p_a),
            p_mid.distance_to(&p_b),
            p_mid.x,
            p_mid.y,
            p_mid.z,
            p_a.x,
            p_a.y,
            p_a.z,
            p_b.x,
            p_b.y,
            p_b.z,
        );
        eprintln!(
            "  pcurves=({},{}) ranges=({:.3},{:.3})/({:.3},{:.3}) forward=({},{}) worst100={:.3e}@{}",
            pcurve_kind(pa),
            pcurve_kind(pb),
            pa_lo,
            pa_hi,
            pb_lo,
            pb_hi,
            pa_forward,
            pb_forward,
            worst,
            worst_sample,
        );
    }
}

#[ignore]
#[test]
fn direct_apex_native_cube_mixed_has_cap_to_cylinder_seams() {
    let result = build_under_direct_apex(false, &[0.08, 0.10, 0.12]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "diagnostic expects native direct-apex to be accepted; got fallback(s): {reasons:#?}"
    );

    assert_native_cap_has_cylinder_seams("cube_mixed", &result, 3);
}

#[ignore]
#[test]
fn direct_apex_5_matrix_dump() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: false,
        },
        Case {
            name: "cube_mixed",
            radii: &[0.08, 0.10, 0.12],
            is_pyramid: false,
        },
        Case {
            name: "cube_extreme",
            radii: &[0.02, 0.10, 0.20],
            is_pyramid: false,
        },
        Case {
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
        },
        Case {
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
            is_pyramid: true,
        },
        Case {
            name: "pyramid4_mixed",
            radii: &[0.04, 0.06, 0.05, 0.07],
            is_pyramid: true,
        },
    ];
    eprintln!("\n=== direct-apex validation validate_all matrix natively under direct-apex ===");
    eprintln!(
        "{:<18}  {:>8}  {:>5}  {:>5}  {:>5}  {:>10}",
        "fixture", "fallback", "struct", "pathA", "other", "worst gap"
    );
    eprintln!("{:-<60}", "");
    for case in cases {
        let result = build_under_direct_apex(case.is_pyramid, case.radii);
        let fallbacks = count_direct_apex_fallbacks(&result);
        let errors = result.brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "{:<18}  {:>8}  {:>5}  {:>5}  {:>5}  {:>10}",
            case.name,
            fallbacks,
            m.structural,
            m.path_a_relevant,
            m.other_geometric,
            if m.worst_path_a_gap > 0.0 {
                format!("{:.3e}", m.worst_path_a_gap)
            } else {
                "-".to_string()
            },
        );
        for error in errors {
            eprintln!("  - {error:?}");
        }
    }
}
