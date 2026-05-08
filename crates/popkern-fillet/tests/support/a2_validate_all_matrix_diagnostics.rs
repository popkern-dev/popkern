use super::*;

#[ignore]
#[test]
fn a2_cube_mixed_cap_outer_per_edge_breakdown() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let result = build_under_a2(false, &[0.08, 0.10, 0.12]);
    let cap_set = result
        .corner_patches
        .iter()
        .map(|p| p.face)
        .collect::<std::collections::BTreeSet<_>>();
    eprintln!(
        "\n=== A2 cube_mixed cap-outer per-edge breakdown: {} cap faces ===",
        cap_set.len()
    );
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
    }
}

#[ignore]
#[test]
fn a2_native_cube_mixed_has_cap_to_cylinder_seams() {
    let result = build_under_a2(false, &[0.08, 0.10, 0.12]);
    let reasons = a2_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "diagnostic expects native A2 to be accepted; got fallback(s): {reasons:#?}"
    );

    assert_native_cap_has_cylinder_seams("cube_mixed", &result, 3);
}

#[ignore]
#[test]
fn a2_5_matrix_dump() {
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
    eprintln!("\n=== A2.5 validate_all matrix natively under A2 ===");
    eprintln!(
        "{:<18}  {:>8}  {:>5}  {:>5}  {:>5}  {:>10}",
        "fixture", "fallback", "struct", "pathA", "other", "worst gap"
    );
    eprintln!("{:-<60}", "");
    for case in cases {
        let result = build_under_a2(case.is_pyramid, case.radii);
        let fallbacks = count_a2_fallbacks(&result);
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
