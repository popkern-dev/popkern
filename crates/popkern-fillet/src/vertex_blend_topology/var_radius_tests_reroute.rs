
#[test]
#[ignore]
fn reroute_full_cube_mixed() {
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
        std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
        std::env::set_var("POPKERN_VAR_RADIUS_REROUTE", "1");
        std::env::set_var("POPKERN_VAR_RADIUS_TRACE", "1");
    }
    let cube =
        popkern_core::topology::BRep::unit_cube(popkern_core::tolerance::Tolerance::default())
            .unwrap();
    let apex = cube
        .vertices
        .keys()
        .find(|&v| cube.half_edges.values().filter(|he| he.origin == v).count() == 3)
        .unwrap();
    let edges: Vec<_> = {
        use std::collections::BTreeSet;
        let s: BTreeSet<_> = cube
            .half_edges
            .values()
            .filter(|he| he.origin == apex)
            .map(|he| he.edge)
            .collect();
        s.into_iter().collect()
    };
    let radii = [0.08_f64, 0.10, 0.12];
    let specs: Vec<crate::types::FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| crate::types::FilletSpec {
            edge: e,
            radius: crate::types::FilletRadius::Constant(r),
        })
        .collect();
    let result = crate::fillet_edges(
        &cube,
        &specs,
        &crate::FilletOptions::default(),
        popkern_core::tolerance::Tolerance::default(),
        None,
    );
    match result {
        Ok(filleted) => {
            eprintln!(
                "OK — final BRep: {} faces, {} edges, {} vertices",
                filleted.brep.faces.len(),
                filleted.brep.edges.len(),
                filleted.brep.vertices.len()
            );
            if let Err(e) = filleted.brep.validate() {
                eprintln!("post-pipeline validate FAILED: {e:?}");
            } else {
                eprintln!("post-pipeline validate OK");
            }
        }
        Err(e) => {
            eprintln!("fillet_edges failed: {e:?}");
        }
    }
}

/// Phase 6c-D Step 2 probe — measure post-reroute G1/G0 for
/// Gregory caps (rung-OFF mode), to size the gap that
/// re-optimisation needs to close.  Mirror of
/// `dump_per_edge_g1_cube_mixed` but with `POPKERN_SPHERE_CAP_RUNG`
/// unset (so cap surfaces stay Gregory after K-P binonic ladder
/// fits them pre-reroute) and `POPKERN_OUTER_FILLET_G1` on
/// (so the optimiser includes outer-fillet G1 residuals).
///
/// Reroute is gated by env: when `POPKERN_VAR_RADIUS_REROUTE=1`
/// (default in this test) we measure post-reroute state.
/// Set `POPKERN_PHASE_6CD_PROBE_NO_REROUTE=1` to instead measure
/// pre-reroute baseline (separates optimiser quality from
/// reroute damage).  `#[ignore]`d because it sets process-level
/// env vars.
#[test]
#[ignore]
fn dump_per_edge_g1_cube_mixed_gregory_rung_off() {
    let no_reroute = std::env::var("POPKERN_PHASE_6CD_PROBE_NO_REROUTE")
        .map(|v| v == "1")
        .unwrap_or(false);
    let outer_g1 = std::env::var("POPKERN_PHASE_6CD_PROBE_OUTER_G1")
        .map(|v| v == "1")
        .unwrap_or(true);
    let sphere_fit = std::env::var("POPKERN_PHASE_6CD_PROBE_SPHERE_FIT")
        .map(|v| v == "1")
        .unwrap_or(true);
    unsafe {
        if sphere_fit {
            std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
        } else {
            std::env::remove_var("POPKERN_SPHERE_FIT_BOUNDARY");
        }
        std::env::remove_var("POPKERN_SPHERE_CAP_RUNG");
        if outer_g1 {
            std::env::set_var("POPKERN_OUTER_FILLET_G1", "1");
        } else {
            std::env::remove_var("POPKERN_OUTER_FILLET_G1");
        }
        if no_reroute {
            std::env::remove_var("POPKERN_VAR_RADIUS_REROUTE");
        } else {
            std::env::set_var("POPKERN_VAR_RADIUS_REROUTE", "1");
        }
    }
    eprintln!(
        "Phase 6c-D probe mode: {} OUTER_FILLET_G1={} SPHERE_FIT={}",
        if no_reroute {
            "PRE-REROUTE"
        } else {
            "POST-REROUTE"
        },
        if outer_g1 { "ON" } else { "OFF" },
        if sphere_fit { "ON" } else { "OFF" }
    );
    let cube =
        popkern_core::topology::BRep::unit_cube(popkern_core::tolerance::Tolerance::default())
            .unwrap();
    let apex = cube
        .vertices
        .keys()
        .find(|&v| cube.half_edges.values().filter(|he| he.origin == v).count() == 3)
        .unwrap();
    let edges: Vec<_> = {
        use std::collections::BTreeSet;
        let s: BTreeSet<_> = cube
            .half_edges
            .values()
            .filter(|he| he.origin == apex)
            .map(|he| he.edge)
            .collect();
        s.into_iter().collect()
    };
    let radii = [0.08_f64, 0.10, 0.12];
    let specs: Vec<crate::types::FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| crate::types::FilletSpec {
            edge: e,
            radius: crate::types::FilletRadius::Constant(r),
        })
        .collect();
    let result = crate::fillet_edges(
        &cube,
        &specs,
        &crate::FilletOptions::default(),
        popkern_core::tolerance::Tolerance::default(),
        None,
    );
    let Ok(filleted) = result else {
        eprintln!("fillet_edges failed");
        return;
    };
    let brep = filleted.brep;

    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
    let cap_set: std::collections::BTreeSet<_> = brep
        .faces
        .iter()
        .filter_map(|(fh, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Gregory(_)) => Some(fh),
            _ => None,
        })
        .collect();

    eprintln!(
        "BRep (rung-OFF Gregory): {} faces, {} edges; {} Gregory-cap faces",
        brep.faces.len(),
        brep.edges.len(),
        cap_set.len()
    );

    let kind = |face: popkern_core::topology::FaceHandle| -> &'static str {
        match brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        {
            Some(SurfaceGeometry::Plane { .. }) => "Plane",
            Some(SurfaceGeometry::Cylinder { .. }) => "Cyl",
            Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
            Some(SurfaceGeometry::Torus { .. }) => "Torus",
            Some(SurfaceGeometry::Cone { .. }) => "Cone",
            Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
            Some(SurfaceGeometry::Gregory(_)) => "Greg",
            None => "?",
        }
    };

    let mut entries: Vec<(
        f64,
        f64,
        popkern_core::topology::EdgeHandle,
        popkern_core::topology::FaceHandle,
        popkern_core::topology::FaceHandle,
        &str,
        &str,
    )> = Vec::new();
    for (eh, edge) in brep.edges.iter() {
        if edge.curve.as_ref().is_none() {
            continue;
        }
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
        let a_label = kind(face_a);
        let b_label = kind(face_b);
        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let coe = edge.curve.as_ref().unwrap();
        let (t_lo, t_hi) = (coe.t_start, coe.t_end);
        let p3_lo = coe.geometry.point_at(t_lo);
        let p3_hi = coe.geometry.point_at(t_hi);
        let pa_at = |t: f64| {
            let uv = pa.point_at(t);
            sa.geometry.point_at(uv.x, uv.y)
        };
        let pb_at = |t: f64| {
            let uv = pb.point_at(t);
            sb.geometry.point_at(uv.x, uv.y)
        };
        let pa_forward = p3_lo.distance_to(&pa_at(pa_lo)) + p3_hi.distance_to(&pa_at(pa_hi))
            <= p3_lo.distance_to(&pa_at(pa_hi)) + p3_hi.distance_to(&pa_at(pa_lo));
        let pb_forward = p3_lo.distance_to(&pb_at(pb_lo)) + p3_hi.distance_to(&pb_at(pb_hi))
            <= p3_lo.distance_to(&pb_at(pb_hi)) + p3_hi.distance_to(&pb_at(pb_lo));

        let n_samples = 7usize;
        let mut max_g0 = 0.0_f64;
        let mut max_g1 = 0.0_f64;
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
            if g0 > max_g0 {
                max_g0 = g0;
            }
            let su_a = sa.geometry.partial_u(uv_a.x, uv_a.y);
            let sv_a = sa.geometry.partial_v(uv_a.x, uv_a.y);
            let su_b = sb.geometry.partial_u(uv_b.x, uv_b.y);
            let sv_b = sb.geometry.partial_v(uv_b.x, uv_b.y);
            let n_a = su_a.cross(&sv_a);
            let n_b = su_b.cross(&sv_b);
            let nn_a = (n_a.dot(&n_a)).sqrt();
            let nn_b = (n_b.dot(&n_b)).sqrt();
            if nn_a > 1e-300 && nn_b > 1e-300 {
                let dot_unit = (n_a.dot(&n_b) / (nn_a * nn_b)).abs();
                let deficit = 1.0 - dot_unit;
                if deficit > max_g1 {
                    max_g1 = deficit;
                }
            }
        }
        entries.push((max_g1, max_g0, eh, face_a, face_b, a_label, b_label));
    }
    entries.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    eprintln!("\n--- per-edge G1/G0 sorted by G1 (top 30) ---");
    for (i, (g1, g0, eh, fa, fb, la, lb)) in entries.iter().take(30).enumerate() {
        let cap_a = cap_set.contains(fa);
        let cap_b = cap_set.contains(fb);
        eprintln!(
            "  [{i:2}] edge={eh:?} g1={g1:.3e} g0={g0:.3e} {fa:?}({la}{}) <-> {fb:?}({lb}{})",
            if cap_a { "*" } else { "" },
            if cap_b { "*" } else { "" }
        );
    }

    // Aggregate stats: cap-vs-non-cap (worst-case post-reroute
    // gap) vs other (sanity baseline).
    let (mut max_cap_g1, mut max_cap_g0) = (0.0_f64, 0.0_f64);
    let (mut max_other_g1, mut max_other_g0) = (0.0_f64, 0.0_f64);
    for (g1, g0, _, fa, fb, _, _) in entries.iter() {
        let cap_a = cap_set.contains(fa);
        let cap_b = cap_set.contains(fb);
        if cap_a ^ cap_b {
            if *g1 > max_cap_g1 {
                max_cap_g1 = *g1;
            }
            if *g0 > max_cap_g0 {
                max_cap_g0 = *g0;
            }
        } else {
            if *g1 > max_other_g1 {
                max_other_g1 = *g1;
            }
            if *g0 > max_other_g0 {
                max_other_g0 = *g0;
            }
        }
    }
    eprintln!(
        "\nPhase 6c-D Step 2 sizing — Gregory rung-OFF, post-reroute:\n\
             cap↔non-cap: max g1={max_cap_g1:.3e} max g0={max_cap_g0:.3e}\n\
             other:       max g1={max_other_g1:.3e} max g0={max_other_g0:.3e}"
    );
}
