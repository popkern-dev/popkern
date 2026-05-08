
#[test]
#[ignore]
fn dump_per_edge_g1_cube_mixed() {
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
        std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
        std::env::set_var("POPKERN_VAR_RADIUS_REROUTE", "1");
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
            Some(SurfaceGeometry::Sphere { .. }) => Some(fh),
            _ => None,
        })
        .collect();

    eprintln!(
        "BRep: {} faces, {} edges; {} sphere-cap faces",
        brep.faces.len(),
        brep.edges.len(),
        cap_set.len()
    );
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
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
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
        // Measure ALL distinct-face edges, label cap vs non-cap.
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
        let a_label = kind(face_a);
        let b_label = kind(face_b);

        let (t_lo, t_hi) = (coe.t_start, coe.t_end);
        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
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
        let _ = (a_is_cap, b_is_cap);
        entries.push((max_g1, max_g0, eh, face_a, face_b, a_label, b_label));
    }
    // Sort by max_g1 descending.
    entries.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    eprintln!("\n--- per-edge G1/G0 sorted by G1 (top 20) ---");
    for (i, (g1, g0, eh, fa, fb, la, lb)) in entries.iter().take(20).enumerate() {
        let cap_a = cap_set.contains(fa);
        let cap_b = cap_set.contains(fb);
        eprintln!(
            "  [{i:2}] edge={eh:?} g1={g1:.3e} g0={g0:.3e} {fa:?}({la}{}) <-> {fb:?}({lb}{})",
            if cap_a { "*" } else { "" },
            if cap_b { "*" } else { "" }
        );
    }

    // For the top-3 cap-vs-Cyl edges (the cut edges), dump the
    // 3D curve point vs each pcurve-evaluated surface point at
    // a few sample params.  This isolates whether the gap is
    // pcurve-fit error, surface-curve mismatch, or both.
    eprintln!("\n--- top-3 cap-vs-Cyl edges: curve vs pcurve-on-surface ---");
    let mut shown = 0usize;
    for (g1, g0, eh, fa, fb, la, lb) in entries.iter() {
        let cap_a = cap_set.contains(fa);
        let cap_b = cap_set.contains(fb);
        if !(cap_a ^ cap_b) {
            continue;
        }
        if !(la.starts_with("Cyl") || lb.starts_with("Cyl")) {
            continue;
        }
        if shown >= 3 {
            break;
        }
        shown += 1;
        let edge = brep.edges.get(*eh).unwrap();
        let coe = edge.curve.as_ref().unwrap();
        let he_a = edge.half_edge;
        let ha = brep.half_edges.get(he_a).unwrap();
        let hb = brep.half_edges.get(ha.twin).unwrap();
        let pa = ha.pcurve.as_ref().unwrap();
        let pb = hb.pcurve.as_ref().unwrap();
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        let sa = brep.faces[face_a].surface.as_ref().unwrap();
        let sb = brep.faces[face_b].surface.as_ref().unwrap();
        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let (t_lo, t_hi) = (coe.t_start, coe.t_end);

        // Endpoint vertex positions and curve degree.
        let v_origin = brep.vertices[ha.origin].position;
        let v_dest = brep.vertices[hb.origin].position;
        let curve_info = match &coe.geometry {
            popkern_core::geometry::CurveGeometry::Nurbs(nc) => {
                format!(
                    "Nurbs deg={} cps={} weights={:?}",
                    nc.degree,
                    nc.control_points.len(),
                    nc.weights
                )
            }
            popkern_core::geometry::CurveGeometry::Line { .. } => "Line".to_string(),
            popkern_core::geometry::CurveGeometry::Arc { .. } => "Arc".to_string(),
            popkern_core::geometry::CurveGeometry::Circle { .. } => "Circle".to_string(),
            popkern_core::geometry::CurveGeometry::Ellipse { .. } => "Ellipse".to_string(),
            _ => "<other>".to_string(),
        };
        eprintln!("edge={eh:?} ({la} on {fa:?} <-> {lb} on {fb:?}) g1={g1:.3e} g0={g0:.3e}");
        eprintln!(
            "  endpoints: ({:.4},{:.4},{:.4}) → ({:.4},{:.4},{:.4})  curve: {curve_info}",
            v_origin.x, v_origin.y, v_origin.z, v_dest.x, v_dest.y, v_dest.z
        );
        // Surface params for context.
        match &sa.geometry {
            popkern_core::geometry::surface::SurfaceGeometry::Sphere { center, radius } => {
                eprintln!(
                    "  surface_a Sphere centre=({:.4},{:.4},{:.4}) r={radius:.4}",
                    center.x, center.y, center.z
                );
            }
            popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                eprintln!(
                    "  surface_a Cylinder origin=({:.4},{:.4},{:.4}) axis=({:.4},{:.4},{:.4}) r={radius:.4}",
                    origin.x, origin.y, origin.z, axis.x, axis.y, axis.z
                );
            }
            _ => {}
        }
        match &sb.geometry {
            popkern_core::geometry::surface::SurfaceGeometry::Sphere { center, radius } => {
                eprintln!(
                    "  surface_b Sphere centre=({:.4},{:.4},{:.4}) r={radius:.4}",
                    center.x, center.y, center.z
                );
            }
            popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                eprintln!(
                    "  surface_b Cylinder origin=({:.4},{:.4},{:.4}) axis=({:.4},{:.4},{:.4}) r={radius:.4}",
                    origin.x, origin.y, origin.z, axis.x, axis.y, axis.z
                );
            }
            _ => {}
        }
        eprintln!(
            "  curve t∈[{t_lo:.3},{t_hi:.3}] pa_param∈[{pa_lo:.3},{pa_hi:.3}] pb_param∈[{pb_lo:.3},{pb_hi:.3}]"
        );
        // Dump pa/pb pcurve type.
        eprintln!("  ha={he_a:?} hb={:?}", ha.twin);
        match &pa.geometry {
            popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(nc) => {
                eprintln!(
                    "  pa Nurbs deg={} n_cps={} cps={:?}",
                    nc.degree,
                    nc.control_points.len(),
                    nc.control_points
                );
            }
            popkern_core::geometry::pcurve::PcurveGeometry::Linear { start, end } => {
                eprintln!(
                    "  pa Linear ({:.4},{:.4})→({:.4},{:.4}) t∈[{:.4},{:.4}]",
                    start.0, start.1, end.0, end.1, pa.t_start, pa.t_end
                );
            }
            _ => eprintln!("  pa <other pcurve geometry>"),
        }
        if let popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(nc) = &pb.geometry {
            eprintln!(
                "  pb Nurbs deg={} n_cps={} cps={:?}",
                nc.degree,
                nc.control_points.len(),
                nc.control_points
            );
        } else {
            eprintln!("  pb <non-nurbs pcurve geometry>");
        }
        // Detect direction.
        let p3_lo = coe.geometry.point_at(t_lo);
        let p3_hi = coe.geometry.point_at(t_hi);
        let pa_at_param = |t: f64| {
            let uv = pa.point_at(t);
            sa.geometry.point_at(uv.x, uv.y)
        };
        let pb_at_param = |t: f64| {
            let uv = pb.point_at(t);
            sb.geometry.point_at(uv.x, uv.y)
        };
        let pa_forward = p3_lo.distance_to(&pa_at_param(pa_lo))
            + p3_hi.distance_to(&pa_at_param(pa_hi))
            <= p3_lo.distance_to(&pa_at_param(pa_hi)) + p3_hi.distance_to(&pa_at_param(pa_lo));
        let pb_forward = p3_lo.distance_to(&pb_at_param(pb_lo))
            + p3_hi.distance_to(&pb_at_param(pb_hi))
            <= p3_lo.distance_to(&pb_at_param(pb_hi)) + p3_hi.distance_to(&pb_at_param(pb_lo));

        for s in [0.0_f64, 0.25, 0.5, 0.75, 1.0] {
            let t = t_lo + s * (t_hi - t_lo);
            let p3 = coe.geometry.point_at(t);
            let s_a = if pa_forward { s } else { 1.0 - s };
            let s_b = if pb_forward { s } else { 1.0 - s };
            let t_a = pa_lo + s_a * (pa_hi - pa_lo);
            let t_b = pb_lo + s_b * (pb_hi - pb_lo);
            let uv_a = pa.point_at(t_a);
            let uv_b = pb.point_at(t_b);
            let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
            let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
            let d_curve_a = p3.distance_to(&p_a);
            let d_curve_b = p3.distance_to(&p_b);
            let d_a_b = p_a.distance_to(&p_b);
            eprintln!(
                "    s={s:.2}  curve=({:.4},{:.4},{:.4}) |curve→a|={d_curve_a:.3e} |curve→b|={d_curve_b:.3e} |a→b|={d_a_b:.3e}",
                p3.x, p3.y, p3.z
            );
        }
    }
}
