//! PK-095: whole-object tessellation of `fillet_edges` output — the
//! per-face diagnostic probe and (eventually) the closure gate.

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn pentagon_blend() -> BRep {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| {
            a.1.position
                .z
                .partial_cmp(&b.1.position.z)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (origin == apex || twin_origin == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.16),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let result = fillet_edges(&pyramid, &specs, &options, tol, None).expect("pipeline");
    result.brep.validate().expect("valid blend");
    result.brep
}

/// Count boundary edges of a face mesh after welding (1e-6).
fn face_mesh_boundary(mesh: &crate::FaceMesh) -> (usize, usize) {
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut remap = Vec::new();
    for v in &mesh.vertices {
        let k = (
            (v.position.x * 1.0e6).round() as i64,
            (v.position.y * 1.0e6).round() as i64,
            (v.position.z * 1.0e6).round() as i64,
        );
        let next = key.len();
        remap.push(*key.entry(k).or_insert(next));
    }
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    let mut degen = 0;
    for t in &mesh.triangles {
        let w = [remap[t[0]], remap[t[1]], remap[t[2]]];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            degen += 1;
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    // A correct face mesh is an open disk: its boundary edges should
    // form exactly one loop tracing the trim boundary; count edges
    // with odd counts as a proxy for internal tearing.
    let torn = counts.values().filter(|&&c| c > 2).count();
    (degen, torn)
}

#[test]
#[ignore = "PK-095 diagnostic: per-face tessellation report for the pentagon blend"]
fn per_face_tessellation_report() {
    let brep = pentagon_blend();
    for (fh, face) in brep.faces.iter() {
        let kind = match face.surface.as_ref().map(|s| &s.geometry) {
            Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
            Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
            Some(popkern_core::geometry::SurfaceGeometry::Gregory(_)) => "Gregory",
            Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "Cylinder",
            _ => "Other",
        };
        let loops = crate::trim::build_trim_loops(
            &brep,
            fh,
            &face.surface.as_ref().unwrap().geometry,
            crate::TessellationOptions::default().chord_height,
        );
        let loop_info: Vec<String> = loops
            .iter()
            .map(|l| {
                format!(
                    "pts={} from_pcurves={} outer={}",
                    l.points.len(),
                    l.from_pcurves,
                    l.is_outer
                )
            })
            .collect();
        match crate::tessellate::tessellate_face(&brep, fh, &crate::TessellationOptions::default())
        {
            Ok(mesh) => {
                let (degen, torn) = face_mesh_boundary(&mesh);
                eprintln!(
                    "face {fh:?} {kind}: tris={} degen={degen} torn_edges={torn} loops=[{}]",
                    mesh.triangles.len(),
                    loop_info.join("; ")
                );
            }
            Err(e) => eprintln!("face {fh:?} {kind}: TESSELLATION FAILED: {e}"),
        }
    }
}

#[test]
#[ignore = "PK-095 diagnostic: half-edge pcurve detail for one shredded plane face"]
fn lateral_face_pcurve_detail() {
    let brep = pentagon_blend();
    // Pick the worst plane face: fewest loop points relative to edges.
    for (fh, face) in brep.faces.iter() {
        if !matches!(
            face.surface.as_ref().map(|s| &s.geometry),
            Some(popkern_core::geometry::SurfaceGeometry::Plane { .. })
        ) {
            continue;
        }
        let loop_ = face.outer_loop;
        let start = brep.loops[loop_].half_edge;
        let mut he = start;
        let mut n_edges = 0;
        loop {
            n_edges += 1;
            he = brep.half_edges[he].next;
            if he == start {
                break;
            }
        }
        eprintln!("plane face {fh:?}: {n_edges} loop half-edges");
        let mut he = start;
        loop {
            let data = &brep.half_edges[he];
            let a = brep.vertices[data.origin].position;
            let b = brep.vertices[brep.half_edges[data.next].origin].position;
            let pdesc = match &data.pcurve {
                Some(pc) => match &pc.geometry {
                    popkern_core::geometry::pcurve::PcurveGeometry::Linear { .. } => {
                        "Linear".to_string()
                    }
                    popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(n) => {
                        format!("Nurbs(deg={},cps={})", n.degree, n.control_points.len())
                    }
                    _ => "OtherPc".to_string(),
                },
                None => "NONE".to_string(),
            };
            let mut samples: Vec<(f64, f64)> = Vec::new();
            if let Some(pc) = &data.pcurve {
                crate::trim::adaptive_sample_pcurve(pc, 0.001, &mut samples);
            }
            let surface = &face.surface.as_ref().unwrap().geometry;
            let true_a = popkern_core::geometry::uv_projection::project_point_to_uv(surface, &a);
            let uv_err = match (samples.first(), true_a) {
                (Some(&(su, sv)), Some((ta, tb))) => {
                    (((su - ta).powi(2) + (sv - tb).powi(2)) as f64).sqrt()
                }
                _ => f64::NAN,
            };
            eprintln!(
                "  he {he:?}: len3d={:.4} pcurve={pdesc} sampled_pts={} pcurve_uv={:?} true_uv={:?} uv_err={uv_err:.3e}",
                a.distance_to(&b),
                samples.len(),
                samples.first(),
                true_a,
            );
            he = data.next;
            if he == start {
                break;
            }
        }
    }
}

#[test]
#[ignore = "PK-095 diagnostic: trim loop winding per face"]
fn trim_loop_winding_report() {
    let brep = pentagon_blend();
    for (fh, face) in brep.faces.iter() {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        let loops = crate::trim::build_trim_loops(
            &brep,
            fh,
            &surface.geometry,
            crate::TessellationOptions::default().chord_height,
        );
        let kind = match &surface.geometry {
            popkern_core::geometry::SurfaceGeometry::Plane { .. } => "Plane",
            popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "Nurbs",
            popkern_core::geometry::SurfaceGeometry::Gregory(_) => "Gregory",
            _ => "Other",
        };
        for l in &loops {
            let mut area2 = 0.0;
            for i in 0..l.points.len() {
                let (x1, y1) = l.points[i];
                let (x2, y2) = l.points[(i + 1) % l.points.len()];
                area2 += x1 * y2 - x2 * y1;
            }
            eprintln!(
                "face {fh:?} {kind} orient={} loop outer={} signed_area={:.4}",
                surface.orientation,
                l.is_outer,
                0.5 * area2
            );
        }
    }
}

#[test]
#[ignore = "PK-095 diagnostic: classification detail for one lateral plane face"]
fn lateral_face_classification_detail() {
    let brep = pentagon_blend();
    for (fh, face) in brep.faces.iter() {
        if !matches!(
            face.surface.as_ref().map(|s| &s.geometry),
            Some(popkern_core::geometry::SurfaceGeometry::Plane { .. })
        ) {
            continue;
        }
        // Lateral faces have exactly 3 loop half-edges.
        let start = brep.loops[face.outer_loop].half_edge;
        let mut n = 0;
        let mut he = start;
        loop {
            n += 1;
            he = brep.half_edges[he].next;
            if he == start {
                break;
            }
        }
        if n != 3 {
            continue;
        }
        let surface = &face.surface.as_ref().unwrap().geometry;
        let loops = crate::trim::build_trim_loops(
            &brep,
            fh,
            surface,
            crate::TessellationOptions::default().chord_height,
        );
        eprintln!("face {fh:?}: trim loop points:");
        for l in &loops {
            for pt in &l.points {
                eprintln!("    ({:.4}, {:.4})", pt.0, pt.1);
            }
        }
        // Reproduce the grid + classification tallies.
        let mesh =
            crate::tessellate::tessellate_face(&brep, fh, &crate::TessellationOptions::default())
                .expect("face mesh");
        // UV extents of the resulting mesh vertices:
        let us: Vec<f64> = mesh.vertices.iter().map(|v| v.uv.0).collect();
        let vs: Vec<f64> = mesh.vertices.iter().map(|v| v.uv.1).collect();
        let minmax = |xs: &[f64]| {
            (
                xs.iter().cloned().fold(f64::INFINITY, f64::min),
                xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            )
        };
        eprintln!(
            "  mesh: {} verts, {} tris, u {:?} v {:?}",
            mesh.vertices.len(),
            mesh.triangles.len(),
            minmax(&us),
            minmax(&vs)
        );
        // Area coverage: sum of output triangle UV areas vs the trim
        // polygon's area. A complete clip covers ~100%.
        let mut tri_area = 0.0;
        for t in &mesh.triangles {
            let a = mesh.vertices[t[0]].uv;
            let b = mesh.vertices[t[1]].uv;
            let c = mesh.vertices[t[2]].uv;
            tri_area += 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
        }
        let mut poly_area = 0.0;
        if let Some(l) = loops.first() {
            for i in 0..l.points.len() {
                let (x1, y1) = l.points[i];
                let (x2, y2) = l.points[(i + 1) % l.points.len()];
                poly_area += 0.5 * (x1 * y2 - x2 * y1);
            }
        }
        eprintln!(
            "  coverage: tris {:.4} / polygon {:.4} = {:.1}%",
            tri_area,
            poly_area.abs(),
            100.0 * tri_area / poly_area.abs().max(1e-12)
        );
        break;
    }
}

#[test]
#[ignore = "PK-095 diagnostic: strip trim polygon extents"]
fn strip_polygon_extents() {
    let brep = pentagon_blend();
    for (fh, face) in brep.faces.iter() {
        if !matches!(
            face.surface.as_ref().map(|s| &s.geometry),
            Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_))
        ) {
            continue;
        }
        let surface = &face.surface.as_ref().unwrap().geometry;
        let loops = crate::trim::build_trim_loops(
            &brep,
            fh,
            surface,
            crate::TessellationOptions::default().chord_height,
        );
        let Some(l) = loops.first() else { continue };
        let umin = l.points.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
        let umax = l
            .points
            .iter()
            .map(|p| p.0)
            .fold(f64::NEG_INFINITY, f64::max);
        let vmin = l.points.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let vmax = l
            .points
            .iter()
            .map(|p| p.1)
            .fold(f64::NEG_INFINITY, f64::max);
        // 3D positions of the polygon's v extremes vs the surface's v=0/1 edges.
        use popkern_core::geometry::surface::SurfaceEval;
        let mid_u = 0.5 * (umin + umax);
        let poly_at_vmin = surface.point_at(mid_u, vmin);
        let true_v0 = surface.point_at(mid_u, 0.0);
        eprintln!(
            "strip {fh:?}: u [{umin:.3},{umax:.3}] v [{vmin:.3},{vmax:.3}] | poly-vmin vs v=0 dist {:.4}",
            poly_at_vmin.distance_to(&true_v0)
        );
    }
}

#[test]
#[ignore = "PK-095 diagnostic: stray face anatomy"]
fn stray_face_anatomy() {
    let brep = pentagon_blend();
    eprintln!("total faces: {}", brep.faces.len());
    for (fh, face) in brep.faces.iter() {
        let start = brep.loops[face.outer_loop].half_edge;
        let mut n = 0;
        let mut zmin = f64::INFINITY;
        let mut zmax = f64::NEG_INFINITY;
        let mut he = start;
        loop {
            n += 1;
            let z = brep.vertices[brep.half_edges[he].origin].position.z;
            zmin = zmin.min(z);
            zmax = zmax.max(z);
            he = brep.half_edges[he].next;
            if he == start {
                break;
            }
        }
        let kind = match face.surface.as_ref().map(|s| &s.geometry) {
            Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
            Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
            Some(popkern_core::geometry::SurfaceGeometry::Gregory(_)) => "Gregory",
            _ => "Other",
        };
        if kind == "Nurbs" {
            eprintln!("face {fh:?} {kind}: {n} edges, z [{zmin:.3},{zmax:.3}]");
        }
    }
}

#[test]
#[ignore = "PK-095 diagnostic: which side is off the true trim"]
fn trim_side_attribution() {
    use popkern_core::geometry::surface::SurfaceEval;
    let brep = pentagon_blend();
    // Find one fillet strip (Nurbs, 10 edges) and the base-adjacent
    // lateral planes; compare each mesh's boundary to the strip's
    // v=0 and v=1 iso curves (the contact lines).
    for (fh, face) in brep.faces.iter() {
        let Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) =
            face.surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        let surface = &face.surface.as_ref().unwrap().geometry;
        let mesh =
            crate::tessellate::tessellate_face(&brep, fh, &crate::TessellationOptions::default())
                .expect("mesh");
        // Contact lines sampled.
        let mut contact: Vec<popkern_core::geometry::Point3> = Vec::new();
        for i in 0..=64 {
            let u = i as f64 / 64.0;
            contact.push(surface.point_at(u, 0.0));
            contact.push(surface.point_at(u, 1.0));
        }
        // Strip mesh boundary points (edge-count method).
        let mut ec = std::collections::HashMap::new();
        for t in &mesh.triangles {
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                *ec.entry((a.min(b), a.max(b))).or_insert(0usize) += 1;
            }
        }
        let mut bpts = std::collections::HashSet::new();
        for (&(a, b), &c) in &ec {
            if c == 1 {
                bpts.insert(a);
                bpts.insert(b);
            }
        }
        let mut worst = 0.0_f64;
        for &i in &bpts {
            let p = mesh.vertices[i].position;
            if p.z < 0.05 || p.z > 0.95 {
                continue;
            }
            let d = contact
                .iter()
                .map(|q| q.distance_to(&p))
                .fold(f64::INFINITY, f64::min);
            worst = worst.max(d);
        }
        eprintln!("strip {fh:?}: worst boundary-to-contact-line dist = {worst:.5}");
        let vmin = mesh
            .vertices
            .iter()
            .map(|v| v.uv.1)
            .fold(f64::INFINITY, f64::min);
        let vmax = mesh
            .vertices
            .iter()
            .map(|v| v.uv.1)
            .fold(f64::NEG_INFINITY, f64::max);
        eprintln!("  mesh uv v-extent: [{vmin:.4},{vmax:.4}]");
        // All count-1 edges with interior endpoints.
        for (&(a, b), &c) in &ec {
            if c != 1 {
                continue;
            }
            let ua = mesh.vertices[a].uv;
            let ub = mesh.vertices[b].uv;
            let interior = |uv: (f64, f64)| uv.1 > 0.02 && uv.1 < 0.98 && uv.0 > 0.17;
            if interior(ua) || interior(ub) {
                eprintln!(
                    "  open edge ({:.3},{:.3})-({:.3},{:.3})",
                    ua.0, ua.1, ub.0, ub.1
                );
            }
        }
        // One worst offender in detail.
        for &i in &bpts {
            let p = mesh.vertices[i].position;
            if p.z < 0.05 || p.z > 0.95 {
                continue;
            }
            let d = contact
                .iter()
                .map(|q| q.distance_to(&p))
                .fold(f64::INFINITY, f64::min);
            if (d - worst).abs() < 1e-9 {
                let uv = mesh.vertices[i].uv;
                let reeval = surface.point_at(uv.0, uv.1);
                eprintln!(
                    "  worst pt: uv=({:.4},{:.4}) pos=({:.4},{:.4},{:.4}) reeval-dist={:.2e}",
                    uv.0,
                    uv.1,
                    p.x,
                    p.y,
                    p.z,
                    reeval.distance_to(&p)
                );
            }
        }
        break;
    }
}

#[test]
#[ignore = "PK-095 diagnostic: strip surface iso vs recorded trim curve"]
fn strip_iso_vs_trim_curve() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let brep = pentagon_blend();
    for (fh, face) in brep.faces.iter() {
        let Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) =
            face.surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        let surface = &face.surface.as_ref().unwrap().geometry;
        // For every loop edge: worst distance from the curve polyline
        // to the surface's nearest of the four iso boundaries.
        let start = brep.loops[face.outer_loop].half_edge;
        let mut he = start;
        loop {
            let data = &brep.half_edges[he];
            let eh = data.edge;
            he = data.next;
            if let Some(curve) = brep.edges[eh].curve.as_ref() {
                let ceval: &dyn CurveEval = &curve.geometry;
                let mut worst = 0.0_f64;
                for i in 0..=16 {
                    let t = curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 16.0;
                    let p = ceval.point_at(t);
                    // nearest of dense samples on the 4 iso borders
                    let mut best = f64::INFINITY;
                    for j in 0..=64 {
                        let w = j as f64 / 64.0;
                        for q in [
                            surface.point_at(w, 0.0),
                            surface.point_at(w, 1.0),
                            surface.point_at(0.0, w),
                            surface.point_at(1.0, w),
                        ] {
                            best = best.min(q.distance_to(&p));
                        }
                    }
                    worst = worst.max(best);
                }
                eprintln!("strip {fh:?} edge {eh:?}: curve-to-iso-border worst {worst:.5}");
            }
            if he == start {
                break;
            }
        }
        break;
    }
}

#[test]
#[ignore = "PK-095 diagnostic: point inversion on the sweep"]
fn sweep_point_inversion_works() {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let brep = pentagon_blend();
    for (fh, face) in brep.faces.iter() {
        let Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) =
            face.surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        let surface = &face.surface.as_ref().unwrap().geometry;
        let p = surface.point_at(0.3, 0.5);
        match project_point_to_uv(surface, &p) {
            Some((u, v)) => {
                let q = surface.point_at(u, v);
                eprintln!(
                    "strip {fh:?}: invert(0.3,0.5) -> ({u:.4},{v:.4}) roundtrip {:.2e}",
                    q.distance_to(&p)
                );
            }
            None => eprintln!("strip {fh:?}: inversion FAILED"),
        }
        // And a trim curve point.
        let start = brep.loops[face.outer_loop].half_edge;
        let eh = brep.half_edges[start].edge;
        if let Some(c) = brep.edges[eh].curve.as_ref() {
            use popkern_core::geometry::curve::CurveEval;
            let ce: &dyn CurveEval = &c.geometry;
            let m = ce.point_at(0.5 * (c.t_start + c.t_end));
            match project_point_to_uv(surface, &m) {
                Some((u, v)) => {
                    let q = surface.point_at(u, v);
                    eprintln!(
                        "  trim mid -> ({u:.4},{v:.4}) surf-dist {:.2e}",
                        q.distance_to(&m)
                    );
                }
                None => eprintln!("  trim mid inversion FAILED"),
            }
        }
        break;
    }
}

#[test]
#[ignore = "PK-095 diagnostic: final strip pcurve agreement"]
fn final_strip_pcurve_agreement() {
    use popkern_core::geometry::surface::SurfaceEval;
    let brep = pentagon_blend();
    for (fh, face) in brep.faces.iter() {
        let Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) =
            face.surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        let surface = &face.surface.as_ref().unwrap().geometry;
        let start = brep.loops[face.outer_loop].half_edge;
        let mut he = start;
        loop {
            let data = &brep.half_edges[he];
            let eh = data.edge;
            let pc_desc = match &data.pcurve {
                Some(pc) => {
                    let kind = match &pc.geometry {
                        popkern_core::geometry::pcurve::PcurveGeometry::Linear { .. } => {
                            "Linear".to_string()
                        }
                        popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(n) => {
                            format!("Nurbs(cps={})", n.control_points.len())
                        }
                        _ => "Other".to_string(),
                    };
                    let mut worst = 0.0_f64;
                    if let Some(c) = brep.edges[eh].curve.as_ref() {
                        use popkern_core::geometry::curve::CurveEval;
                        let ce: &dyn CurveEval = &c.geometry;
                        for f in [0.25_f64, 0.5, 0.75] {
                            let tp = pc.t_start + (pc.t_end - pc.t_start) * f;
                            let uv = pc.point_at(tp);
                            let sp = surface.point_at(uv.x, uv.y);
                            let d = (0..=32)
                                .map(|i| {
                                    let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 32.0;
                                    ce.point_at(t).distance_to(&sp)
                                })
                                .fold(f64::INFINITY, f64::min);
                            worst = worst.max(d);
                        }
                    }
                    format!("{kind} agree={worst:.2e}")
                }
                None => "NONE".to_string(),
            };
            let a = brep.vertices[data.origin].position;
            let b = brep.vertices[brep.half_edges[data.next].origin].position;
            eprintln!(
                "  strip {fh:?} he {he:?}: {pc_desc} [({:.2},{:.2},{:.2})->({:.2},{:.2},{:.2})]",
                a.x, a.y, a.z, b.x, b.y, b.z
            );
            he = data.next;
            if he == start {
                break;
            }
        }
        break;
    }
}

#[test]
#[ignore = "PK-095 diagnostic: staleness inputs for the 2-CP trim pcurves"]
fn staleness_inputs() {
    let brep = pentagon_blend();
    eprintln!(
        "max_imprecision={:.3e} linear_tol={:.3e}",
        brep.max_imprecision,
        brep.tolerance.linear()
    );
    for (he, data) in brep.half_edges.iter() {
        let Some(pc) = data.pcurve.as_ref() else {
            continue;
        };
        let popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(n) = &pc.geometry else {
            continue;
        };
        if n.control_points.len() != 2 {
            continue;
        }
        let face = brep.loops[data.loop_].face;
        let Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) =
            brep.faces[face].surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        let dev = brep.measure_half_edge_pcurve_deviation(he);
        eprintln!(
            "  2-CP he {he:?}: deviation={:?}",
            dev.map(|r| r.max_deviation)
        );
    }
}

#[test]
#[ignore = "PK-095 diagnostic: does fix_stale_pcurves clear the 2-CP trims"]
fn fix_stale_clears_trim_pcurves() {
    let mut brep = pentagon_blend();
    brep.fix_stale_pcurves();
    let mut remaining = 0;
    let mut worst = 0.0_f64;
    for (he, data) in brep.half_edges.iter() {
        let Some(pc) = data.pcurve.as_ref() else {
            continue;
        };
        let face = brep.loops[data.loop_].face;
        let Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) =
            brep.faces[face].surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        if let popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(n) = &pc.geometry
            && n.control_points.len() == 2
        {
            remaining += 1;
            if let Some(r) = brep.measure_half_edge_pcurve_deviation(he) {
                worst = worst.max(r.max_deviation);
            }
        }
    }
    eprintln!("after explicit fix_stale: {remaining} 2-CP sweep pcurves, worst dev {worst:.3e}");
}

#[test]
#[ignore = "PK-095 diagnostic: 2-CP trim half-edge ground truth"]
fn trim_2cp_ground_truth() {
    let brep = pentagon_blend();
    for (he, data) in brep.half_edges.iter() {
        let Some(pc) = data.pcurve.as_ref() else {
            continue;
        };
        let face = brep.loops[data.loop_].face;
        let Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) =
            brep.faces[face].surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        let popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(n) = &pc.geometry else {
            continue;
        };
        if n.control_points.len() != 2 {
            continue;
        }
        let eh = data.edge;
        let has_curve = brep.edges[eh].curve.is_some();
        let dev = brep
            .measure_half_edge_pcurve_deviation(he)
            .map(|r| r.max_deviation);
        eprintln!(
            "  he {he:?} edge {eh:?} has_curve={has_curve} deviation={dev:?} pc_cps={:?}",
            n.control_points
        );
    }
}
