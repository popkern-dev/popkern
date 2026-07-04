//! PK-096: per-edge curve-vs-adjacent-surface deviation audit — what a
//! STEP consumer's stitching tolerance actually sees.
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
#[ignore = "PK-096 diagnostic: edge-surface stitch gaps"]
fn edge_surface_gap_audit() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.partial_cmp(&b.1.position.z).unwrap())
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let o = pyramid.half_edges[he].origin;
            let t = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (o == apex || t == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<popkern_fillet::FilletSpec> = apex_edges
        .iter()
        .map(|&edge| popkern_fillet::FilletSpec {
            edge,
            radius: popkern_fillet::FilletRadius::Constant(0.16),
        })
        .collect();
    let options = popkern_fillet::FilletOptions {
        continuity: popkern_fillet::ContinuityLevel::G2,
        ..Default::default()
    };
    let result =
        popkern_fillet::fillet_edges(&pyramid, &specs, &options, tol, None).expect("blend");
    let brep = &result.brep;

    // For each edge: max over samples of distance(curve(t), nearest
    // point on each adjacent surface).
    let mut rows: Vec<(f64, String)> = Vec::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            rows.push((f64::NAN, format!("{eh:?}: NO CURVE")));
            continue;
        };
        let ceval: &dyn CurveEval = &curve.geometry;
        let he = edge.half_edge;
        let faces = [
            brep.loops[brep.half_edges[he].loop_].face,
            brep.loops[brep.half_edges[brep.half_edges[he].twin].loop_].face,
        ];
        let mut worst = 0.0_f64;
        let mut kinds = Vec::new();
        for fh in faces {
            let Some(sof) = brep.faces.get(fh).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            kinds.push(match &sof.geometry {
                popkern_core::geometry::SurfaceGeometry::Plane { .. } => "P",
                popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "N",
                popkern_core::geometry::SurfaceGeometry::Gregory(_) => "G",
                _ => "?",
            });
            let mut seed: Option<(f64, f64)> = None;
            for i in 0..=16 {
                let t = curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 16.0;
                let p = ceval.point_at(t);
                let uv = match seed {
                    Some((u, v)) => project_point_to_uv_near(&sof.geometry, &p, u, v)
                        .or_else(|| project_point_to_uv(&sof.geometry, &p)),
                    None => project_point_to_uv(&sof.geometry, &p),
                };
                if let Some((u, v)) = uv {
                    seed = Some((u, v));
                    let d = sof.geometry.point_at(u, v).distance_to(&p);
                    worst = worst.max(d);
                }
            }
        }
        rows.push((worst, format!("{eh:?} [{}]", kinds.join(""))));
    }
    rows.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    for (d, label) in rows.iter().take(20) {
        eprintln!("{d:.3e}  {label}");
    }
    let over: usize = rows.iter().filter(|(d, _)| *d > 1.0e-5).count();
    eprintln!("edges over 1e-5: {over} / {}", rows.len());

    // Pcurve fidelity: what a STEP consumer's SURFACE_CURVE check
    // sees — max deviation between each half-edge's pcurve image on
    // its surface and the edge's 3D curve.
    let mut worst_pc: Vec<(f64, String)> = Vec::new();
    for (he, data) in brep.half_edges.iter() {
        if let Some(r) = brep.measure_half_edge_pcurve_deviation(he) {
            let face = brep.loops[data.loop_].face;
            let kind = match brep.faces[face].surface.as_ref().map(|s| &s.geometry) {
                Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "P",
                Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "N",
                Some(popkern_core::geometry::SurfaceGeometry::Gregory(_)) => "G",
                _ => "?",
            };
            let pc_kind = match data.pcurve.as_ref().map(|p| &p.geometry) {
                Some(popkern_core::geometry::pcurve::PcurveGeometry::Linear { .. }) => {
                    "Lin".to_string()
                }
                Some(popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(n)) => {
                    format!("N{}", n.control_points.len())
                }
                Some(_) => "Other".to_string(),
                None => "None".to_string(),
            };
            worst_pc.push((r.max_deviation, format!("{he:?} on {kind} pc={pc_kind}")));
        }
    }
    worst_pc.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    eprintln!("--- pcurve deviations (top 8) ---");
    for (d, label) in worst_pc.iter().take(8) {
        eprintln!("{d:.3e}  {label}");
    }
    let over_pc: usize = worst_pc.iter().filter(|(d, _)| *d > 1.0e-5).count();
    eprintln!("half-edges over 1e-5: {over_pc} / {}", worst_pc.len());

    // Direct check: sample the cap faces' grid rim rows against the
    // strip sweeps.
    let sweeps: Vec<_> = brep
        .faces
        .iter()
        .filter_map(|(_, f)| {
            let sof = f.surface.as_ref()?;
            matches!(
                sof.geometry,
                popkern_core::geometry::SurfaceGeometry::Nurbs(_)
            )
            .then(|| sof.geometry.clone())
        })
        .collect();
    let mut worst_row = 0.0_f64;
    for (_, f) in brep.faces.iter() {
        let Some(popkern_core::geometry::SurfaceGeometry::Gregory(
            popkern_core::geometry::gregory::GregoryPatch::RectangularBiNonic(bp),
        )) = f.surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        for side in 0..4 {
            let row: Vec<popkern_core::geometry::Point3> = (0..10)
                .map(|k| match side {
                    0 => bp.grid[0][k],
                    1 => bp.grid[9][k],
                    2 => bp.grid[k][0],
                    _ => bp.grid[k][9],
                })
                .collect();
            // Is this a rim side? every CP within 2e-3 of a sweep.
            let near = |p: &popkern_core::geometry::Point3| -> Option<f64> {
                let mut best = f64::INFINITY;
                for sw in &sweeps {
                    if let Some((u, v)) = project_point_to_uv(sw, p) {
                        best = best.min(sw.point_at(u, v).distance_to(p));
                    }
                }
                (best < 2.0e-3).then_some(best)
            };
            if row.iter().all(|p| near(p).is_some()) {
                // Sample the CURVE, not the CPs.
                let mut side_worst = 0.0_f64;
                for i in 0..=16 {
                    let t = i as f64 / 16.0;
                    let mut work = row.clone();
                    for level in (1..10).rev() {
                        for j in 0..level {
                            work[j] = popkern_core::geometry::Point3::new(
                                work[j].x + (work[j + 1].x - work[j].x) * t,
                                work[j].y + (work[j + 1].y - work[j].y) * t,
                                work[j].z + (work[j + 1].z - work[j].z) * t,
                            );
                        }
                    }
                    if let Some(d) = near(&work[0]) {
                        side_worst = side_worst.max(d);
                    }
                }
                if side_worst > 1.0e-6 {
                    eprintln!(
                        "  row side {side} worst {side_worst:.3e} start ({:.3},{:.3},{:.3})",
                        row[0].x, row[0].y, row[0].z
                    );
                }
                worst_row = worst_row.max(side_worst);
            }
        }
    }
    eprintln!("worst reconciled rim row vs sweep: {worst_row:.3e}");
}

#[test]
#[ignore = "PK-096 diagnostic: Gregory closest_point convergence"]
fn gregory_projection_convergence() {
    use popkern_core::geometry::surface::SurfaceEval;
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.partial_cmp(&b.1.position.z).unwrap())
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let o = pyramid.half_edges[he].origin;
            let t = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (o == apex || t == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<popkern_fillet::FilletSpec> = apex_edges
        .iter()
        .map(|&edge| popkern_fillet::FilletSpec {
            edge,
            radius: popkern_fillet::FilletRadius::Constant(0.16),
        })
        .collect();
    let options = popkern_fillet::FilletOptions {
        continuity: popkern_fillet::ContinuityLevel::G2,
        ..Default::default()
    };
    let result =
        popkern_fillet::fillet_edges(&pyramid, &specs, &options, tol, None).expect("blend");
    for (_, f) in result.brep.faces.iter() {
        let Some(popkern_core::geometry::SurfaceGeometry::Gregory(_)) =
            f.surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        let sfc = &f.surface.as_ref().unwrap().geometry;
        // Points ON the surface at various params, incl. boundaries.
        for (u, v) in [(0.3, 0.0), (0.7, 0.0), (0.0, 0.4), (0.5, 0.5), (1.0, 0.6)] {
            let p = sfc.point_at(u, v);
            match project_point_to_uv(sfc, &p) {
                Some((pu, pv)) => {
                    let q = sfc.point_at(pu, pv);
                    eprintln!(
                        "  ({u},{v}) -> ({pu:.4},{pv:.4}) roundtrip {:.3e}",
                        q.distance_to(&p)
                    );
                }
                None => eprintln!("  ({u},{v}) -> FAILED"),
            }
        }
        break;
    }
}

#[test]
#[ignore = "PK-096 diagnostic: round-trip the exported STEP and audit what a consumer sees"]
fn step_roundtrip_gap_audit() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.partial_cmp(&b.1.position.z).unwrap())
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let o = pyramid.half_edges[he].origin;
            let t = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (o == apex || t == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<popkern_fillet::FilletSpec> = apex_edges
        .iter()
        .map(|&edge| popkern_fillet::FilletSpec {
            edge,
            radius: popkern_fillet::FilletRadius::Constant(0.16),
        })
        .collect();
    let options = popkern_fillet::FilletOptions {
        continuity: popkern_fillet::ContinuityLevel::G2,
        ..Default::default()
    };
    let result =
        popkern_fillet::fillet_edges(&pyramid, &specs, &options, tol, None).expect("blend");

    // Export exactly as the showcase does.
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: "roundtrip".into(),
            id: "RT-001".into(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids: vec![result.brep.clone()],
                context: popkern_io::RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata::default(),
        }],
        roots: vec![0],
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    let data =
        popkern_io::export_step(&doc, &popkern_io::StepExportOptions::default()).expect("export");
    let imported =
        popkern_io::import_step(&data, &popkern_io::StepImportOptions::default()).expect("import");
    let mut audited = 0;
    for product in &imported.document.products {
        for shape in &product.shapes {
            for brep in &shape.solids {
                audited += 1;
                eprintln!(
                    "imported: faces={} edges={} vertices={} validate={:?}",
                    brep.faces.len(),
                    brep.edges.len(),
                    brep.vertices.len(),
                    brep.validate().err().map(|e| e.to_string())
                );
                // Per-edge curve-vs-adjacent-surface gaps on the
                // imported model.
                let mut rows: Vec<(f64, String)> = Vec::new();
                for (eh, edge) in brep.edges.iter() {
                    let Some(curve) = edge.curve.as_ref() else {
                        rows.push((f64::NAN, format!("{eh:?}: NO CURVE")));
                        continue;
                    };
                    let ceval: &dyn CurveEval = &curve.geometry;
                    let he = edge.half_edge;
                    let Some(he_data) = brep.half_edges.get(he) else {
                        rows.push((f64::NAN, format!("{eh:?}: STALE half_edge")));
                        continue;
                    };
                    let Some(twin) = brep.half_edges.get(he_data.twin) else {
                        rows.push((f64::NAN, format!("{eh:?}: STALE twin")));
                        continue;
                    };
                    let f1 = brep.loops.get(he_data.loop_).map(|l| l.face);
                    let f2 = brep.loops.get(twin.loop_).map(|l| l.face);
                    let faces: Vec<_> = [f1, f2].into_iter().flatten().collect();
                    let mut worst = 0.0_f64;
                    for fh in faces {
                        let Some(sof) = brep.faces.get(fh).and_then(|f| f.surface.as_ref()) else {
                            continue;
                        };
                        let mut seed: Option<(f64, f64)> = None;
                        for i in 0..=16 {
                            let t = curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 16.0;
                            let p = ceval.point_at(t);
                            let uv = match seed {
                                Some((u, v)) => project_point_to_uv_near(&sof.geometry, &p, u, v)
                                    .or_else(|| project_point_to_uv(&sof.geometry, &p)),
                                None => project_point_to_uv(&sof.geometry, &p),
                            };
                            if let Some((u, v)) = uv {
                                seed = Some((u, v));
                                let d = sof.geometry.point_at(u, v).distance_to(&p);
                                worst = worst.max(d);
                            }
                        }
                    }
                    rows.push((worst, format!("{eh:?}")));
                }
                rows.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
                for (d, label) in rows.iter().take(10) {
                    eprintln!("  {d:.3e}  {label}");
                }
                let over: usize = rows.iter().filter(|(d, _)| *d > 1.0e-5).count();
                eprintln!("  imported edges over 1e-5: {over} / {}", rows.len());
            }
        }
    }
    assert!(audited > 0, "no solids imported");
}
