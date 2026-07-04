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
    for (he, _) in brep.half_edges.iter() {
        if let Some(r) = brep.measure_half_edge_pcurve_deviation(he) {
            worst_pc.push((r.max_deviation, format!("{he:?}")));
        }
    }
    worst_pc.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    eprintln!("--- pcurve deviations (top 8) ---");
    for (d, label) in worst_pc.iter().take(8) {
        eprintln!("{d:.3e}  {label}");
    }
    let over_pc: usize = worst_pc.iter().filter(|(d, _)| *d > 1.0e-5).count();
    eprintln!("half-edges over 1e-5: {over_pc} / {}", worst_pc.len());
}
