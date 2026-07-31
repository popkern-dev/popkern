//! PK-097 Phase A gate: the PnS bi-4 cap assembled on the SETBACK
//! border of a G1 pentagon-pyramid apex.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn pns_setback_candidate_on_pentagon_g1() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.partial_cmp(&b.1.position.z).unwrap())
        .map(|(h, v)| (h, v.position))
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let o = pyramid.half_edges[he].origin;
            let t = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (o == apex.0 || t == apex.0).then_some(eh)
        })
        .collect();
    let specs: Vec<popkern_fillet::FilletSpec> = apex_edges
        .iter()
        .map(|&edge| popkern_fillet::FilletSpec {
            edge,
            radius: popkern_fillet::FilletRadius::Constant(0.16),
        })
        .collect();
    // G1 sweeps — the arc cross-section path.
    let preview = popkern_fillet::preview_fillet_with_continuity(
        &pyramid,
        &specs,
        tol,
        popkern_fillet::ContinuityLevel::G1,
    )
    .expect("G1 preview");
    let fillet_surfaces: Vec<popkern_fillet::fillet_surface::FilletSurfaceResult> = preview
        .per_edge
        .iter()
        .map(|eg| popkern_fillet::fillet_surface::FilletSurfaceResult {
            surface: eg.blend_surface.clone(),
            trim_curve_1: eg.trim_curve_1.clone(),
            trim_curve_2: eg.trim_curve_2.clone(),
            continuity_downgrade: None,
        })
        .collect();
    let report =
        popkern_fillet::construct_test_hooks::pns_setback_candidate(&fillet_surfaces, apex.1)
            .expect("candidate");
    eprintln!("{report:?}");
    // Radial seam check on the SURFACES (not just control polygons):
    // sample each adjacent pair along their shared side.
    {
        use popkern_core::geometry::surface::SurfaceEval;
        use popkern_core::geometry::uv_projection::project_point_to_uv;
        let full = popkern_fillet::construct_test_hooks::pns_setback_candidate_full(
            &fillet_surfaces,
            apex.1,
        )
        .expect("candidate");
        let np = full.surfaces.len();
        let mut worst = 0.0_f64;
        for k in 0..np {
            let a = &full.surfaces[k];
            let b = &full.surfaces[(k + 1) % np];
            let mut best_side_worst = f64::INFINITY;
            for side in 0..4 {
                let mut side_worst = 0.0_f64;
                for i in 0..=8 {
                    let t = i as f64 / 8.0;
                    let (u, v) = match side {
                        0 => (0.0, t),
                        1 => (1.0, t),
                        2 => (t, 0.0),
                        _ => (t, 1.0),
                    };
                    let pnt = a.point_at(u, v);
                    let d = match project_point_to_uv(b, &pnt) {
                        Some((pu, pv)) => b.point_at(pu, pv).distance_to(&pnt),
                        None => f64::INFINITY,
                    };
                    side_worst = side_worst.max(d);
                }
                best_side_worst = best_side_worst.min(side_worst);
            }
            worst = worst.max(best_side_worst);
        }
        eprintln!("surface radial seam worst: {worst:.4e}");
        for k in 0..2 {
            let p2 = &full.patches[k];
            eprintln!(
                "patch {k} corners [0][0]=({:.3},{:.3},{:.3}) [0][4]=({:.3},{:.3},{:.3}) [4][0]=({:.3},{:.3},{:.3}) [4][4]=({:.3},{:.3},{:.3}) pole=({:.3},{:.3},{:.3})",
                p2[0][0].x,
                p2[0][0].y,
                p2[0][0].z,
                p2[0][4].x,
                p2[0][4].y,
                p2[0][4].z,
                p2[4][0].x,
                p2[4][0].y,
                p2[4][0].z,
                p2[4][4].x,
                p2[4][4].y,
                p2[4][4].z,
                full.pole.x,
                full.pole.y,
                full.pole.z
            );
        }
        // Side layout matrix for patches 0..2: each side's midpoint
        // distance to (prev patch, next patch).
        for k in 0..3 {
            let np2 = full.surfaces.len();
            let prev = &full.surfaces[(k + np2 - 1) % np2];
            let next = &full.surfaces[(k + 1) % np2];
            let patch = &full.patches[k];
            for side in 0..4 {
                let side_cps: [popkern_core::geometry::Point3; 5] = match side {
                    0 => core::array::from_fn(|j| patch[0][j]),
                    1 => core::array::from_fn(|j| patch[4][j]),
                    2 => core::array::from_fn(|i| patch[i][0]),
                    _ => core::array::from_fn(|i| patch[i][4]),
                };
                let mid = {
                    let mut w = side_cps;
                    for level in (1..5).rev() {
                        for i in 0..level {
                            w[i] = popkern_core::geometry::Point3::new(
                                w[i].x + (w[i + 1].x - w[i].x) * 0.5,
                                w[i].y + (w[i + 1].y - w[i].y) * 0.5,
                                w[i].z + (w[i + 1].z - w[i].z) * 0.5,
                            );
                        }
                    }
                    w[0]
                };
                let dist = |sfc: &popkern_core::geometry::SurfaceGeometry| match project_point_to_uv(
                    sfc, &mid,
                ) {
                    Some((u, v)) => sfc.point_at(u, v).distance_to(&mid),
                    None => f64::NAN,
                };
                // Distance of the side mid to the nearest half-arc.
                let mut arc_d = f64::INFINITY;
                for (cubic, _) in &full.half_cubics {
                    for i in 0..=32 {
                        let t = i as f64 / 32.0;
                        let mut w = *cubic;
                        for level in (1..4).rev() {
                            for j in 0..level {
                                w[j] = popkern_core::geometry::Point3::new(
                                    w[j].x + (w[j + 1].x - w[j].x) * t,
                                    w[j].y + (w[j + 1].y - w[j].y) * t,
                                    w[j].z + (w[j + 1].z - w[j].z) * t,
                                );
                            }
                        }
                        arc_d = arc_d.min(w[0].distance_to(&mid));
                    }
                }
                eprintln!(
                    "patch {k} side {side}: prev {:.2e} next {:.2e} arc {:.2e}",
                    dist(prev),
                    dist(next),
                    arc_d
                );
            }
        }
    }
    // Phase A gates (measured on the pentagon at radius 0.16):
    // border exact by construction, seams exact by eq (2) + shared
    // deltas, G1 within the production PnS coplanarity gate (5e-2)
    // with margin.
    assert!(
        report.g0_border_after < 1.0e-12,
        "border must equal the setback arcs: {:.3e}",
        report.g0_border_after
    );
    assert!(
        report.seam_g0 < 1.0e-12,
        "radial seams must be exact: {:.3e}",
        report.seam_g0
    );
    assert!(
        report.g1_angle_after < 5.0e-2,
        "border G1 must beat the production coplanarity gate: {:.3e}",
        report.g1_angle_after
    );
}

/// PK-200 acceptance: the setback candidate emits first-class
/// PolyhedralSpline sectors with the cap's true valence (2n — this
/// ring splits every arc at its midpoint) and dense sector indices.
#[test]
fn pns_setback_candidate_emits_provenance() {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.partial_cmp(&b.1.position.z).unwrap())
        .map(|(h, v)| (h, v.position))
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let o = pyramid.half_edges[he].origin;
            let t = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (o == apex.0 || t == apex.0).then_some(eh)
        })
        .collect();
    let specs: Vec<popkern_fillet::FilletSpec> = apex_edges
        .iter()
        .map(|&edge| popkern_fillet::FilletSpec {
            edge,
            radius: popkern_fillet::FilletRadius::Constant(0.16),
        })
        .collect();
    let preview = popkern_fillet::preview_fillet_with_continuity(
        &pyramid,
        &specs,
        tol,
        popkern_fillet::ContinuityLevel::G1,
    )
    .expect("G1 preview");
    let fillet_surfaces: Vec<popkern_fillet::fillet_surface::FilletSurfaceResult> = preview
        .per_edge
        .iter()
        .map(|eg| popkern_fillet::fillet_surface::FilletSurfaceResult {
            surface: eg.blend_surface.clone(),
            trim_curve_1: eg.trim_curve_1.clone(),
            trim_curve_2: eg.trim_curve_2.clone(),
            continuity_downgrade: None,
        })
        .collect();
    let full =
        popkern_fillet::construct_test_hooks::pns_setback_candidate_full(&fillet_surfaces, apex.1)
            .expect("candidate");
    let count = full.surfaces.len();
    assert_eq!(count, full.patches.len());
    assert_eq!(count, 2 * 5, "split ring: two sectors per arc");
    for (k, surface) in full.surfaces.iter().enumerate() {
        match surface {
            SurfaceGeometry::PolyhedralSpline(patch) => {
                assert_eq!(patch.provenance.valence, count, "sector {k} valence");
                assert_eq!(patch.provenance.sector, k, "sector {k} index");
            }
            other => panic!("sector {k}: expected PolyhedralSpline, got {other:?}"),
        }
    }
}
