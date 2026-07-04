//! PK-097 Phase A gate: the PnS bi-4 cap assembled on the SETBACK
//! border of a G1 pentagon-pyramid apex.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
#[ignore = "PK-097 Phase A: measurement + gate for the G1 setback candidate"]
fn pns_setback_candidate_on_pentagon_g1() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
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
