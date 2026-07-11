//! Audit exactly what the end-runoff pass fires on.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
#[ignore = "PK-095 diagnostic: runoff firing audit"]
fn runoff_firing_audit() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
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
    let mut n = 0;
    for d in &result.repair_report.apex_diagnostics {
        if let popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(e) = d
            && e.stage == "end_runoff"
        {
            n += 1;
            eprintln!("runoff: {} overshoot={:?}", e.sector_name, e.residual);
        }
    }
    eprintln!("total runoff firings: {n}");
    eprintln!(
        "faces={} edges={} vertices={}",
        result.brep.faces.len(),
        result.brep.edges.len(),
        result.brep.vertices.len()
    );
}
