use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionRoute, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

#[test]
fn mixed_radius_gregory_records_default_legacy_route() {
    let result = mixed_radius_gregory_cube(FilletOptions::default());
    let routes = &result.repair_report.apex_routes;

    assert_eq!(routes.len(), 1, "fixture should have one smooth apex route");
    assert_eq!(routes[0].route, ApexConstructionRoute::Legacy);
    assert_eq!(routes[0].incident_edges.len(), 3);
    assert!(routes[0].built_faces.is_empty());
    assert!(routes[0].fallback_reason.is_none());
}

#[test]
fn mixed_radius_gregory_records_a2_requested_route() {
    let result = mixed_radius_gregory_cube(FilletOptions {
        use_a2_apex_construction: true,
        ..FilletOptions::default()
    });
    let routes = &result.repair_report.apex_routes;

    assert_eq!(routes.len(), 1, "fixture should have one smooth apex route");
    assert_eq!(routes[0].route, ApexConstructionRoute::A2);
    assert_eq!(routes[0].incident_edges.len(), 3);
    assert_eq!(routes[0].built_faces.len(), 3);
    assert!(routes[0].fallback_reason.is_none());
}

fn mixed_radius_gregory_cube(options: FilletOptions) -> popkern_fillet::BlendResult {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges = cube
        .half_edges
        .values()
        .filter(|half_edge| half_edge.origin == vertex)
        .map(|half_edge| half_edge.edge)
        .collect::<Vec<_>>();
    let radii = [0.08, 0.10, 0.12];
    let specs = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect::<Vec<_>>();
    fillet_edges(&cube, &specs, &options, Tolerance::default(), None)
        .expect("mixed-radius Gregory cube must build")
}
