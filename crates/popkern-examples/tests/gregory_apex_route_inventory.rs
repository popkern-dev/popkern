use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, BrepEdgeSurfaceCompatibilityOptions, BrepTopologyOutputPolicy, BrepTrimRepairOptions,
    analyze_brep_edge_surface_compatibility, analyze_brep_export_quality,
};
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

#[test]
#[ignore = "A2 migration promotion gate: native route is topologically complete but not STEP-fidelity equivalent yet"]
fn mixed_radius_gregory_a2_quality_gap_is_explicit() {
    let legacy = mixed_radius_gregory_cube(FilletOptions::default());
    let a2 = mixed_radius_gregory_cube(FilletOptions {
        use_a2_apex_construction: true,
        ..FilletOptions::default()
    });
    let legacy_quality = quality(&legacy.brep);
    let a2_quality = quality(&a2.brep);

    assert_eq!(
        a2.repair_report.apex_routes[0].route,
        ApexConstructionRoute::A2
    );
    assert!(
        a2_quality.pcurve_blockers <= legacy_quality.pcurve_blockers
            && a2_quality.edge_surface_blockers <= legacy_quality.edge_surface_blockers
            && a2_quality.max_pcurve_deviation <= legacy_quality.max_pcurve_deviation
            && a2_quality.max_edge_surface_deviation <= legacy_quality.max_edge_surface_deviation,
        "native A2 is not STEP-fidelity promotable yet: legacy={legacy_quality:?}, a2={a2_quality:?}"
    );
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

#[derive(Debug)]
struct RouteQuality {
    pcurve_blockers: usize,
    max_pcurve_deviation: f64,
    edge_surface_blockers: usize,
    max_edge_surface_deviation: f64,
}

fn quality(brep: &BRep) -> RouteQuality {
    let export = analyze_brep_export_quality(
        brep,
        popkern_core::topology::BrepExportQualityOptions {
            topology_output: BrepTopologyOutputPolicy::OpenSheet,
            trim: BrepTrimRepairOptions {
                endpoint_tolerance: 1e-6,
                sampled_tolerance: 1e-6,
                samples: 9,
                validate_output: true,
            },
            require_edge_curves: true,
            allow_collapsed_trim_edge_curves: true,
            max_pcurve_deviation: Some(1e-6),
            max_edge_surface_deviation: Some(1e-6),
            validate_brep: true,
        },
    );
    let edge_surface = analyze_brep_edge_surface_compatibility(
        brep,
        BrepEdgeSurfaceCompatibilityOptions {
            tolerance: 1e-6,
            max_anchor_motion: 0.0,
            samples: 9,
            require_all_compatible: false,
            validate_output: false,
        },
    );
    RouteQuality {
        pcurve_blockers: export.pcurve_blockers,
        max_pcurve_deviation: export.max_pcurve_deviation,
        edge_surface_blockers: edge_surface.incompatible_half_edges,
        max_edge_surface_deviation: edge_surface.max_deviation,
    }
}
