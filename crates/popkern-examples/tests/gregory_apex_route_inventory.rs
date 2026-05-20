use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, BrepEdgeSurfaceCompatibilityOptions, BrepTopologyOutputPolicy, BrepTrimRepairOptions,
    analyze_brep_edge_surface_compatibility, analyze_brep_export_quality,
};
use popkern_fillet::{
    ApexConstructionRoute, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

#[test]
fn mixed_radius_gregory_records_default_direct_apex_route() {
    let result = mixed_radius_gregory_cube(FilletOptions::default());
    let routes = &result.repair_report.apex_routes;

    assert_eq!(routes.len(), 1, "fixture should have one smooth apex route");
    assert_eq!(routes[0].route, ApexConstructionRoute::DirectApex);
    assert_eq!(routes[0].incident_edges.len(), 3);
    assert_eq!(routes[0].built_faces.len(), 3);
    assert_eq!(routes[0].exact_boundary_sector_count, 3);
    assert_eq!(routes[0].exact_boundary_prepared_curve_count, 6);
    assert!(routes[0].fallback_reason.is_none());
}

#[test]
fn mixed_radius_gregory_records_direct_apex_requested_route() {
    let result = mixed_radius_gregory_cube(FilletOptions {
        use_direct_apex_construction: true,
        ..FilletOptions::default()
    });
    let routes = &result.repair_report.apex_routes;

    assert_eq!(routes.len(), 1, "fixture should have one smooth apex route");
    assert_eq!(routes[0].route, ApexConstructionRoute::DirectApex);
    assert_eq!(routes[0].incident_edges.len(), 3);
    assert_eq!(routes[0].built_faces.len(), 3);
    assert_eq!(routes[0].exact_boundary_sector_count, 3);
    assert_eq!(routes[0].exact_boundary_prepared_curve_count, 6);
    assert!(routes[0].fallback_reason.is_none());
}

#[test]
fn mixed_radius_gregory_default_matches_requested_direct_apex_quality() {
    let default = mixed_radius_gregory_cube(FilletOptions::default());
    let requested_direct_apex = mixed_radius_gregory_cube(FilletOptions {
        use_direct_apex_construction: true,
        ..FilletOptions::default()
    });
    let default_quality = quality(&default.brep);
    let requested_direct_apex_quality = quality(&requested_direct_apex.brep);
    let _evidence_rows =
        default_quality.evidence_rows() + requested_direct_apex_quality.evidence_rows();

    assert_eq!(
        default.repair_report.apex_routes[0].route,
        ApexConstructionRoute::DirectApex
    );
    assert_eq!(
        requested_direct_apex.repair_report.apex_routes[0].route,
        ApexConstructionRoute::DirectApex
    );
    assert!(
        default_quality.pcurve_blockers == requested_direct_apex_quality.pcurve_blockers
            && default_quality.edge_surface_blockers
                == requested_direct_apex_quality.edge_surface_blockers
            && default_quality.max_pcurve_deviation
                == requested_direct_apex_quality.max_pcurve_deviation
            && default_quality.max_edge_surface_deviation
                == requested_direct_apex_quality.max_edge_surface_deviation,
        "default mixed-radius Gregory should use the same STEP-fidelity path as requested DirectApex: default={default_quality:?}, requested_direct_apex={requested_direct_apex_quality:?}"
    );
}

impl RouteQuality {
    fn evidence_rows(&self) -> usize {
        let pcurve_rows = self
            .worst_pcurve_edges
            .iter()
            .map(|edge| {
                edge.edge.len()
                    + (edge.max_deviation > 0.0) as usize
                    + edge
                        .half_edges
                        .iter()
                        .map(|half_edge| {
                            half_edge.half_edge.len()
                                + half_edge.face.len()
                                + (half_edge.max_deviation > 0.0) as usize
                        })
                        .sum::<usize>()
            })
            .sum::<usize>();
        let edge_surface_rows = self
            .worst_edge_surface_edges
            .iter()
            .map(|edge| {
                edge.edge.len()
                    + edge.half_edge.len()
                    + edge.face.len()
                    + (edge.max_deviation > 0.0) as usize
            })
            .sum::<usize>();
        pcurve_rows + edge_surface_rows
    }
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
    worst_pcurve_edges: Vec<EdgePcurveEvidence>,
    worst_edge_surface_edges: Vec<EdgeSurfaceEvidence>,
}

#[derive(Debug)]
struct EdgePcurveEvidence {
    edge: String,
    max_deviation: f64,
    half_edges: Vec<HalfEdgePcurveEvidence>,
}

#[derive(Debug)]
struct HalfEdgePcurveEvidence {
    half_edge: String,
    face: String,
    max_deviation: f64,
}

#[derive(Debug)]
struct EdgeSurfaceEvidence {
    edge: String,
    half_edge: String,
    face: String,
    max_deviation: f64,
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
        worst_pcurve_edges: worst_pcurve_edges(brep, 30),
        worst_edge_surface_edges: edge_surface
            .findings
            .iter()
            .take(30)
            .map(|finding| EdgeSurfaceEvidence {
                edge: format!("{:?}", finding.edge),
                half_edge: format!("{:?}", finding.half_edge),
                face: format!("{:?}", finding.face),
                max_deviation: finding.max_deviation,
            })
            .collect(),
    }
}

fn worst_pcurve_edges(brep: &BRep, limit: usize) -> Vec<EdgePcurveEvidence> {
    let mut reports = brep
        .edges
        .keys()
        .filter_map(|edge| brep.measure_edge_pcurve_deviation(edge))
        .filter(|report| report.max_deviation > 1e-12)
        .collect::<Vec<_>>();
    reports.sort_by(|a, b| b.max_deviation.total_cmp(&a.max_deviation));
    reports
        .into_iter()
        .take(limit)
        .map(|report| EdgePcurveEvidence {
            edge: format!("{:?}", report.edge),
            max_deviation: report.max_deviation,
            half_edges: report
                .half_edges
                .iter()
                .map(|half_edge| HalfEdgePcurveEvidence {
                    half_edge: format!("{:?}", half_edge.half_edge),
                    face: format!("{:?}", half_edge.face),
                    max_deviation: half_edge.max_deviation,
                })
                .collect(),
        })
        .collect()
}
