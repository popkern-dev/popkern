use super::boundary::{ExactBoundaryCurveBudget, PreparedExactBoundaryCurve};
use super::ribbon::{
    ExactBoundaryRibbonBudget, SupportRibbonSampleRequest, compute_support_ribbon_sample,
    prepare_exact_boundary_ribbons,
};
use super::sector::{ExactBoundaryOuterSide, extract_exact_boundary_cap_sectors};
use super::support::ExactBoundarySupportSurfaceFamily;
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{
    NurbsCurve, NurbsSurface, Point3, SurfaceEval, SurfaceGeometry, Vec3,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};
use std::f64::consts::{PI, TAU};

#[test]
fn plane_support_ribbons_use_support_pcurves_and_are_tangent() {
    let (brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &trihedral_sectors(&edges))
            .unwrap();
    let boundaries = super::boundary::prepare_exact_boundary_curves(
        &brep,
        &sectors,
        ExactBoundaryCurveBudget::default(),
    )
    .unwrap();

    let ribbons = prepare_exact_boundary_ribbons(
        &brep,
        &sectors,
        &boundaries,
        ExactBoundaryRibbonBudget::default(),
    )
    .unwrap();

    assert_eq!(ribbons.ribbons.len(), 6);
    for ribbon in &ribbons.ribbons {
        assert!(ribbon.max_support_residual <= 1.0e-6);
        for sample in &ribbon.samples {
            assert!(sample.boundary_point.distance_to(&sample.support_point) <= 1.0e-6);
            assert!(sample.support_normal.dot(&sample.boundary_tangent).abs() <= 1.0e-8);
            assert!(sample.support_normal.dot(&sample.cross_derivative).abs() <= 1.0e-8);
            assert!(sample.boundary_tangent.dot(&sample.cross_derivative).abs() <= 1.0e-8);
        }
    }
}

#[test]
fn missing_support_pcurve_uses_projected_candidate_when_valid() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &trihedral_sectors(&edges))
            .unwrap();
    let missing = sectors[0].support_half_edges[0];
    brep.half_edges[missing].pcurve = None;
    let boundaries = super::boundary::prepare_exact_boundary_curves(
        &brep,
        &sectors,
        ExactBoundaryCurveBudget::default(),
    )
    .unwrap();

    let ribbons = prepare_exact_boundary_ribbons(
        &brep,
        &sectors,
        &boundaries,
        ExactBoundaryRibbonBudget::default(),
    )
    .expect("missing support pcurve should be repaired by validated projection");

    let ribbon = ribbons
        .ribbons
        .iter()
        .find(|ribbon| ribbon.support_half_edge == missing)
        .expect("projected ribbon should preserve support half-edge ownership");
    assert!(ribbon.max_support_residual <= 1.0e-6);
}

#[test]
fn cylinder_ribbon_sample_lies_in_support_tangent_plane() {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let dummy_edge = BRep::unit_cube(Tolerance::default())
        .unwrap()
        .edges
        .keys()
        .next()
        .unwrap();
    let theta = 0.4;
    let boundary = cylinder_line_boundary(&cylinder, theta, 0.2, 0.8, dummy_edge);
    let other_boundary = offset_boundary(&boundary, Vec3::new(0.0, 0.0, 0.2));
    let pcurve = Pcurve::linear((theta, 0.2), (theta, 0.8), 0.0, 1.0);

    let sample = compute_support_ribbon_sample(SupportRibbonSampleRequest {
        boundary: &boundary,
        other_boundary: &other_boundary,
        surface: &cylinder,
        pcurve: &pcurve,
        parameter: 0.5,
        pcurve_reversed: false,
        derivative_scale: 0.2,
        budget: ExactBoundaryRibbonBudget::default(),
    })
    .unwrap();

    assert!(sample.boundary_point.distance_to(&sample.support_point) <= 1.0e-6);
    assert!(sample.support_normal.dot(&sample.boundary_tangent).abs() <= 1.0e-8);
    assert!(sample.support_normal.dot(&sample.cross_derivative).abs() <= 1.0e-8);
    assert!(sample.boundary_tangent.dot(&sample.cross_derivative).abs() <= 1.0e-8);
}

#[test]
fn sphere_cone_and_nurbs_support_samples_are_tangent() {
    let dummy_edge = dummy_edge();
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
    };
    assert_support_sample_is_tangent(
        &sphere,
        Pcurve::linear((0.4, 0.8), (0.4, 0.8001), 0.0, 1.0),
        surface_line_boundary(&sphere, 0.4, 0.8, 0.4, 0.8001, dummy_edge),
        ExactBoundarySupportSurfaceFamily::Sphere,
    );

    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.35,
    };
    assert_support_sample_is_tangent(
        &cone,
        Pcurve::linear((0.3, 0.8), (0.3, 1.4), 0.0, 1.0),
        surface_line_boundary(&cone, 0.3, 0.8, 0.3, 1.4, dummy_edge),
        ExactBoundarySupportSurfaceFamily::Cone,
    );

    let nurbs = SurfaceGeometry::Nurbs(bilinear_nurbs_sheet());
    assert_support_sample_is_tangent(
        &nurbs,
        Pcurve::linear((0.2, 0.3), (0.8, 0.3), 0.0, 1.0),
        surface_line_boundary(&nurbs, 0.2, 0.3, 0.8, 0.3, dummy_edge),
        ExactBoundarySupportSurfaceFamily::Nurbs,
    );
}

#[test]
fn periodic_and_high_curvature_support_samples_remain_stable() {
    let dummy_edge = dummy_edge();
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert_support_sample_is_tangent(
        &cylinder,
        Pcurve::linear((TAU + 0.25, 0.1), (TAU + 0.25, 0.9), 0.0, 1.0),
        surface_line_boundary(&cylinder, TAU + 0.25, 0.1, TAU + 0.25, 0.9, dummy_edge),
        ExactBoundarySupportSurfaceFamily::Cylinder,
    );

    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.08,
    };
    let sample = support_sample(
        &sphere,
        Pcurve::linear((PI / 3.0, 0.9), (PI / 3.0, 0.9001), 0.0, 1.0),
        surface_line_boundary(&sphere, PI / 3.0, 0.9, PI / 3.0, 0.9001, dummy_edge),
        0.01,
    );
    assert!(sample.boundary_point.distance_to(&sample.support_point) <= 1.0e-6);
    assert!(sample.cross_derivative.length() <= 0.0100001);
}

fn hex_cap_fixture() -> (BRep, FaceHandle, Vec<EdgeHandle>) {
    let brep = BRep::n_gonal_pyramid(6, 1.0, 0.8, Tolerance::default()).unwrap();
    let cap_face = brep
        .faces
        .keys()
        .find(|&face| loop_edges(&brep, face).len() == 6)
        .unwrap();
    let edges = loop_edges(&brep, cap_face);
    (brep, cap_face, edges)
}

fn loop_edges(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut cursor = start;
    let mut edges = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = &brep.half_edges[cursor];
        edges.push(half_edge.edge);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    edges
}

fn trihedral_sectors(edges: &[EdgeHandle]) -> Vec<GregoryCapSector> {
    (0..3)
        .map(|sector_index| GregoryCapSector {
            sector_index,
            fillet_index: sector_index,
            fillet_end_edge: edges[2 * sector_index],
            face_edge: edges[2 * sector_index + 1],
        })
        .collect()
}

fn trihedral_plan(brep: &BRep, edges: &[EdgeHandle]) -> ApexPlan {
    ApexPlan {
        apex_vertex: brep.vertices.keys().next().unwrap(),
        cap_kind: CapKind::Gregory {
            per_fillet_terminations: vec![
                [Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)];
                3
            ],
            per_fillet_caps: vec![(Point3::new(0.0, 0.0, 0.0), 1.0); 3],
        },
        target_topology: TargetTopology {
            apex_vertices: Vec::new(),
            apex_edges: (0..6)
                .map(|index| TargetEdge {
                    origin_vertex_idx: index,
                    dest_vertex_idx: (index + 1) % 6,
                    curve: dummy_curve(),
                    role: if index % 2 == 0 {
                        TargetEdgeRole::CapVsCyl {
                            fillet_idx: index / 2,
                        }
                    } else {
                        TargetEdgeRole::CapVsPlane {
                            shared_face_idx: index / 2,
                        }
                    },
                })
                .collect(),
            cap_face: TargetCapFace::Gregory {
                outer_loop_edge_indices: (0..edges.len()).collect(),
            },
            face_boundary_updates: Vec::new(),
        },
        telemetry: PlanTelemetry::default(),
    }
}

fn dummy_curve() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

fn dummy_edge() -> EdgeHandle {
    BRep::unit_cube(Tolerance::default())
        .unwrap()
        .edges
        .keys()
        .next()
        .unwrap()
}

fn assert_support_sample_is_tangent(
    surface: &SurfaceGeometry,
    pcurve: Pcurve,
    boundary: PreparedExactBoundaryCurve,
    family: ExactBoundarySupportSurfaceFamily,
) {
    assert_eq!(super::support::support_surface_family(surface), family);
    let sample = support_sample(surface, pcurve, boundary, 0.2);
    assert!(sample.boundary_point.distance_to(&sample.support_point) <= 1.0e-6);
    assert!(sample.support_normal.dot(&sample.boundary_tangent).abs() <= 1.0e-8);
    assert!(sample.support_normal.dot(&sample.cross_derivative).abs() <= 1.0e-8);
    assert!(sample.boundary_tangent.dot(&sample.cross_derivative).abs() <= 1.0e-8);
}

fn support_sample(
    surface: &SurfaceGeometry,
    pcurve: Pcurve,
    boundary: PreparedExactBoundaryCurve,
    derivative_scale: f64,
) -> super::ribbon::ExactBoundaryRibbonSample {
    let other_boundary = offset_boundary(&boundary, Vec3::new(0.0, 0.0, 0.2));
    compute_support_ribbon_sample(SupportRibbonSampleRequest {
        boundary: &boundary,
        other_boundary: &other_boundary,
        surface,
        pcurve: &pcurve,
        parameter: 0.5,
        pcurve_reversed: false,
        derivative_scale,
        budget: ExactBoundaryRibbonBudget::default(),
    })
    .unwrap()
}

fn surface_line_boundary(
    surface: &SurfaceGeometry,
    u0: f64,
    v0: f64,
    u1: f64,
    v1: f64,
    edge: EdgeHandle,
) -> PreparedExactBoundaryCurve {
    let start = surface.point_at(u0, v0);
    let end = surface.point_at(u1, v1);
    PreparedExactBoundaryCurve {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        side: ExactBoundaryOuterSide::FilletEnd,
        edge,
        curve: NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![start, end],
            weights: vec![1.0, 1.0],
        },
    }
}

fn bilinear_nurbs_sheet() -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    }
}

fn cylinder_line_boundary(
    cylinder: &SurfaceGeometry,
    theta: f64,
    start_height: f64,
    end_height: f64,
    edge: EdgeHandle,
) -> PreparedExactBoundaryCurve {
    let start = cylinder.point_at(theta, start_height);
    let end = cylinder.point_at(theta, end_height);
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    };
    PreparedExactBoundaryCurve {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        side: ExactBoundaryOuterSide::FilletEnd,
        edge,
        curve,
    }
}

fn offset_boundary(
    boundary: &PreparedExactBoundaryCurve,
    offset: Vec3,
) -> PreparedExactBoundaryCurve {
    PreparedExactBoundaryCurve {
        sector_index: boundary.sector_index,
        sector_name: boundary.sector_name.clone(),
        side: ExactBoundaryOuterSide::FaceEdge,
        edge: boundary.edge,
        curve: NurbsCurve {
            degree: boundary.curve.degree,
            knots: boundary.curve.knots.clone(),
            control_points: boundary
                .curve
                .control_points
                .iter()
                .map(|point| *point + offset)
                .collect(),
            weights: boundary.curve.weights.clone(),
        },
    }
}
