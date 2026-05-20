use super::boundary::{ExactBoundaryCurveBudget, PreparedExactBoundaryCurve};
use super::ribbon::{
    ExactBoundaryRibbonBudget, ExactBoundaryRibbonDeclineReason, compute_support_ribbon_sample,
    prepare_exact_boundary_ribbons,
};
use super::sector::{ExactBoundaryOuterSide, extract_exact_boundary_cap_sectors};
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{NurbsCurve, Point3, SurfaceEval, SurfaceGeometry, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

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
fn missing_support_pcurve_declines_with_stable_diagnostic() {
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

    let diagnostic = prepare_exact_boundary_ribbons(
        &brep,
        &sectors,
        &boundaries,
        ExactBoundaryRibbonBudget::default(),
    )
    .expect_err("missing support pcurve must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryRibbonDeclineReason::SupportPcurveMissing
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.ribbon.support_pcurve_missing"
    );
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

    let sample = compute_support_ribbon_sample(
        &boundary,
        &other_boundary,
        &cylinder,
        &pcurve,
        0.5,
        false,
        0.2,
        ExactBoundaryRibbonBudget::default(),
    )
    .unwrap();

    assert!(sample.boundary_point.distance_to(&sample.support_point) <= 1.0e-6);
    assert!(sample.support_normal.dot(&sample.boundary_tangent).abs() <= 1.0e-8);
    assert!(sample.support_normal.dot(&sample.cross_derivative).abs() <= 1.0e-8);
    assert!(sample.boundary_tangent.dot(&sample.cross_derivative).abs() <= 1.0e-8);
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
