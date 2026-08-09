use super::boundary::{ExactBoundaryCurveBudget, prepare_exact_boundary_curves};
use super::corner::{ExactBoundaryInnerSpokeBudget, prepare_exact_boundary_inner_spokes_with_pole};
use super::ribbon::{ExactBoundaryRibbonBudget, prepare_exact_boundary_ribbons};
use super::sector::extract_exact_boundary_cap_sectors;
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

#[test]
fn support_pcurve_adversarial_matrix_has_stable_outcomes() {
    assert_support_pcurve_is_repaired(SupportPcurveAdversary::Missing);
    assert_support_pcurve_is_repaired(SupportPcurveAdversary::Stale);
    assert_reversed_support_pcurve_is_accepted();
}

#[test]
fn geometry_budget_adversarial_matrix_has_stable_declines() {
    let (brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &trihedral_sectors(&edges))
            .unwrap();

    let boundary_diagnostic = prepare_exact_boundary_curves(
        &brep,
        &sectors,
        ExactBoundaryCurveBudget {
            target_degree: 5,
            max_control_points: 2,
            max_knots: 4,
        },
    )
    .expect_err("over-constrained boundary budget must decline");
    assert_eq!(
        boundary_diagnostic.code,
        "exact_boundary_cap.knot_budget_exceeded"
    );

    let boundaries =
        prepare_exact_boundary_curves(&brep, &sectors, ExactBoundaryCurveBudget::default())
            .unwrap();
    let spoke_diagnostic = prepare_exact_boundary_inner_spokes_with_pole(
        &sectors,
        &boundaries,
        ExactBoundaryInnerSpokeBudget {
            min_spoke_length: 10.0,
            ..ExactBoundaryInnerSpokeBudget::default()
        },
        None,
    )
    .expect_err("collapsed generated spoke must decline");
    assert_eq!(
        spoke_diagnostic.code,
        "exact_boundary_cap.inner_spoke.collapsed"
    );
}

#[test]
fn rational_imported_boundary_is_prepared_without_losing_weights() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let rational_edge = edges[0];
    install_rational_quadratic_edge_curve(&mut brep, rational_edge);
    let plan = trihedral_plan(&brep, &edges);
    let sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &trihedral_sectors(&edges))
            .unwrap();

    let prepared =
        prepare_exact_boundary_curves(&brep, &sectors, ExactBoundaryCurveBudget::default())
            .expect("rational imported boundary should remain an exact-boundary candidate");
    let rational_boundary = prepared
        .curves
        .iter()
        .find(|curve| curve.edge == rational_edge)
        .expect("rational edge should produce a prepared boundary");

    assert!(
        rational_boundary
            .curve
            .weights
            .iter()
            .any(|weight| (weight - 1.0).abs() > 1.0e-9),
        "rational boundary weights must survive compatibility preparation"
    );
}

enum SupportPcurveAdversary {
    Missing,
    Stale,
}

fn assert_support_pcurve_is_repaired(adversary: SupportPcurveAdversary) {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let (sectors, support_half_edge) = sectors_and_first_support(&brep, cap_face, &edges);
    match adversary {
        SupportPcurveAdversary::Missing => {
            brep.half_edges[support_half_edge].pcurve = None;
        }
        SupportPcurveAdversary::Stale => {
            brep.half_edges[support_half_edge].pcurve =
                Some(Pcurve::linear((100.0, 100.0), (101.0, 100.0), 0.0, 1.0));
        }
    }
    let boundaries =
        prepare_exact_boundary_curves(&brep, &sectors, ExactBoundaryCurveBudget::default())
            .unwrap();

    let ribbons = prepare_exact_boundary_ribbons(
        &brep,
        &sectors,
        &boundaries,
        ExactBoundaryRibbonBudget::default(),
    )
    .expect("support pcurve adversary should be repaired by validated projection");

    let repaired = ribbons
        .ribbons
        .iter()
        .find(|ribbon| ribbon.support_half_edge == support_half_edge)
        .expect("repaired support half-edge should produce a ribbon");
    assert!(repaired.max_support_residual <= 1.0e-6);
}

fn assert_reversed_support_pcurve_is_accepted() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let (sectors, support_half_edge) = sectors_and_first_support(&brep, cap_face, &edges);
    reverse_support_pcurve(&mut brep, support_half_edge);
    let boundaries =
        prepare_exact_boundary_curves(&brep, &sectors, ExactBoundaryCurveBudget::default())
            .unwrap();

    let ribbons = prepare_exact_boundary_ribbons(
        &brep,
        &sectors,
        &boundaries,
        ExactBoundaryRibbonBudget::default(),
    )
    .expect("orientation-reversed support pcurve should be accepted by reversal detection");

    let reversed = ribbons
        .ribbons
        .iter()
        .find(|ribbon| ribbon.support_half_edge == support_half_edge)
        .expect("reversed-storage support half-edge should produce a ribbon");
    assert!(reversed.max_support_residual <= 1.0e-6);
}

fn sectors_and_first_support(
    brep: &BRep,
    cap_face: FaceHandle,
    edges: &[EdgeHandle],
) -> (
    Vec<super::sector::ExactBoundaryCapSectorInput>,
    HalfEdgeHandle,
) {
    let plan = trihedral_plan(brep, edges);
    let sectors =
        extract_exact_boundary_cap_sectors(brep, cap_face, &plan, &trihedral_sectors(edges))
            .unwrap();
    let support_half_edge = sectors[0].fillet_seam.support_half_edge;
    (sectors, support_half_edge)
}

fn reverse_support_pcurve(brep: &mut BRep, half_edge: HalfEdgeHandle) {
    let pcurve = brep.half_edges[half_edge]
        .pcurve
        .as_ref()
        .expect("fixture support pcurve exists")
        .clone();
    let start = pcurve.point_at(pcurve.t_start);
    let end = pcurve.point_at(pcurve.t_end);
    brep.half_edges[half_edge].pcurve = Some(Pcurve::linear(
        (end.x, end.y),
        (start.x, start.y),
        pcurve.t_start,
        pcurve.t_end,
    ));
}

fn install_rational_quadratic_edge_curve(brep: &mut BRep, edge: EdgeHandle) {
    let edge_data = brep.edges[edge].clone();
    let start = brep.vertices[brep.half_edges[edge_data.half_edge].origin].position;
    let end =
        brep.vertices[brep.half_edges[brep.half_edges[edge_data.half_edge].twin].origin].position;
    let midpoint = Point3::new(
        0.5 * (start.x + end.x),
        0.5 * (start.y + end.y),
        0.5 * (start.z + end.z),
    );
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![start, midpoint, end],
            weights: vec![1.0, 0.5, 1.0],
        }),
        t_start: 0.0,
        t_end: 1.0,
    });
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
            face_edge: Some(edges[2 * sector_index + 1]),
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
            corner_kinds: Vec::new(),
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
