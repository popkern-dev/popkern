use super::boundary::{
    ExactBoundaryCurveBudget, ExactBoundaryCurveDeclineReason, prepare_exact_boundary_curves,
};
use super::sector::extract_exact_boundary_cap_sectors;
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

#[test]
fn prepares_boundaries_with_exact_degree_and_knot_operations() {
    let (brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &trihedral_sectors(&edges))
            .unwrap();

    let prepared =
        prepare_exact_boundary_curves(&brep, &sectors, ExactBoundaryCurveBudget::default())
            .unwrap();

    assert_eq!(prepared.curves.len(), 6);
    for curve in &prepared.curves {
        assert_eq!(curve.curve.degree, 9);
        assert!(curve.curve.control_points.len() <= 64);
        assert!(curve.curve.knots.len() <= 70);
    }
}

#[test]
fn budget_overflow_declines_with_stable_diagnostic() {
    let (brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &trihedral_sectors(&edges))
            .unwrap();

    let diagnostic = prepare_exact_boundary_curves(
        &brep,
        &sectors,
        ExactBoundaryCurveBudget {
            target_degree: 5,
            max_control_points: 2,
            max_knots: 4,
        },
    )
    .expect_err("budget must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryCurveDeclineReason::KnotBudgetExceeded
    );
    assert_eq!(diagnostic.code, "exact_boundary_cap.knot_budget_exceeded");
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
