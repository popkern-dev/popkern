use super::boundary::{ExactBoundaryCurveBudget, prepare_exact_boundary_curves};
use super::corner::{ExactBoundaryInnerSpokeBudget, prepare_exact_boundary_inner_spokes_with_pole};
use super::sector::extract_exact_boundary_cap_sectors;
use super::topology::materialize_exact_boundary_sector_topology;
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
fn materializes_sector_faces_around_shared_pole() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &trihedral_sectors(&edges))
            .unwrap();
    let boundaries =
        prepare_exact_boundary_curves(&brep, &sectors, ExactBoundaryCurveBudget::default())
            .unwrap();
    let spokes = prepare_exact_boundary_inner_spokes_with_pole(
        &sectors,
        &boundaries,
        ExactBoundaryInnerSpokeBudget::default(),
        None,
    )
    .unwrap();
    let before = (
        brep.vertices.len(),
        brep.edges.len(),
        brep.half_edges.len(),
        brep.faces.len(),
    );

    let topology =
        materialize_exact_boundary_sector_topology(&mut brep, cap_face, &sectors, &spokes)
            .expect("hex cap can be split into exact-boundary sectors");

    assert_eq!(topology.sectors.len(), 3);
    assert_eq!(topology.spoke_edges.len(), 3);
    assert!(brep.vertices.contains_key(topology.pole_vertex));
    assert_eq!(brep.vertices.len(), before.0 + 1);
    assert_eq!(brep.edges.len(), before.1 + 3);
    assert_eq!(brep.half_edges.len(), before.2 + 6);
    assert_eq!(brep.faces.len(), before.3 + 2);
    assert_eq!(topology.sectors[0].cap_face, cap_face);
    assert!(
        topology.sectors[1..]
            .iter()
            .all(|sector| sector.cap_face != cap_face)
    );
    brep.validate()
        .expect("materialized sector topology validates");
}

#[test]
fn topology_decline_rolls_back_on_stale_cap_face() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors =
        extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &trihedral_sectors(&edges))
            .unwrap();
    let boundaries =
        prepare_exact_boundary_curves(&brep, &sectors, ExactBoundaryCurveBudget::default())
            .unwrap();
    let spokes = prepare_exact_boundary_inner_spokes_with_pole(
        &sectors,
        &boundaries,
        ExactBoundaryInnerSpokeBudget::default(),
        None,
    )
    .unwrap();
    let before = (
        brep.vertices.len(),
        brep.edges.len(),
        brep.half_edges.len(),
        brep.faces.len(),
    );
    brep.faces.remove(cap_face);

    let diagnostic =
        materialize_exact_boundary_sector_topology(&mut brep, cap_face, &sectors, &spokes)
            .expect_err("stale cap face must decline");

    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.topology.cap_face_missing"
    );
    assert_eq!(
        (
            brep.vertices.len(),
            brep.edges.len(),
            brep.half_edges.len(),
            brep.faces.len()
        ),
        (before.0, before.1, before.2, before.3 - 1)
    );
}

fn hex_cap_fixture() -> (BRep, FaceHandle, Vec<EdgeHandle>) {
    let brep = BRep::n_gonal_pyramid(6, 1.0, 0.8, Tolerance::default()).unwrap();
    let cap_face = brep
        .faces
        .keys()
        .find(|&face| loop_edges(&brep, face).len() == 6)
        .expect("hexagonal base face exists");
    let edges = loop_edges(&brep, cap_face);
    (brep, cap_face, edges)
}

fn loop_edges(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let loop_handle = brep.faces[face].outer_loop;
    let start = brep.loops[loop_handle].half_edge;
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
