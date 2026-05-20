use super::sector::{
    ExactBoundaryOuterSide, ExactBoundarySectorDeclineReason, extract_exact_boundary_cap_sectors,
};
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
fn extracts_trihedral_sectors_in_planner_order() {
    let (brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors = trihedral_sectors(&edges);

    let extracted = extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors).unwrap();

    assert_eq!(extracted.len(), 3);
    for (index, sector) in extracted.iter().enumerate() {
        assert_eq!(sector.sector_index, index);
        assert_eq!(
            sector.sector_name,
            format!("direct_apex.cap_sector.{index}")
        );
        assert_eq!(sector.fillet_index, index);
        assert_eq!(sector.fillet_end_edge, edges[2 * index]);
        assert_eq!(sector.face_edge, edges[2 * index + 1]);
        assert_eq!(
            sector.outer_seams[0].side,
            ExactBoundaryOuterSide::FilletEnd
        );
        assert_eq!(sector.outer_seams[1].side, ExactBoundaryOuterSide::FaceEdge);
        assert_ne!(sector.outer_seams[0].edge, sector.outer_seams[1].edge);
        assert!(brep.faces.contains_key(sector.fillet_support_face));
        assert!(brep.faces.contains_key(sector.bridge_support_face));
    }
}

#[test]
fn stale_sector_edge_declines_with_stable_diagnostic() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors = trihedral_sectors(&edges);
    brep.edges.remove(edges[2]);

    let diagnostic = extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors)
        .expect_err("stale edge declines");

    assert_eq!(
        diagnostic.reason,
        ExactBoundarySectorDeclineReason::StaleEdge
    );
    assert_eq!(diagnostic.code, "exact_boundary_cap.sector.stale_edge");
    assert_eq!(diagnostic.sector_index, Some(1));
    assert_eq!(
        diagnostic.sector_name.as_deref(),
        Some("direct_apex.cap_sector.1")
    );
    assert_eq!(diagnostic.edge, Some(edges[2]));
}

fn hex_cap_fixture() -> (BRep, FaceHandle, Vec<EdgeHandle>) {
    let brep = BRep::n_gonal_pyramid(6, 1.0, 0.8, Tolerance::default()).unwrap();
    let cap_face = brep
        .faces
        .keys()
        .find(|&face| loop_edges(&brep, face).len() == 6)
        .expect("hexagonal base face exists");
    let edges = loop_edges(&brep, cap_face);
    assert_eq!(edges.len(), 6);
    for edge in &edges {
        assert!(brep.edges[*edge].curve.is_some());
    }
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
