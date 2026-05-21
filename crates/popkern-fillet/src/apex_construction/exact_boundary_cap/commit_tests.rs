use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{NurbsSurface, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, ShellHandle, TopologyRole,
    VertexHandle,
};

use super::commit::{ExactBoundaryCommitDeclineReason, commit_exact_boundary_sector_surfaces};
use super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide};
use super::solve::{ExactBoundarySectorSurface, ExactBoundarySolvedSurfaces};

#[test]
fn commit_declines_when_sector_topology_is_not_split() {
    let mut brep = BRep::new(Tolerance::default());
    let aggregate_cap = insert_face(&mut brep, ShellHandle::default());
    let cap_half_edges = insert_half_edge_pair(&mut brep);
    let support_half_edges = insert_half_edge_pair(&mut brep);
    let sector = sector_input(
        0,
        "sector_0",
        aggregate_cap,
        cap_half_edges,
        support_half_edges,
    );

    let diagnostic = commit_exact_boundary_sector_surfaces(
        &mut brep,
        aggregate_cap,
        &[sector],
        &solved_surfaces(),
        &prepared_pcurves(),
    )
    .expect_err("aggregate cap topology must not be accepted as sector topology");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryCommitDeclineReason::SectorTopologyMissing
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.commit.sector_topology_missing"
    );
    assert!(brep.faces[aggregate_cap].surface.is_none());
}

#[test]
fn commit_writes_sector_surfaces_and_pcurves_for_split_topology() {
    let mut brep = BRep::new(Tolerance::default());
    let aggregate_cap = insert_face(&mut brep, ShellHandle::default());
    let sector_face = insert_face(&mut brep, ShellHandle::default());
    let cap_half_edges = insert_half_edge_pair(&mut brep);
    let support_half_edges = insert_half_edge_pair(&mut brep);
    let sector = sector_input(
        0,
        "sector_0",
        sector_face,
        cap_half_edges,
        support_half_edges,
    );

    let report = commit_exact_boundary_sector_surfaces(
        &mut brep,
        aggregate_cap,
        &[sector],
        &solved_surfaces(),
        &prepared_pcurves(),
    )
    .expect("pre-split sector topology can receive exact-boundary candidates");

    assert_eq!(report.committed_sector_count, 1);
    assert_eq!(report.committed_pcurve_pair_count, 2);
    assert!(brep.faces[sector_face].surface.is_some());
    assert!(brep.half_edges[cap_half_edges[0]].pcurve.is_some());
    assert!(brep.half_edges[cap_half_edges[1]].pcurve.is_some());
    assert!(brep.half_edges[support_half_edges[0]].pcurve.is_some());
    assert!(brep.half_edges[support_half_edges[1]].pcurve.is_some());
}

#[test]
fn commit_declines_when_pcurve_pairs_do_not_cover_each_sector_side() {
    let mut brep = BRep::new(Tolerance::default());
    let aggregate_cap = insert_face(&mut brep, ShellHandle::default());
    let sector_face = insert_face(&mut brep, ShellHandle::default());
    let cap_half_edges = insert_half_edge_pair(&mut brep);
    let support_half_edges = insert_half_edge_pair(&mut brep);
    let sector = sector_input(
        0,
        "sector_0",
        sector_face,
        cap_half_edges,
        support_half_edges,
    );
    let pcurves = PreparedExactBoundaryPcurves {
        pairs: vec![
            pcurve_pair(ExactBoundaryOuterSide::FilletEnd),
            pcurve_pair(ExactBoundaryOuterSide::FilletEnd),
        ],
    };

    let diagnostic = commit_exact_boundary_sector_surfaces(
        &mut brep,
        aggregate_cap,
        &[sector],
        &solved_surfaces(),
        &pcurves,
    )
    .expect_err("duplicated side pairs must not satisfy per-sector pcurve ownership");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryCommitDeclineReason::PcurveInventoryMismatch
    );
    assert!(brep.faces[sector_face].surface.is_none());
}

#[test]
fn commit_diagnostic_promotes_to_apex_diagnostic() {
    let diagnostic = super::commit::ExactBoundaryCommitDiagnostic {
        reason: ExactBoundaryCommitDeclineReason::SectorTopologyMissing,
        code: "exact_boundary_cap.commit.sector_topology_missing",
        sector_name: Some("sector_0".to_string()),
        face: Some(FaceHandle::default()),
    };

    let crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(apex) =
        diagnostic.into_apex_diagnostic()
    else {
        panic!("commit diagnostic must promote through the exact-boundary diagnostic channel");
    };
    assert_eq!(apex.stage, "exact_boundary_cap.commit");
    assert_eq!(
        apex.code,
        "exact_boundary_cap.commit.sector_topology_missing"
    );
    assert_eq!(apex.sector_name, "sector_0");
}

fn insert_face(brep: &mut BRep, shell: ShellHandle) -> FaceHandle {
    let face = brep.faces.insert(Face {
        outer_loop: LoopHandle::default(),
        inner_loops: Vec::new(),
        shell,
        surface: None,
        name: None,
        tolerance: None,
        role: TopologyRole::Production,
    });
    let loop_ = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face,
    });
    brep.faces[face].outer_loop = loop_;
    face
}

fn insert_half_edge_pair(brep: &mut BRep) -> [HalfEdgeHandle; 2] {
    let loop_ = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: FaceHandle::default(),
    });
    let first = brep.half_edges.insert(half_edge(loop_));
    let second = brep.half_edges.insert(half_edge(loop_));
    brep.half_edges[first].twin = second;
    brep.half_edges[first].next = first;
    brep.half_edges[first].prev = first;
    brep.half_edges[second].twin = first;
    brep.half_edges[second].next = second;
    brep.half_edges[second].prev = second;
    brep.loops[loop_].half_edge = first;
    [first, second]
}

fn half_edge(loop_: LoopHandle) -> HalfEdge {
    HalfEdge {
        origin: VertexHandle::default(),
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge: Default::default(),
        loop_,
        pcurve: None,
    }
}

fn sector_input(
    sector_index: usize,
    sector_name: &str,
    cap_face: FaceHandle,
    cap_half_edges: [HalfEdgeHandle; 2],
    support_half_edges: [HalfEdgeHandle; 2],
) -> ExactBoundaryCapSectorInput {
    ExactBoundaryCapSectorInput {
        sector_index,
        sector_name: sector_name.to_string(),
        fillet_index: sector_index,
        cap_face,
        fillet_end_edge: Default::default(),
        face_edge: Default::default(),
        fillet_support_face: FaceHandle::default(),
        bridge_support_face: FaceHandle::default(),
        cap_half_edges,
        support_half_edges,
        outer_seams: [
            outer_seam(
                ExactBoundaryOuterSide::FilletEnd,
                cap_half_edges[0],
                support_half_edges[0],
            ),
            outer_seam(
                ExactBoundaryOuterSide::FaceEdge,
                cap_half_edges[1],
                support_half_edges[1],
            ),
        ],
    }
}

fn outer_seam(
    side: ExactBoundaryOuterSide,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
) -> ExactBoundaryOuterSeam {
    ExactBoundaryOuterSeam {
        side,
        edge: Default::default(),
        cap_half_edge,
        support_half_edge,
        support_face: FaceHandle::default(),
        cap_pcurve_present: false,
        support_pcurve_present: false,
    }
}

fn solved_surfaces() -> ExactBoundarySolvedSurfaces {
    ExactBoundarySolvedSurfaces {
        surfaces: vec![ExactBoundarySectorSurface {
            sector_index: 0,
            sector_name: "sector_0".to_string(),
            surface: unit_square_surface(),
            max_corner_gap: 0.0,
        }],
    }
}

fn unit_square_surface() -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        weights: vec![1.0; 4],
        num_u: 2,
        num_v: 2,
    }
}

fn prepared_pcurves() -> PreparedExactBoundaryPcurves {
    PreparedExactBoundaryPcurves {
        pairs: vec![
            pcurve_pair(ExactBoundaryOuterSide::FilletEnd),
            pcurve_pair(ExactBoundaryOuterSide::FaceEdge),
        ],
    }
}

fn pcurve_pair(side: ExactBoundaryOuterSide) -> PreparedExactBoundaryPcurvePair {
    PreparedExactBoundaryPcurvePair {
        sector_index: 0,
        sector_name: "sector_0".to_string(),
        side,
        edge: Default::default(),
        support_half_edge: HalfEdgeHandle::default(),
        support_face: FaceHandle::default(),
        cap_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
        support_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
        max_cap_residual: 0.0,
        max_support_residual: 0.0,
        max_pair_residual: 0.0,
    }
}
