use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
    ShellHandle, TopologyRole, Vertex, VertexHandle,
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
        [EdgeHandle::default(), EdgeHandle::default()],
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
    let support_face = insert_support_face(&mut brep, ShellHandle::default());
    let seams = insert_shared_trim_seams(&mut brep, sector_face, support_face);
    let sector = sector_input(
        0,
        "sector_0",
        sector_face,
        seams.edges,
        seams.cap_half_edges,
        seams.support_half_edges,
    );

    let report = commit_exact_boundary_sector_surfaces(
        &mut brep,
        aggregate_cap,
        &[sector],
        &solved_surfaces(),
        &prepared_pcurves_for(&seams),
    )
    .expect("pre-split sector topology can receive exact-boundary candidates");

    assert_eq!(report.committed_sector_count, 1);
    assert_eq!(report.committed_pcurve_pair_count, 2);
    assert!(brep.faces[sector_face].surface.is_some());
    assert!(brep.half_edges[seams.cap_half_edges[0]].pcurve.is_some());
    assert!(brep.half_edges[seams.cap_half_edges[1]].pcurve.is_some());
    assert!(
        brep.half_edges[seams.support_half_edges[0]]
            .pcurve
            .is_some()
    );
    assert!(
        brep.half_edges[seams.support_half_edges[1]]
            .pcurve
            .is_some()
    );
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
        [EdgeHandle::default(), EdgeHandle::default()],
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
fn commit_materializes_support_subdivision_for_unsurfaced_exact_seam() {
    let mut brep = BRep::new(Tolerance::default());
    let aggregate_cap = insert_face(&mut brep, ShellHandle::default());
    let sector_face = insert_face(&mut brep, ShellHandle::default());
    let support_face = insert_face(&mut brep, ShellHandle::default());
    let seams = insert_shared_trim_seams(&mut brep, sector_face, support_face);
    let sector = sector_input(
        0,
        "sector_0",
        sector_face,
        seams.edges,
        seams.cap_half_edges,
        seams.support_half_edges,
    );

    let report = commit_exact_boundary_sector_surfaces(
        &mut brep,
        aggregate_cap,
        &[sector],
        &solved_surfaces(),
        &prepared_pcurves_for(&seams),
    )
    .expect("exact-boundary commit should materialize a support-owned seam patch");

    assert_eq!(report.committed_sector_count, 1);
    assert_eq!(report.committed_pcurve_pair_count, 2);
    assert!(brep.faces[sector_face].surface.is_some());
    assert!(brep.half_edges[seams.cap_half_edges[0]].pcurve.is_some());
    let support_half_edge = brep.half_edges[seams.cap_half_edges[0]].twin;
    assert_ne!(support_half_edge, seams.support_half_edges[0]);
    assert_eq!(
        brep.half_edges[support_half_edge].twin,
        seams.cap_half_edges[0]
    );
    let support_patch_loop = brep.half_edges[support_half_edge].loop_;
    let support_patch_face = brep.loops[support_patch_loop].face;
    assert_ne!(support_patch_face, support_face);
    assert!(brep.faces[support_patch_face].surface.is_some());
    let previous_neighbor_twin = brep.half_edges[seams.support_half_edges[0]].twin;
    assert_ne!(previous_neighbor_twin, seams.cap_half_edges[0]);
    assert_eq!(
        brep.half_edges[previous_neighbor_twin].twin,
        seams.support_half_edges[0]
    );
    assert!(
        !brep
            .boundary_edges
            .contains(&brep.half_edges[seams.support_half_edges[0]].edge)
    );
}

#[test]
fn commit_assigns_prepared_support_surface_when_whole_loop_fits() {
    let mut brep = BRep::new(Tolerance::default());
    let aggregate_cap = insert_face(&mut brep, ShellHandle::default());
    let sector_face = insert_face(&mut brep, ShellHandle::default());
    let support_face = insert_face(&mut brep, ShellHandle::default());
    let seams = insert_shared_trim_seams(&mut brep, sector_face, support_face);
    set_line_edge_curves(&mut brep, &seams);
    let sector = sector_input(
        0,
        "sector_0",
        sector_face,
        seams.edges,
        seams.cap_half_edges,
        seams.support_half_edges,
    );

    let report = commit_exact_boundary_sector_surfaces(
        &mut brep,
        aggregate_cap,
        &[sector],
        &solved_surfaces(),
        &prepared_pcurves_for(&seams),
    )
    .expect("loop-valid prepared support surface can be committed to unsurfaced support topology");

    assert_eq!(report.committed_sector_count, 1);
    assert!(brep.faces[support_face].surface.is_some());
    assert!(
        brep.half_edges[seams.support_half_edges[0]]
            .pcurve
            .is_some()
    );
    assert!(
        brep.half_edges[seams.support_half_edges[1]]
            .pcurve
            .is_some()
    );
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

fn insert_support_face(brep: &mut BRep, shell: ShellHandle) -> FaceHandle {
    let face = insert_face(brep, shell);
    brep.faces[face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(unit_square_surface()),
        orientation: true,
    });
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

struct SharedTrimSeams {
    edges: [EdgeHandle; 2],
    cap_half_edges: [HalfEdgeHandle; 2],
    support_half_edges: [HalfEdgeHandle; 2],
}

fn insert_shared_trim_seams(
    brep: &mut BRep,
    cap_face: FaceHandle,
    support_face: FaceHandle,
) -> SharedTrimSeams {
    let first = insert_shared_trim_seam(
        brep,
        cap_face,
        support_face,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    let second = insert_shared_trim_seam(
        brep,
        cap_face,
        support_face,
        Point3::new(0.0, 0.25, 0.0),
        Point3::new(1.0, 0.25, 0.0),
    );
    SharedTrimSeams {
        edges: [first.0, second.0],
        cap_half_edges: [first.1, second.1],
        support_half_edges: [first.2, second.2],
    }
}

fn insert_shared_trim_seam(
    brep: &mut BRep,
    cap_face: FaceHandle,
    support_face: FaceHandle,
    start: Point3,
    end: Point3,
) -> (EdgeHandle, HalfEdgeHandle, HalfEdgeHandle) {
    let start_vertex = insert_vertex(brep, start);
    let end_vertex = insert_vertex(brep, end);
    let edge = brep.edges.insert(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: TopologyRole::Production,
    });
    let cap_loop = brep.faces[cap_face].outer_loop;
    let support_loop = brep.faces[support_face].outer_loop;
    let cap_half_edge = brep.half_edges.insert(half_edge(cap_loop));
    let support_half_edge = brep.half_edges.insert(half_edge(support_loop));
    brep.half_edges[cap_half_edge].origin = start_vertex;
    brep.half_edges[cap_half_edge].twin = support_half_edge;
    brep.half_edges[cap_half_edge].edge = edge;
    brep.half_edges[cap_half_edge].next = cap_half_edge;
    brep.half_edges[cap_half_edge].prev = cap_half_edge;
    brep.half_edges[support_half_edge].origin = end_vertex;
    brep.half_edges[support_half_edge].twin = cap_half_edge;
    brep.half_edges[support_half_edge].edge = edge;
    brep.half_edges[support_half_edge].next = support_half_edge;
    brep.half_edges[support_half_edge].prev = support_half_edge;
    brep.edges[edge].half_edge = cap_half_edge;
    brep.loops[cap_loop].half_edge = cap_half_edge;
    brep.loops[support_loop].half_edge = support_half_edge;
    brep.vertices[start_vertex].outgoing = cap_half_edge;
    brep.vertices[end_vertex].outgoing = support_half_edge;
    (edge, cap_half_edge, support_half_edge)
}

fn set_line_edge_curves(brep: &mut BRep, seams: &SharedTrimSeams) {
    for edge in seams.edges {
        let half_edge = brep.edges[edge].half_edge;
        let start = brep.vertices[brep.half_edges[half_edge].origin].position;
        let end_vertex = brep.half_edges[brep.half_edges[half_edge].twin].origin;
        let end = brep.vertices[end_vertex].position;
        brep.edges[edge].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: start,
                direction: end - start,
            },
            t_start: 0.0,
            t_end: 1.0,
        });
    }
}

fn insert_vertex(brep: &mut BRep, position: Point3) -> VertexHandle {
    brep.vertices.insert(Vertex {
        position,
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    })
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
    edges: [EdgeHandle; 2],
    cap_half_edges: [HalfEdgeHandle; 2],
    support_half_edges: [HalfEdgeHandle; 2],
) -> ExactBoundaryCapSectorInput {
    ExactBoundaryCapSectorInput {
        sector_index,
        sector_name: sector_name.to_string(),
        fillet_index: sector_index,
        cap_face,
        fillet_end_edge: edges[0],
        face_edge: Some(edges[1]),
        fillet_support_face: FaceHandle::default(),
        bridge_support_face: Some(FaceHandle::default()),
        fillet_seam: outer_seam(
            ExactBoundaryOuterSide::FilletEnd,
            edges[0],
            cap_half_edges[0],
            support_half_edges[0],
        ),
        bridge_seam: Some(outer_seam(
            ExactBoundaryOuterSide::FaceEdge,
            edges[1],
            cap_half_edges[1],
            support_half_edges[1],
        )),
    }
}

fn outer_seam(
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
) -> ExactBoundaryOuterSeam {
    ExactBoundaryOuterSeam {
        side,
        edge,
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
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
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

fn prepared_pcurves_for(seams: &SharedTrimSeams) -> PreparedExactBoundaryPcurves {
    let mut pcurves = prepared_pcurves();
    pcurves.pairs[0].edge = seams.edges[0];
    pcurves.pairs[0].support_half_edge = seams.support_half_edges[0];
    pcurves.pairs[0].support_pcurve = Pcurve::linear((1.0, 0.0), (0.0, 0.0), 0.0, 1.0);
    pcurves.pairs[1].edge = seams.edges[1];
    pcurves.pairs[1].support_half_edge = seams.support_half_edges[1];
    pcurves.pairs[1].cap_pcurve = Pcurve::linear((0.0, 0.25), (1.0, 0.25), 0.0, 1.0);
    pcurves.pairs[1].support_pcurve = Pcurve::linear((1.0, 0.25), (0.0, 0.25), 0.0, 1.0);
    pcurves
}

fn pcurve_pair(side: ExactBoundaryOuterSide) -> PreparedExactBoundaryPcurvePair {
    PreparedExactBoundaryPcurvePair {
        sector_index: 0,
        sector_name: "sector_0".to_string(),
        side,
        edge: Default::default(),
        support_half_edge: HalfEdgeHandle::default(),
        support_face: FaceHandle::default(),
        support_surface: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        cap_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
        support_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
        max_cap_residual: 0.0,
        max_support_residual: 0.0,
        max_pair_residual: 0.0,
    }
}
