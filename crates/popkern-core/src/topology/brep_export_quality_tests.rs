use super::*;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::{Point3, SurfaceGeometry, SurfaceOnFace, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell,
    ShellHandle, Solid, Vertex,
};

fn options(topology_output: BrepTopologyOutputPolicy) -> BrepExportQualityOptions {
    BrepExportQualityOptions {
        topology_output,
        trim: BrepTrimRepairOptions {
            endpoint_tolerance: 1e-6,
            sampled_tolerance: 1e-6,
            samples: 5,
            validate_output: true,
        },
        require_edge_curves: false,
        allow_collapsed_trim_edge_curves: false,
        max_pcurve_deviation: Some(1e-6),
        max_edge_surface_deviation: Some(1e-6),
        validate_brep: true,
    }
}

#[test]
fn export_quality_accepts_non_pns_open_sheet_policy() {
    let brep = single_square_sheet();
    let report = analyze_brep_export_quality(&brep, options(BrepTopologyOutputPolicy::OpenSheet));

    assert_eq!(report.readiness, BrepExportReadiness::ReadyOpenSheet);
    assert!(report.is_ready());
    assert!(report.boundary_edges > 0);
    assert_eq!(report.topology_blockers, 0);
    assert_eq!(report.trim_blockers, 0);
}

#[test]
fn export_quality_rejects_open_sheet_for_closed_shell_policy() {
    let brep = single_square_sheet();
    let report = analyze_brep_export_quality(&brep, options(BrepTopologyOutputPolicy::ClosedShell));

    assert_eq!(report.readiness, BrepExportReadiness::Blocked);
    assert!(report.boundary_edges > 0);
    assert_eq!(report.topology_blockers, 1);
}

#[test]
fn export_quality_reports_missing_edge_curves_when_required() {
    let brep = single_square_sheet();
    let mut opts = options(BrepTopologyOutputPolicy::OpenSheet);
    opts.require_edge_curves = true;

    let report = analyze_brep_export_quality(&brep, opts);

    assert_eq!(report.readiness, BrepExportReadiness::Blocked);
    assert_eq!(report.missing_edge_curves, brep.edges.len());
    assert_eq!(report.edge_curve_blockers, brep.edges.len());
}

#[test]
fn export_quality_allows_collapsed_trim_missing_curve_when_enabled() {
    let brep = collapsed_trim_sheet();
    let mut opts = options(BrepTopologyOutputPolicy::OpenSheet);
    opts.require_edge_curves = true;
    opts.allow_collapsed_trim_edge_curves = true;
    opts.validate_brep = false;

    let report = analyze_brep_export_quality(&brep, opts);
    let artifact = build_brep_export_quality_artifact(&brep, &report);

    assert_eq!(report.missing_edge_curves, 1);
    assert_eq!(report.collapsed_trim_missing_edge_curves, 1);
    assert_eq!(report.edge_curve_blockers, 0);
    assert_eq!(report.authority.degenerate_topology, 1);
    assert_eq!(artifact.authority.degenerate_topology, 1);
    assert_eq!(report.readiness, BrepExportReadiness::ReadyOpenSheet);
}

#[test]
fn export_quality_collapsed_trim_authority_respects_policy_tolerance() {
    let mut brep = collapsed_trim_sheet();
    let face = brep.faces.keys().next().expect("fixture has one face");
    brep.faces[face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1e-4, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    });
    let mut opts = options(BrepTopologyOutputPolicy::OpenSheet);
    opts.require_edge_curves = true;
    opts.allow_collapsed_trim_edge_curves = true;
    opts.validate_brep = false;
    opts.trim.endpoint_tolerance = 1e-6;

    let report = analyze_brep_export_quality(&brep, opts);

    assert_eq!(report.collapsed_trim_missing_edge_curves, 0);
    assert_eq!(report.edge_curve_blockers, 1);
    assert_eq!(report.authority.degenerate_topology, 0);
    assert_eq!(report.readiness, BrepExportReadiness::Blocked);
}

#[test]
fn export_quality_artifact_lists_boundary_edges_deterministically() {
    let brep = single_square_sheet();
    let report = analyze_brep_export_quality(&brep, options(BrepTopologyOutputPolicy::OpenSheet));
    let artifact = build_brep_export_quality_artifact(&brep, &report);

    assert_eq!(artifact.summary.readiness, "ready_open_sheet");
    assert_eq!(artifact.summary.topology_output, "open_sheet");
    assert_eq!(artifact.boundary_edges.len(), brep.boundary_edges.len());
    assert!(
        artifact
            .boundary_edges
            .windows(2)
            .all(|pair| pair[0].edge <= pair[1].edge)
    );
    assert!(artifact.boundary_edges.iter().all(|edge| !edge.has_curve));
    assert_eq!(artifact.blockers.topology, 0);
    assert_eq!(artifact.blockers.trims, 0);
    assert_eq!(artifact.blockers.pcurves, 0);
    assert_eq!(artifact.blockers.edge_surfaces, 0);
}

#[test]
fn export_quality_accepts_closed_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let report = analyze_brep_export_quality(&brep, options(BrepTopologyOutputPolicy::ClosedShell));

    assert_eq!(report.readiness, BrepExportReadiness::Ready);
    assert_eq!(report.boundary_edges, 0);
    assert_eq!(report.topology_blockers, 0);
    assert_eq!(report.trim_blockers, 0);
    assert_eq!(report.pcurve_blockers, 0);
    assert_eq!(report.edge_surface_blockers, 0);
    assert!(report.max_pcurve_deviation < 1e-6);
    assert!(report.max_edge_surface_deviation < 1e-6);
}

#[test]
fn export_quality_blocks_stale_pcurves_when_tolerance_is_required() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let half_edge = brep.half_edges.keys().next().unwrap();
    brep.half_edges[half_edge].pcurve = Some(Pcurve::linear((10.0, 10.0), (11.0, 10.0), 0.0, 1.0));

    let report = analyze_brep_export_quality(&brep, options(BrepTopologyOutputPolicy::ClosedShell));
    let artifact = build_brep_export_quality_artifact(&brep, &report);

    assert_eq!(report.readiness, BrepExportReadiness::Blocked);
    assert_eq!(report.pcurve_blockers, 1);
    assert!(report.authority.stale_pcurves >= 1);
    assert!(report.max_pcurve_deviation > 1.0);
    assert_eq!(artifact.blockers.pcurves, 1);
    assert!(artifact.authority.stale_pcurves >= 1);
    assert!(artifact.summary.max_pcurve_deviation > 1.0);
}

#[test]
fn export_quality_blocks_edge_curve_off_adjacent_surface() {
    use crate::geometry::curve::{CurveEval, CurveOnEdge};

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let curve = brep.edges[edge].curve.as_ref().unwrap();
    let offset = Vec3::new(0.11, 0.07, 0.13);
    let start = curve.geometry.point_at(curve.t_start) + offset;
    let end = curve.geometry.point_at(curve.t_end) + offset;
    brep.edges[edge].curve = Some(
        CurveOnEdge::line_segment(start, end, 1e-12).expect("offset edge remains non-degenerate"),
    );

    let report = analyze_brep_export_quality(&brep, options(BrepTopologyOutputPolicy::ClosedShell));
    let artifact = build_brep_export_quality_artifact(&brep, &report);

    assert_eq!(report.readiness, BrepExportReadiness::Blocked);
    assert!(report.edge_surface_blockers >= 1);
    assert!(report.max_edge_surface_deviation > 1e-2);
    assert_eq!(
        artifact.blockers.edge_surfaces,
        report.edge_surface_blockers
    );
    assert!(artifact.summary.max_edge_surface_deviation > 1e-2);
}

fn single_square_sheet() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let square = add_square(
        &mut brep,
        shell,
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
    );
    brep.shells[shell].face = square.face;
    brep
}

fn collapsed_trim_sheet() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let vertex = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    });
    let face = brep.faces.insert(Face {
        outer_loop: LoopHandle::default(),
        inner_loops: Vec::new(),
        shell,
        surface: Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::zero(),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        }),
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    let loop_ = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face,
    });
    let half_edge = brep.half_edges.insert(HalfEdge {
        origin: vertex,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge: EdgeHandle::default(),
        loop_,
        pcurve: Some(Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0)),
    });
    let edge = brep.edges.insert(Edge {
        half_edge,
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });
    brep.half_edges[half_edge].edge = edge;
    brep.half_edges[half_edge].next = half_edge;
    brep.half_edges[half_edge].prev = half_edge;
    brep.vertices[vertex].outgoing = half_edge;
    brep.loops[loop_].half_edge = half_edge;
    brep.faces[face].outer_loop = loop_;
    brep.shells[shell].face = face;
    brep.boundary_edges.insert(edge);
    brep
}

struct SquareHandles {
    face: FaceHandle,
}

fn add_square(brep: &mut BRep, shell: ShellHandle, points: [Point3; 4]) -> SquareHandles {
    let vertices = points.map(|position| {
        brep.vertices.insert(Vertex {
            position,
            outgoing: HalfEdgeHandle::default(),
            name: None,
            tolerance: None,
        })
    });
    let face = brep.faces.insert(Face {
        outer_loop: LoopHandle::default(),
        inner_loops: Vec::new(),
        shell,
        surface: None,
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    let loop_ = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face,
    });
    brep.faces[face].outer_loop = loop_;

    let half_edges = vertices.map(|origin| {
        brep.half_edges.insert(HalfEdge {
            origin,
            twin: HalfEdgeHandle::default(),
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge: EdgeHandle::default(),
            loop_,
            pcurve: None,
        })
    });
    for half_edge in half_edges {
        let edge = brep.edges.insert(Edge {
            half_edge,
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        brep.half_edges[half_edge].edge = edge;
        brep.boundary_edges.insert(edge);
    }

    for index in 0..4 {
        let prev = (index + 3) % 4;
        let next = (index + 1) % 4;
        brep.half_edges[half_edges[index]].prev = half_edges[prev];
        brep.half_edges[half_edges[index]].next = half_edges[next];
        brep.vertices[vertices[index]].outgoing = half_edges[index];
    }
    brep.loops[loop_].half_edge = half_edges[0];
    SquareHandles { face }
}
