use super::*;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::{Point3, SurfaceGeometry, SurfaceOnFace, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{
    Edge, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle, Solid,
    Vertex,
};

fn options() -> BrepEdgeCurveCompletionOptions {
    BrepEdgeCurveCompletionOptions {
        tolerance: 1e-6,
        samples: 5,
        allow_vertex_segments_without_surface: true,
        require_all_completed: true,
        validate_output: true,
    }
}

#[test]
fn completes_non_pns_open_sheet_edges_from_surface_trims() {
    let brep = square_sheet(true, true);
    let result = complete_missing_edge_curves_transactionally(&brep, options()).unwrap();

    assert!(
        brep.edges.iter().all(|(_, edge)| edge.curve.is_none()),
        "input BRep is untouched"
    );
    assert_eq!(result.report.inspected_missing_edges, 4);
    assert_eq!(result.report.completed_edges, 4);
    assert_eq!(result.report.declined_edges, 0);
    assert!(
        !result.report.shared_fact_ids.is_empty(),
        "edge curve completion report must expose source shared fact ids"
    );
    assert!(
        result
            .report
            .decisions
            .iter()
            .all(|decision| !decision.shared_fact_ids.is_empty())
    );
    assert!(
        result
            .brep
            .edges
            .iter()
            .all(|(_, edge)| edge.curve.is_some()),
        "repaired clone has edge curves"
    );
    assert!(result.report.decisions.iter().all(|decision| matches!(
        decision.outcome,
        BrepEdgeCurveCompletionOutcome::CompletedFromSurfaceTrims
    )));
}

#[test]
fn completes_non_surfaced_edges_from_vertices_when_policy_allows() {
    let brep = square_sheet(false, false);
    let result = complete_missing_edge_curves_transactionally(&brep, options()).unwrap();

    assert_eq!(result.report.completed_edges, 4);
    assert!(result.report.decisions.iter().all(|decision| matches!(
        decision.outcome,
        BrepEdgeCurveCompletionOutcome::CompletedFromSurfaceTrims
            | BrepEdgeCurveCompletionOutcome::CompletedFromVertices
    )));
}

#[test]
fn declines_surfaced_edges_without_trim_sources() {
    let brep = square_sheet(true, false);
    let failure = complete_missing_edge_curves_transactionally(&brep, options()).unwrap_err();

    assert_eq!(failure.report.inspected_missing_edges, 4);
    assert_eq!(failure.report.completed_edges, 0);
    assert_eq!(failure.report.declined_edges, 4);
    assert!(
        !failure.report.shared_fact_ids.is_empty(),
        "declined completion report must expose source shared fact ids"
    );
    assert!(failure.report.decisions.iter().all(|decision| matches!(
        decision.outcome,
        BrepEdgeCurveCompletionOutcome::Declined(
            BrepEdgeCurveCompletionDeclineReason::MissingSurfaceTrimSource
        )
    )));
}

#[test]
fn completes_from_vertices_when_surface_trims_are_stale_but_segment_lies_on_surfaces() {
    let mut brep = square_sheet(true, true);
    for (_, half_edge) in &mut brep.half_edges {
        half_edge.pcurve = Some(Pcurve::linear((0.25, 0.25), (1.25, 0.25), 0.0, 1.0));
    }

    let result = complete_missing_edge_curves_transactionally(&brep, options()).unwrap();

    assert_eq!(result.report.completed_edges, 4);
    assert_eq!(result.report.declined_edges, 0);
    assert!(result.report.decisions.iter().all(|decision| matches!(
        decision.outcome,
        BrepEdgeCurveCompletionOutcome::CompletedFromSurfaceTrims
            | BrepEdgeCurveCompletionOutcome::CompletedFromVertices
    )));
    assert!(
        result
            .brep
            .edges
            .iter()
            .all(|(_, edge)| edge.curve.is_some())
    );
}

fn square_sheet(with_surface: bool, with_pcurves: bool) -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let face = add_square(&mut brep, shell, with_surface, with_pcurves);
    brep.shells[shell].face = face;
    brep
}

fn add_square(
    brep: &mut BRep,
    shell: ShellHandle,
    with_surface: bool,
    with_pcurves: bool,
) -> FaceHandle {
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
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
        surface: with_surface.then_some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
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
    brep.faces[face].outer_loop = loop_;
    let half_edges = vertices.map(|origin| {
        brep.half_edges.insert(HalfEdge {
            origin,
            twin: HalfEdgeHandle::default(),
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge: Default::default(),
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
        if with_pcurves {
            let start = points[index];
            let end = points[next];
            brep.half_edges[half_edges[index]].pcurve =
                Some(Pcurve::linear((start.x, start.y), (end.x, end.y), 0.0, 1.0));
        }
    }
    brep.loops[loop_].half_edge = half_edges[0];
    face
}
