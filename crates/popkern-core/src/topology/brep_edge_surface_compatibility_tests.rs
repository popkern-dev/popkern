use super::*;
use crate::geometry::Vec3;
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::edge_triple::BoundarySide;
use crate::geometry::gregory::{GregoryPatch, RectangularGregory};
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::{Point3, SurfaceGeometry, SurfaceOnFace};
use crate::tolerance::Tolerance;
use crate::topology::brep_edge_surface_compatibility_helpers::adjacent_surfaced_pair;
use crate::topology::{
    BRep, Edge, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle,
    Solid, SurfaceBoundaryHintDeclineReason, Vertex,
};

fn options() -> BrepEdgeSurfaceCompatibilityOptions {
    BrepEdgeSurfaceCompatibilityOptions {
        tolerance: 1e-6,
        max_anchor_motion: 0.0,
        samples: 9,
        require_all_compatible: true,
        validate_output: true,
    }
}

#[test]
fn clean_unit_cube_has_no_edge_surface_findings() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let report = analyze_brep_edge_surface_compatibility(&brep, options());

    assert_eq!(report.incompatible_half_edges, 0);
    assert!(report.max_deviation <= 1e-6);
    assert_eq!(report.inspected_half_edges, 24);
}

#[test]
fn detects_edge_curve_that_does_not_lie_on_adjacent_surfaces() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = first_repairable_edge(&brep);
    corrupt_edge_curve(&mut brep, edge, Vec3::new(0.11, 0.07, 0.13));

    let report = analyze_brep_edge_surface_compatibility(&brep, options());

    assert!(report.incompatible_half_edges >= 1);
    assert!(report.findings.iter().any(|finding| finding.edge == edge));
    assert!(report.max_deviation > 1e-2);
}

#[test]
fn repair_is_transactional_and_rederives_plane_plane_edge() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = first_repairable_edge(&brep);
    let original_start = brep.edges[edge]
        .curve
        .as_ref()
        .unwrap()
        .geometry
        .point_at(brep.edges[edge].curve.as_ref().unwrap().t_start);
    corrupt_edge_curve(&mut brep, edge, Vec3::new(0.11, 0.07, 0.13));

    let result = repair_brep_edge_surface_compatibility_transactionally(&brep, options()).unwrap();
    let repaired_report = analyze_brep_edge_surface_compatibility(&result.brep, options());
    let input_report = analyze_brep_edge_surface_compatibility(&brep, options());

    assert_eq!(result.report.repaired_edges, 1);
    assert_eq!(result.report.declined_edges, 0);
    assert!(result.report.decisions[0].evidence.residual_deviation_a <= 1e-6);
    assert!(result.report.decisions[0].evidence.residual_deviation_b <= 1e-6);
    assert_eq!(repaired_report.incompatible_half_edges, 0);
    assert!(result.report.max_endpoint_gap <= 1e-6);
    assert!(
        input_report.incompatible_half_edges > 0,
        "input BRep remains corrupted"
    );

    let repaired_curve = result.brep.edges[edge].curve.as_ref().unwrap();
    assert!(
        repaired_curve
            .geometry
            .point_at(repaired_curve.t_start)
            .distance_to(&original_start)
            <= 1e-6
    );
}

#[test]
fn repair_can_reanchor_vertices_to_verified_surface_intersection_within_policy() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = first_repairable_edge(&brep);
    let half_edge = brep.edges[edge].half_edge;
    let start_vertex = brep.half_edges[half_edge].origin;
    let end_vertex = brep.half_edges[brep.half_edges[half_edge].next].origin;
    let original_start = brep.vertices[start_vertex].position;
    let original_end = brep.vertices[end_vertex].position;
    let offset = Vec3::new(0.0, 0.0, 0.02);
    brep.vertices[start_vertex].position = original_start + offset;
    brep.vertices[end_vertex].position = original_end + offset;
    corrupt_edge_curve(&mut brep, edge, offset);
    let mut options = options();
    options.max_anchor_motion = 0.08;

    let result = repair_brep_edge_surface_compatibility_transactionally(&brep, options).unwrap();

    assert_eq!(result.report.repaired_edges, 1);
    assert_eq!(result.report.after.incompatible_half_edges, 0);
    assert!(result.report.decisions[0].evidence.attempted_anchor_repair);
    assert!(result.report.decisions[0].evidence.anchor_motion <= options.max_anchor_motion);
    assert!(
        result.brep.max_imprecision >= result.report.decisions[0].evidence.anchor_motion,
        "BRep precision envelope must cover accepted vertex motion"
    );
    assert!(
        result.brep.vertices[start_vertex]
            .position
            .distance_to(&original_start)
            <= 1e-6
    );
    assert!(
        result.brep.vertices[end_vertex]
            .position
            .distance_to(&original_end)
            <= 1e-6
    );
    assert_eq!(
        brep.vertices[start_vertex].position,
        original_start + offset,
        "input BRep is untouched"
    );
}

#[test]
fn single_surfaced_plane_trim_repairs_without_adjacent_surface() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = first_repairable_edge(&brep);
    corrupt_edge_curve(&mut brep, edge, Vec3::new(0.11, 0.07, 0.13));
    let half_edge = brep.edges[edge].half_edge;
    let face = brep.loops[brep.half_edges[half_edge].loop_].face;
    brep.faces[face].surface = None;

    let result = repair_brep_edge_surface_compatibility_transactionally(&brep, options()).unwrap();

    assert_eq!(result.report.repaired_edges, 1);
    assert_eq!(result.report.declined_edges, 0);
    assert!(matches!(
        result.report.decisions[0].outcome,
        BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary
    ));
    assert_eq!(result.report.after.incompatible_half_edges, 0);
    assert!(result.report.decisions[0].evidence.surface_b.is_none());
}

#[test]
fn repairs_gregory_boundary_edge_from_surface_boundary_hint() {
    let brep = gregory_plane_boundary_fixture();
    let result =
        repair_brep_edge_surface_compatibility_transactionally(&brep, no_validation_options())
            .unwrap();

    assert_eq!(result.report.repaired_edges, 1);
    assert!(matches!(
        result.report.decisions[0].outcome,
        BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary
    ));
    assert!(
        !result.report.shared_fact_ids.is_empty(),
        "repair report must aggregate shared fact ids"
    );
    assert!(
        !result.report.decisions[0].shared_fact_ids.is_empty(),
        "edge repair decision must expose shared fact ids"
    );
    assert_eq!(result.report.vertex_stars.inspected_vertices, 2);
    assert_eq!(result.report.vertex_stars.coherent_vertices, 2);
    assert_eq!(result.report.after.incompatible_half_edges, 0);
}

#[test]
fn repairs_gregory_boundary_edge_against_nurbs_neighbor() {
    let mut brep = gregory_plane_boundary_fixture();
    let face_b =
        brep.faces
            .keys()
            .find(|&face| {
                brep.faces[face].surface.as_ref().is_some_and(|surface| {
                    matches!(surface.geometry, SurfaceGeometry::Plane { .. })
                })
            })
            .expect("fixture has analytic neighbor face");
    brep.faces[face_b].surface = Some(SurfaceOnFace {
        geometry: nurbs_patch(-1.0),
        orientation: true,
    });

    let result =
        repair_brep_edge_surface_compatibility_transactionally(&brep, no_validation_options())
            .unwrap();

    assert_eq!(result.report.repaired_edges, 1);
    assert_eq!(result.report.after.incompatible_half_edges, 0);
    assert_eq!(result.report.vertex_stars.coherent_vertices, 2);
    assert!(matches!(
        result.report.decisions[0].outcome,
        BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary
    ));
}

#[test]
fn invalid_authored_boundary_hint_is_downgraded_before_repair() {
    let mut brep = gregory_plane_boundary_fixture();
    let edge = brep.edges.keys().next().expect("fixture has one edge");
    let he = brep.edges[edge].half_edge;
    brep.half_edges[he].pcurve = Some(Pcurve::linear((0.25, 0.0), (0.25, 1.0), 0.0, 1.0));
    assert!(brep.set_surface_boundary_hint(he, BoundarySide::UMin));

    let result =
        repair_brep_edge_surface_compatibility_transactionally(&brep, no_validation_options())
            .unwrap();

    assert_eq!(result.report.repaired_edges, 1);
    assert!(result.brep.surface_boundary_hint(he).is_none());
    let decline = result
        .brep
        .surface_boundary_decline(he)
        .expect("invalid authored hint must be preserved as a structured decline");
    assert_eq!(
        decline.reason,
        SurfaceBoundaryHintDeclineReason::AuthoredHintBoundaryMismatch
    );
    assert_eq!(decline.best_boundary, Some(BoundarySide::UMin));
    assert!(decline.best_boundary_deviation > 0.0);
    assert!(matches!(
        result.report.decisions[0].outcome,
        BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary
    ));
}

#[test]
fn repairs_nurbs_boundary_edge_against_nurbs_neighbor() {
    let brep = nurbs_nurbs_boundary_fixture();
    let result =
        repair_brep_edge_surface_compatibility_transactionally(&brep, no_validation_options())
            .unwrap();

    assert_eq!(result.report.repaired_edges, 1);
    assert_eq!(result.report.after.incompatible_half_edges, 0);
    assert!(matches!(
        result.report.decisions[0].outcome,
        BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary
    ));
}

fn first_repairable_edge(brep: &BRep) -> EdgeHandle {
    brep.edges
        .keys()
        .find(|&edge| adjacent_surfaced_pair(brep, edge).is_some())
        .expect("unit cube has repairable edge")
}

fn no_validation_options() -> BrepEdgeSurfaceCompatibilityOptions {
    BrepEdgeSurfaceCompatibilityOptions {
        validate_output: false,
        ..options()
    }
}

fn nurbs_nurbs_boundary_fixture() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let start = Point3::new(0.0, 0.0, 0.0);
    let end = Point3::new(0.0, 1.0, 0.0);
    let va = brep.vertices.insert(vertex(start));
    let vb = brep.vertices.insert(vertex(end));
    let face_a = brep.faces.insert(face(
        shell,
        SurfaceOnFace {
            geometry: nurbs_patch(1.0),
            orientation: true,
        },
    ));
    let face_b = brep.faces.insert(face(
        shell,
        SurfaceOnFace {
            geometry: nurbs_patch(-1.0),
            orientation: true,
        },
    ));
    let loop_a = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: face_a,
    });
    let loop_b = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: face_b,
    });
    let he_a = brep.half_edges.insert(half_edge(va, loop_a));
    let he_b = brep.half_edges.insert(half_edge(vb, loop_b));
    let edge = brep.edges.insert(Edge {
        half_edge: he_a,
        curve: CurveOnEdge::line_segment(
            start + Vec3::new(0.0, 0.0, 0.1),
            end + Vec3::new(0.0, 0.0, 0.1),
            1e-12,
        ),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });
    brep.half_edges[he_a].edge = edge;
    brep.half_edges[he_b].edge = edge;
    brep.half_edges[he_a].twin = he_b;
    brep.half_edges[he_b].twin = he_a;
    brep.half_edges[he_a].next = he_b;
    brep.half_edges[he_b].next = he_a;
    brep.half_edges[he_a].prev = he_b;
    brep.half_edges[he_b].prev = he_a;
    brep.half_edges[he_a].pcurve = Some(Pcurve::linear((0.0, 0.0), (0.0, 1.0), 0.0, 1.0));
    brep.half_edges[he_b].pcurve = Some(Pcurve::linear((0.0, 0.0), (0.0, 1.0), 0.0, 1.0));
    assert!(brep.set_surface_boundary_hint(he_a, BoundarySide::UMin));
    brep.loops[loop_a].half_edge = he_a;
    brep.loops[loop_b].half_edge = he_b;
    brep.faces[face_a].outer_loop = loop_a;
    brep.faces[face_b].outer_loop = loop_b;
    brep.vertices[va].outgoing = he_a;
    brep.vertices[vb].outgoing = he_b;
    brep.shells[shell].face = face_a;
    brep
}

fn nurbs_patch(x_extent: f64) -> SurfaceGeometry {
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(x_extent, 0.0, 0.0),
            Point3::new(x_extent, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    })
}

fn gregory_plane_boundary_fixture() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let start = Point3::new(0.0, 0.0, 0.0);
    let end = Point3::new(0.0, 1.0, 0.0);
    let va = brep.vertices.insert(vertex(start));
    let vb = brep.vertices.insert(vertex(end));
    let face_a = brep.faces.insert(face(shell, gregory_surface(start, end)));
    let face_b = brep.faces.insert(face(shell, plane_surface()));
    let loop_a = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: face_a,
    });
    let loop_b = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: face_b,
    });
    let he_a = brep.half_edges.insert(half_edge(va, loop_a));
    let he_b = brep.half_edges.insert(half_edge(vb, loop_b));
    let edge = brep.edges.insert(Edge {
        half_edge: he_a,
        curve: CurveOnEdge::line_segment(
            start + Vec3::new(0.0, 0.0, 0.1),
            end + Vec3::new(0.0, 0.0, 0.1),
            1e-12,
        ),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });
    brep.half_edges[he_a].edge = edge;
    brep.half_edges[he_b].edge = edge;
    brep.half_edges[he_a].twin = he_b;
    brep.half_edges[he_b].twin = he_a;
    brep.half_edges[he_a].next = he_b;
    brep.half_edges[he_b].next = he_a;
    brep.half_edges[he_a].prev = he_b;
    brep.half_edges[he_b].prev = he_a;
    brep.half_edges[he_a].pcurve = Some(Pcurve::linear((0.0, 0.0), (0.0, 1.0), 0.0, 1.0));
    brep.half_edges[he_b].pcurve = Some(Pcurve::linear((0.0, 0.0), (0.0, 1.0), 0.0, 1.0));
    brep.loops[loop_a].half_edge = he_a;
    brep.loops[loop_b].half_edge = he_b;
    brep.faces[face_a].outer_loop = loop_a;
    brep.faces[face_b].outer_loop = loop_b;
    brep.vertices[va].outgoing = he_a;
    brep.vertices[vb].outgoing = he_b;
    brep.shells[shell].face = face_a;
    brep
}

fn vertex(position: Point3) -> Vertex {
    Vertex {
        position,
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    }
}

fn half_edge(origin: crate::topology::VertexHandle, loop_: LoopHandle) -> HalfEdge {
    HalfEdge {
        origin,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge: Default::default(),
        loop_,
        pcurve: None,
    }
}

fn face(shell: ShellHandle, surface: SurfaceOnFace) -> Face {
    Face {
        outer_loop: LoopHandle::default(),
        inner_loops: Vec::new(),
        shell,
        surface: Some(surface),
        name: None,
        tolerance: None,
        role: Default::default(),
    }
}

fn gregory_surface(start: Point3, end: Point3) -> SurfaceOnFace {
    let lerp = |s: f64| start + (end - start) * s;
    let offset = Vec3::new(1.0, 0.0, 0.0);
    let u0 = [lerp(0.0), lerp(1.0 / 3.0), lerp(2.0 / 3.0), lerp(1.0)];
    let row = |scale: f64| {
        [
            u0[0] + offset * scale,
            u0[1] + offset * scale,
            u0[2] + offset * scale,
            u0[3] + offset * scale,
        ]
    };
    SurfaceOnFace {
        geometry: SurfaceGeometry::Gregory(GregoryPatch::Rectangular(Box::new(
            RectangularGregory {
                boundary_cps: [u0, row(1.0 / 3.0), row(2.0 / 3.0), row(1.0)],
                inner_cps: [(row(1.0 / 3.0)[1], row(1.0 / 3.0)[1]); 4],
            },
        ))),
        orientation: true,
    }
}

fn plane_surface() -> SurfaceOnFace {
    SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    }
}

fn corrupt_edge_curve(brep: &mut BRep, edge: EdgeHandle, offset: Vec3) {
    let curve = brep.edges[edge].curve.as_ref().unwrap();
    let start = curve.geometry.point_at(curve.t_start) + offset;
    let end = curve.geometry.point_at(curve.t_end) + offset;
    brep.edges[edge].curve = Some(
        CurveOnEdge::line_segment(start, end, 1e-12)
            .expect("offset keeps cube edge non-degenerate"),
    );
}

#[path = "brep_edge_surface_compatibility_anchor_tests.rs"]
mod anchor_tests;
#[path = "brep_edge_surface_compatibility_radial_tests.rs"]
mod radial_tests;
#[path = "brep_edge_surface_compatibility_trim_tests.rs"]
mod trim_tests;
