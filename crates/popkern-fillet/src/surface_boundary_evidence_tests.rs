use super::*;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregory};
use popkern_core::geometry::{Point3, SurfaceOnFace, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    Edge, Face, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle, Solid,
    SurfaceBoundaryAuditOutcome, Vertex,
};

#[test]
fn audit_authors_exact_triple_for_canonical_gregory_boundary() {
    let mut brep = gregory_plane_boundary_fixture();
    let edge = brep.edges.keys().next().unwrap();
    let half_edge = brep.edges[edge].half_edge;
    brep.half_edges[half_edge].pcurve = Some(Pcurve::linear((0.25, 0.0), (0.25, 1.0), 0.0, 1.0));

    audit_gregory_analytic_boundary_evidence(&mut brep);

    assert_eq!(
        brep.surface_boundary_hint(half_edge),
        Some(BoundarySide::UMin)
    );
    assert!(brep.surface_boundary_decline(half_edge).is_none());
    let pcurve = brep.half_edges[half_edge].pcurve.as_ref().unwrap();
    let (t0, t1) = pcurve.parameter_range();
    assert_eq!(pcurve.point_at(t0).x, 0.0);
    assert_eq!(pcurve.point_at(t1).x, 0.0);
    assert!(brep.surface_boundary_audit_events.iter().any(|event| {
        event.half_edge == half_edge
            && event.outcome == SurfaceBoundaryAuditOutcome::SharedTrimAccepted
    }));
}

fn gregory_plane_boundary_fixture() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: Default::default(),
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

    let he_a_return = brep.half_edges.insert(half_edge(vb, loop_a));
    let he_b_return = brep.half_edges.insert(half_edge(va, loop_b));
    let return_edge = brep.edges.insert(Edge {
        half_edge: he_a_return,
        curve: CurveOnEdge::line_segment(end, start, 1e-12),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });
    brep.half_edges[he_a_return].edge = return_edge;
    brep.half_edges[he_b_return].edge = return_edge;
    brep.half_edges[he_a_return].twin = he_b_return;
    brep.half_edges[he_b_return].twin = he_a_return;

    brep.half_edges[he_a].next = he_a_return;
    brep.half_edges[he_a].prev = he_a_return;
    brep.half_edges[he_a_return].next = he_a;
    brep.half_edges[he_a_return].prev = he_a;
    brep.half_edges[he_b].next = he_b_return;
    brep.half_edges[he_b].prev = he_b_return;
    brep.half_edges[he_b_return].next = he_b;
    brep.half_edges[he_b_return].prev = he_b;
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

fn half_edge(origin: popkern_core::topology::VertexHandle, loop_: LoopHandle) -> HalfEdge {
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
