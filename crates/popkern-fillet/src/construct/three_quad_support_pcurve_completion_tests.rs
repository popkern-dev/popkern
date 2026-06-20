use std::collections::BTreeSet;

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregory};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, Edge, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle,
    Solid, Vertex, VertexHandle,
};

#[test]
fn rejected_cap_boundary_hint_falls_back_to_sample_repair() {
    let (mut brep, cap_face, cap_half_edge) = hinted_cap_support_fixture();
    let mut cap_faces = BTreeSet::new();
    cap_faces.insert(cap_face);

    let report =
        super::complete_missing_external_cap_support_pcurves(&mut brep, cap_face, &cap_faces);

    assert_eq!(report.eligible_edges, 1);
    assert_eq!(report.declined_cap_boundary_pcurves, 1);
    assert_eq!(report.repaired_cap_boundary_pcurves, 1);
    assert_eq!(report.authored_cap_boundary_pcurves, 0);

    let pcurve = brep.half_edges[cap_half_edge]
        .pcurve
        .as_ref()
        .expect("cap pcurve repaired from samples");
    let cap_surface = &brep.faces[cap_face].surface.as_ref().unwrap().geometry;
    let curve = brep.edges[brep.half_edges[cap_half_edge].edge]
        .curve
        .as_ref()
        .unwrap();
    let (t0, t1) = pcurve.parameter_range();
    let uv0 = pcurve.point_at(t0);
    let uv1 = pcurve.point_at(t1);
    assert!(
        uv0.x.abs() < 1.0e-9 && uv1.x.abs() < 1.0e-9,
        "repaired cap pcurve should follow the sampled UMin seam, got {uv0:?} -> {uv1:?}"
    );

    for index in 0..=4 {
        let s = index as f64 / 4.0;
        let uv = pcurve.point_at(t0 + s * (t1 - t0));
        let cap_point = cap_surface.point_at(uv.x, uv.y);
        let curve_point = curve
            .geometry
            .point_at(curve.t_start + s * (curve.t_end - curve.t_start));
        assert!(
            cap_point.distance_to(&curve_point) < 1.0e-8,
            "repaired cap pcurve should evaluate on the edge curve"
        );
    }
}

fn hinted_cap_support_fixture() -> (BRep, FaceHandle, HalfEdgeHandle) {
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
    let cap_face = brep.faces.insert(face(shell, gregory_surface(start, end)));
    let support_face = brep.faces.insert(face(shell, plane_surface()));
    let cap_loop = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: cap_face,
    });
    let support_loop = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: support_face,
    });
    let cap_half_edge = brep.half_edges.insert(half_edge(va, cap_loop));
    let support_half_edge = brep.half_edges.insert(half_edge(vb, support_loop));
    let edge = brep.edges.insert(Edge {
        half_edge: cap_half_edge,
        curve: CurveOnEdge::line_segment(start, end, 1.0e-12),
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });

    brep.half_edges[cap_half_edge].edge = edge;
    brep.half_edges[support_half_edge].edge = edge;
    brep.half_edges[cap_half_edge].twin = support_half_edge;
    brep.half_edges[support_half_edge].twin = cap_half_edge;
    brep.half_edges[cap_half_edge].next = support_half_edge;
    brep.half_edges[support_half_edge].next = cap_half_edge;
    brep.half_edges[cap_half_edge].prev = support_half_edge;
    brep.half_edges[support_half_edge].prev = cap_half_edge;
    brep.half_edges[cap_half_edge].pcurve = Some(Pcurve::linear((1.0, 0.0), (1.0, 1.0), 0.0, 1.0));
    brep.half_edges[support_half_edge].pcurve =
        Some(Pcurve::linear((0.0, 1.0), (0.0, 0.0), 0.0, 1.0));
    assert!(brep.set_surface_boundary_hint(cap_half_edge, BoundarySide::UMax));

    brep.loops[cap_loop].half_edge = cap_half_edge;
    brep.loops[support_loop].half_edge = support_half_edge;
    brep.faces[cap_face].outer_loop = cap_loop;
    brep.faces[support_face].outer_loop = support_loop;
    brep.vertices[va].outgoing = cap_half_edge;
    brep.vertices[vb].outgoing = support_half_edge;
    brep.shells[shell].face = cap_face;
    (brep, cap_face, cap_half_edge)
}

fn vertex(position: Point3) -> Vertex {
    Vertex {
        position,
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    }
}

fn half_edge(origin: VertexHandle, loop_: LoopHandle) -> HalfEdge {
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
