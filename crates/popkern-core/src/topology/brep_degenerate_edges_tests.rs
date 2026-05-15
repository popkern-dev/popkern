use super::*;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::{Point3, SurfaceGeometry, SurfaceOnFace, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, Edge, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell, ShellHandle,
    Solid, Vertex,
};

fn options() -> BrepDegenerateEdgeOptions {
    BrepDegenerateEdgeOptions {
        tolerance: 1e-6,
        samples: 5,
    }
}

#[test]
fn classifies_collapsed_surface_trim_without_mutation() {
    let brep = single_edge_sheet(true, true, true);
    let report = analyze_brep_degenerate_edges(&brep, options());

    assert_eq!(report.degenerate_edges, 1);
    assert_eq!(report.collapsed_surface_trims, 1);
    assert_eq!(report.zero_length_vertex_edges, 0);
    assert_eq!(report.ambiguous_edges, 0);
    let finding = &report.findings[0];
    assert_eq!(finding.kind, BrepDegenerateEdgeKind::CollapsedSurfaceTrim);
    assert_eq!(finding.source_trims, 1);
    assert!(finding.surface_span <= 1e-6);
    assert!(finding.pcurve_span > 0.0);
}

#[test]
fn classifies_zero_length_vertex_edge_without_surface_evidence() {
    let brep = single_edge_sheet(false, false, true);
    let report = analyze_brep_degenerate_edges(&brep, options());

    assert_eq!(report.degenerate_edges, 1);
    assert_eq!(report.collapsed_surface_trims, 0);
    assert_eq!(report.zero_length_vertex_edges, 1);
    assert_eq!(
        report.findings[0].kind,
        BrepDegenerateEdgeKind::ZeroLengthVertexEdge
    );
}

#[test]
fn classifies_ambiguous_trim_when_surface_span_is_not_collapsed() {
    let brep = single_edge_sheet(true, true, false);
    let report = analyze_brep_degenerate_edges(&brep, options());

    assert_eq!(report.degenerate_edges, 1);
    assert_eq!(report.ambiguous_edges, 1);
    assert_eq!(
        report.findings[0].kind,
        BrepDegenerateEdgeKind::AmbiguousCollapsedTrim
    );
    assert!(report.findings[0].surface_span > 1e-6);
}

fn single_edge_sheet(with_surface: bool, with_pcurve: bool, collapsed_surface: bool) -> BRep {
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
        surface: with_surface.then_some(surface(collapsed_surface)),
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
        edge: Default::default(),
        loop_,
        pcurve: with_pcurve.then_some(Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0)),
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

fn surface(collapsed: bool) -> SurfaceOnFace {
    let u_axis = if collapsed {
        Vec3::new(0.0, 0.0, 0.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis,
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    }
}
