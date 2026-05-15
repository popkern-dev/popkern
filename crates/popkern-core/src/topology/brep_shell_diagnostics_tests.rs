use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, Shell,
    ShellHandle, Solid, Vertex,
};

#[test]
fn diagnostics_accept_closed_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let report = analyze_brep_loop_shell_validity(&brep);

    assert!(report.has_no_errors());
    assert_eq!(report.shells, 1);
    assert_eq!(report.face_components, 1);
    assert_eq!(report.boundary_edges, 0);
    assert_eq!(report.boundary_cycles, 0);
    assert_eq!(report.boundary_chains, 0);
}

#[test]
fn diagnostics_reports_open_box_boundary_cycle() {
    let brep = BRep::open_box(Tolerance::default()).unwrap();
    let report = analyze_brep_loop_shell_validity(&brep);

    assert!(report.has_no_errors());
    assert_eq!(report.boundary_edges, 4);
    assert_eq!(report.boundary_cycles, 1);
    assert_eq!(report.boundary_chains, 0);
    assert_eq!(report.face_components, 1);
}

#[test]
fn diagnostics_reports_fragmented_shell_components() {
    let brep = two_square_fragmented_shell();
    let report = analyze_brep_loop_shell_validity(&brep);

    assert!(report.has_no_errors());
    assert_eq!(report.face_components, 2);
    assert_eq!(report.boundary_cycles, 2);
    assert!(report.findings.iter().any(|finding| {
        finding.kind == BrepTopologyDiagnosticKind::FragmentedShell
            && finding.severity == BrepTopologyDiagnosticSeverity::Warning
    }));
}

#[test]
fn diagnostics_reports_corrupted_loop_and_radial_entry() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let loop_handle = brep.loops.keys().next().unwrap();
    brep.loops[loop_handle].half_edge = HalfEdgeHandle::default();
    brep.nonmanifold_radial
        .insert(HalfEdgeHandle::default(), HalfEdgeHandle::default());

    let report = analyze_brep_loop_shell_validity(&brep);

    assert!(!report.has_no_errors());
    assert!(report.findings.iter().any(|finding| {
        finding.kind == BrepTopologyDiagnosticKind::MissingHandle
            && finding.loop_ == Some(loop_handle)
    }));
    assert!(report.findings.iter().any(|finding| {
        finding.kind == BrepTopologyDiagnosticKind::BrokenLoopCycle
            && finding.loop_ == Some(loop_handle)
    }));
    assert!(report.findings.iter().any(|finding| {
        finding.kind == BrepTopologyDiagnosticKind::RadialConflict
            && finding.severity == BrepTopologyDiagnosticSeverity::Error
    }));
}

fn two_square_fragmented_shell() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    brep.solids[solid].shell = shell;

    let first = add_square(
        &mut brep,
        shell,
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
    );
    let _second = add_square(
        &mut brep,
        shell,
        [
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ],
    );
    brep.shells[shell].face = first.face;
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
    for index in 0..4 {
        let edge = brep.edges.insert(Edge {
            half_edge: half_edges[index],
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        brep.half_edges[half_edges[index]].edge = edge;
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
