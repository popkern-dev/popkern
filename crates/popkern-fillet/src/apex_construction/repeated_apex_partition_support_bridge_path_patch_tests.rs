use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

use super::*;

#[test]
fn bridge_boundary_path_detects_multi_loop_regrouping_need() {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: Default::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let vertices = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ]
    .map(|position| {
        brep.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        })
    });
    let first_loop = add_boundary_loop(&mut brep, shell, &[vertices[0], vertices[1], vertices[3]]);
    let second_loop = add_boundary_loop(&mut brep, shell, &[vertices[1], vertices[2], vertices[3]]);

    let path = bridge_boundary_path_for_cylinder_side(
        &brep,
        &[first_loop, second_loop],
        vertices[0],
        vertices[2],
    )
    .expect("multi-loop boundary path is detected");

    assert_eq!(path.len(), 2);
    assert_eq!(path[0].start, vertices[0]);
    assert_eq!(path[0].end, vertices[1]);
    assert_eq!(path[1].start, vertices[1]);
    assert_eq!(path[1].end, vertices[2]);
}

#[test]
fn bridge_path_patch_authors_support_owned_trim_without_rewiring_source_loops() {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: Default::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let vertices = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ]
    .map(|position| {
        brep.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        })
    });
    let support_face =
        add_boundary_loop_face(&mut brep, shell, &[vertices[0], vertices[2], vertices[3]]);
    let first_loop = add_boundary_loop(&mut brep, shell, &[vertices[0], vertices[1], vertices[3]]);
    let second_loop = add_boundary_loop(&mut brep, shell, &[vertices[1], vertices[2], vertices[3]]);
    let first_loop_before = collect_loop_half_edges(&brep, first_loop).expect("first loop");
    let second_loop_before = collect_loop_half_edges(&brep, second_loop).expect("second loop");
    let curve = CurveOnEdge::line_segment(
        brep.vertices[vertices[0]].position,
        brep.vertices[vertices[2]].position,
        brep.tolerance.linear(),
    )
    .expect("nondegenerate trim curve");

    let (cylinder_half_edge, support_patch_face) =
        author_side_support_bridge_path_patch_for_cylinder_trim(
            &mut brep,
            SideSupportBridgePathPatchRequest {
                support_face,
                bridge_loops: &[first_loop, second_loop],
                start: vertices[0],
                end: vertices[2],
                curve,
                label: "test generator",
                fillet_index: 0,
            },
        )
        .expect("patch authoring succeeds")
        .expect("path patch applies");

    let support_half_edge = brep.half_edges[cylinder_half_edge].twin;
    let support_loop = brep.half_edges[support_half_edge].loop_;
    assert_eq!(brep.loops[support_loop].face, support_patch_face);
    assert_eq!(
        half_edge_end_vertex(&brep, cylinder_half_edge),
        Some(vertices[2])
    );
    assert_eq!(
        support_split_half_edge_end_vertex(&brep, support_half_edge),
        Some(vertices[0])
    );
    assert_eq!(
        collect_loop_half_edges(&brep, first_loop).expect("first loop after"),
        first_loop_before
    );
    assert_eq!(
        collect_loop_half_edges(&brep, second_loop).expect("second loop after"),
        second_loop_before
    );
    validate_replacement_loop_support_side_ownership(
        &brep,
        FaceHandle::default(),
        0,
        &[cylinder_half_edge],
    )
    .expect("cylinder trim has a live non-cylinder support side");
}

#[test]
fn bridge_path_patch_accepts_single_span_boundary_without_rewiring_source_loop() {
    let mut brep = BRep::new(Tolerance::default());
    let solid = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell = brep.shells.insert(Shell {
        face: Default::default(),
        solid,
    });
    brep.solids[solid].shell = shell;
    let vertices = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ]
    .map(|position| {
        brep.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        })
    });
    let support_face =
        add_boundary_loop_face(&mut brep, shell, &[vertices[0], vertices[1], vertices[2]]);
    let bridge_loop = add_boundary_loop(&mut brep, shell, &[vertices[0], vertices[1], vertices[2]]);
    let bridge_loop_before = collect_loop_half_edges(&brep, bridge_loop).expect("bridge loop");
    let curve = CurveOnEdge::line_segment(
        brep.vertices[vertices[0]].position,
        brep.vertices[vertices[1]].position,
        brep.tolerance.linear(),
    )
    .expect("nondegenerate trim curve");

    let (cylinder_half_edge, support_patch_face) =
        author_side_support_bridge_path_patch_for_cylinder_trim(
            &mut brep,
            SideSupportBridgePathPatchRequest {
                support_face,
                bridge_loops: &[bridge_loop],
                start: vertices[0],
                end: vertices[1],
                curve,
                label: "single span generator",
                fillet_index: 0,
            },
        )
        .expect("patch authoring succeeds")
        .expect("single span path patch applies");

    let support_half_edge = brep.half_edges[cylinder_half_edge].twin;
    assert_eq!(
        half_edge_end_vertex(&brep, cylinder_half_edge),
        Some(vertices[1])
    );
    assert_eq!(
        support_split_half_edge_end_vertex(&brep, support_half_edge),
        Some(vertices[0])
    );
    assert_eq!(
        brep.loops[brep.half_edges[support_half_edge].loop_].face,
        support_patch_face
    );
    assert_eq!(
        collect_loop_half_edges(&brep, bridge_loop).expect("bridge loop after"),
        bridge_loop_before
    );
    validate_replacement_loop_support_side_ownership(
        &brep,
        FaceHandle::default(),
        0,
        &[cylinder_half_edge],
    )
    .expect("single-span cylinder trim has a live non-cylinder support side");
}

fn add_boundary_loop(
    brep: &mut BRep,
    shell: popkern_core::topology::ShellHandle,
    vertices: &[VertexHandle],
) -> LoopHandle {
    let face = add_boundary_loop_face(brep, shell, vertices);
    brep.faces[face].outer_loop
}

fn add_boundary_loop_face(
    brep: &mut BRep,
    shell: popkern_core::topology::ShellHandle,
    vertices: &[VertexHandle],
) -> FaceHandle {
    let face = brep.faces.insert(Face {
        outer_loop: Default::default(),
        inner_loops: vec![],
        shell,
        surface: None,
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    let loop_ = brep.loops.insert(Loop {
        half_edge: Default::default(),
        face,
    });
    brep.faces[face].outer_loop = loop_;
    let half_edges = vertices
        .iter()
        .copied()
        .map(|vertex| {
            let edge = brep.edges.insert(Edge {
                half_edge: Default::default(),
                curve: None,
                name: None,
                tolerance: None,
                canonical: Default::default(),
                role: Default::default(),
            });
            let half_edge = brep.half_edges.insert(HalfEdge {
                origin: vertex,
                twin: Default::default(),
                next: Default::default(),
                prev: Default::default(),
                edge,
                loop_,
                pcurve: None,
            });
            brep.edges[edge].half_edge = half_edge;
            brep.boundary_edges.insert(edge);
            if brep.vertices[vertex].outgoing == HalfEdgeHandle::default() {
                brep.vertices[vertex].outgoing = half_edge;
            }
            half_edge
        })
        .collect::<Vec<_>>();
    for (index, half_edge) in half_edges.iter().copied().enumerate() {
        brep.half_edges[half_edge].next = half_edges[(index + 1) % half_edges.len()];
        brep.half_edges[half_edge].prev =
            half_edges[(index + half_edges.len() - 1) % half_edges.len()];
    }
    brep.loops[loop_].half_edge = half_edges[0];
    face
}
