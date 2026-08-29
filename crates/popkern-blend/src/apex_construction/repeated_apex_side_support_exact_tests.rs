use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, Edge, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, Vertex, VertexHandle,
};

use super::repeated_apex_side_support::{
    SideSupportRegionContinuationOutcome, SideSupportRegionDiagnosticSnapshot,
    apply_selected_side_support_continuation,
};
use super::repeated_apex_trim_subcomplex::{
    TrimRewriteSubcomplexBuilder, TrimRewriteSubcomplexPlan,
};

#[test]
fn exact_authorable_merge_continuation_applies_curves_and_pcurves() {
    let (mut brep, plan) = exact_authorable_shared_edge_merge_fixture();
    let before_edges = brep.edges.len();
    let before_faces = brep.faces.len();
    let before_loops = brep.loops.len();

    let snapshot = SideSupportRegionDiagnosticSnapshot::build_read_only(&brep, &plan);

    let labels = snapshot.diagnostic_label();
    assert!(
        labels.contains("support-region-merge primitive_ready=true"),
        "{labels}"
    );
    assert!(labels.contains("exact_curve_authoring=ready"), "{labels}");
    assert!(labels.contains("execution ready_for_executor"), "{labels}");

    let applied = apply_selected_side_support_continuation(&mut brep, &plan)
        .expect("exact-authorable side-support merge applies");

    assert!(matches!(
        applied,
        SideSupportRegionContinuationOutcome::Applied(_)
    ));
    assert_eq!(brep.edges.len(), before_edges + 1);
    assert_eq!(brep.faces.len(), before_faces + 1);
    assert_eq!(brep.loops.len(), before_loops + 1);
    assert_eq!(
        brep.edges
            .values()
            .filter(|edge| edge.curve.is_some())
            .count(),
        2
    );
    assert_eq!(
        brep.half_edges
            .values()
            .filter(|half_edge| half_edge.pcurve.is_some())
            .count(),
        4
    );
    brep.validate_half_edge_pairing().unwrap();
    brep.validate_half_edge_cycles().unwrap();
    brep.validate_manifold_edges().unwrap();
}

fn exact_authorable_shared_edge_merge_fixture() -> (BRep, TrimRewriteSubcomplexPlan) {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    for face in brep.faces.values_mut() {
        face.surface = None;
    }
    for edge in brep.edges.values_mut() {
        edge.curve = None;
    }
    for half_edge in brep.half_edges.values_mut() {
        half_edge.pcurve = None;
    }

    let support_face = brep.faces.keys().next().expect("cube has a face");
    let shared_half_edge = first_shared_half_edge(&brep, support_face);
    let shared_twin = brep.half_edges[shared_half_edge].twin;
    let outside_face = brep.loops[brep.half_edges[shared_twin].loop_].face;
    let cylinder_face = add_isolated_cylinder_carrier_face(&mut brep, support_face);
    let support_start = brep.half_edges[shared_half_edge].origin;
    let support_end = brep.half_edges[shared_twin].origin;
    let outside_vertices = outer_loop_vertices(&brep, outside_face)
        .into_iter()
        .filter(|vertex| *vertex != support_start && *vertex != support_end)
        .collect::<Vec<_>>();
    assert_eq!(outside_vertices.len(), 2);

    assign_exact_merge_surfaces(&mut brep, support_face, outside_face, cylinder_face);

    let mut plan = TrimRewriteSubcomplexBuilder::new();
    plan.push_existing_seam(
        0,
        cylinder_face,
        outside_face,
        outside_vertices[0],
        outside_vertices[1],
    );
    plan.push_new_generator_trim(
        0,
        cylinder_face,
        support_face,
        support_start,
        outside_vertices[0],
    );
    plan.push_new_generator_trim(
        1,
        cylinder_face,
        support_face,
        outside_vertices[1],
        support_end,
    );
    (brep, plan.finish())
}

fn assign_exact_merge_surfaces(
    brep: &mut BRep,
    support_face: FaceHandle,
    outside_face: FaceHandle,
    cylinder_face: FaceHandle,
) {
    let outside_points = outer_loop_vertices(brep, outside_face)
        .into_iter()
        .map(|vertex| brep.vertices[vertex].position)
        .collect::<Vec<_>>();
    assert_eq!(outside_points.len(), 4);
    let plane = square_plane(&outside_points);
    brep.faces[support_face].surface = Some(SurfaceOnFace {
        geometry: plane.clone(),
        orientation: true,
    });
    brep.faces[outside_face].surface = Some(SurfaceOnFace {
        geometry: plane.clone(),
        orientation: true,
    });
    brep.faces[cylinder_face].surface = Some(SurfaceOnFace {
        geometry: square_vertex_cylinder(&outside_points, &plane),
        orientation: true,
    });
}

fn add_isolated_cylinder_carrier_face(brep: &mut BRep, support_face: FaceHandle) -> FaceHandle {
    let shell = brep.faces[support_face].shell;
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
    let vertices = [
        Point3::new(10.0, 10.0, 0.0),
        Point3::new(11.0, 10.0, 0.0),
        Point3::new(11.0, 11.0, 0.0),
        Point3::new(10.0, 11.0, 0.0),
    ]
    .map(|position| {
        brep.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        })
    });
    let half_edges = vertices.map(|vertex| {
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
            twin: HalfEdgeHandle::default(),
            next: Default::default(),
            prev: Default::default(),
            edge,
            loop_,
            pcurve: None,
        });
        brep.edges[edge].half_edge = half_edge;
        brep.vertices[vertex].outgoing = half_edge;
        brep.boundary_edges.insert(edge);
        half_edge
    });
    for i in 0..4 {
        brep.half_edges[half_edges[i]].next = half_edges[(i + 1) % 4];
        brep.half_edges[half_edges[i]].prev = half_edges[(i + 3) % 4];
    }
    brep.loops[loop_].half_edge = half_edges[0];
    face
}

fn first_shared_half_edge(brep: &BRep, face: FaceHandle) -> HalfEdgeHandle {
    let loop_ = brep.faces[face].outer_loop;
    let start = brep.loops[loop_].half_edge;
    let mut cursor = start;
    loop {
        let twin = brep.half_edges[cursor].twin;
        if brep.half_edges[twin].loop_ != loop_ {
            return cursor;
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            panic!("fixture face has a shared edge");
        }
    }
}

fn square_plane(points: &[Point3]) -> SurfaceGeometry {
    let u_axis = (points[1] - points[0])
        .normalized()
        .expect("fixture has non-degenerate square edge");
    let normal = (points[1] - points[0])
        .cross(&(points[2] - points[0]))
        .normalized()
        .expect("fixture has non-degenerate square plane");
    let v_axis = normal
        .cross(&u_axis)
        .normalized()
        .expect("fixture has non-degenerate square plane");
    SurfaceGeometry::Plane {
        origin: points[0],
        normal,
        u_axis,
        v_axis,
    }
}

fn square_vertex_cylinder(points: &[Point3], plane: &SurfaceGeometry) -> SurfaceGeometry {
    let SurfaceGeometry::Plane { normal, .. } = plane else {
        panic!("fixture plane expected");
    };
    let center = Point3::new(
        points.iter().map(|point| point.x).sum::<f64>() / points.len() as f64,
        points.iter().map(|point| point.y).sum::<f64>() / points.len() as f64,
        points.iter().map(|point| point.z).sum::<f64>() / points.len() as f64,
    );
    SurfaceGeometry::Cylinder {
        origin: center,
        axis: Vec3::new(normal.x, normal.y, normal.z),
        radius: center.distance_to(&points[0]),
    }
}

fn outer_loop_vertices(brep: &BRep, face: FaceHandle) -> Vec<VertexHandle> {
    let loop_handle = brep.faces[face].outer_loop;
    let start = brep.loops[loop_handle].half_edge;
    let mut vertices = Vec::new();
    let mut cursor = start;
    loop {
        let half_edge = &brep.half_edges[cursor];
        vertices.push(half_edge.origin);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    vertices
}
