use crate::constructors::builder::BRepBuilder;
use crate::euler::{
    CylinderMediatedSupportContinuationDelta, CylinderMediatedSupportContinuationRequest,
};
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

#[test]
fn preflight_accepts_cylinder_mediated_operands() {
    let fixture = cylinder_mediated_fixture();

    let plan = fixture
        .brep
        .preflight_cylinder_mediated_support_continuation(fixture.request())
        .expect("valid operands preflight");

    assert_eq!(
        plan.expected_delta,
        CylinderMediatedSupportContinuationDelta::Unknown {
            reason: "cylinder-mediated-support-continuation-operator-required"
        }
    );
    assert_eq!(
        plan.outside_interval,
        vec![
            fixture.outside_cylinder_half_edge,
            fixture.outside_end_half_edge
        ]
    );
    assert_eq!(plan.outside_complement_interval.len(), 4);
}

#[test]
fn preflight_refuses_misdirected_cylinder_interval() {
    let fixture = cylinder_mediated_fixture();
    let mut request = fixture.request();
    request.cylinder_interval.reverse();

    let err = fixture
        .brep
        .preflight_cylinder_mediated_support_continuation(request)
        .expect_err("reversed interval is not a valid directed operand");

    assert!(err.to_string().contains("endpoints mismatch"), "{err}");
}

#[test]
fn preflight_refuses_wrong_support_cylinder_adjacency() {
    let fixture = cylinder_mediated_fixture();
    let mut request = fixture.request();
    request.support_cylinder_half_edge = fixture.outside_cylinder_half_edge;

    let err = fixture
        .brep
        .preflight_cylinder_mediated_support_continuation(request)
        .expect_err("outside half-edge is not on support loop");

    assert!(
        err.to_string().contains("support/cylinder adjacency"),
        "{err}"
    );
}

#[test]
fn preflight_refuses_outside_adjacency_off_endpoint_intervals() {
    let fixture = cylinder_mediated_fixture();
    let mut request = fixture.request();
    request.outside_cylinder_half_edge = fixture.support_cylinder_half_edge;

    let err = fixture
        .brep
        .preflight_cylinder_mediated_support_continuation(request)
        .expect_err("outside adjacency must belong to outside endpoint interval");

    assert!(
        err.to_string().contains("outside/cylinder adjacency"),
        "{err}"
    );
}

#[test]
fn preflight_accepts_cube_extreme_style_long_interval_requests() {
    let fixture = long_interval_cylinder_mediated_fixture();

    let first = fixture
        .brep
        .preflight_cylinder_mediated_support_continuation(fixture.first_request())
        .expect("first long-interval request preflights");
    let second = fixture
        .brep
        .preflight_cylinder_mediated_support_continuation(fixture.second_request())
        .expect("second long-interval request preflights");

    assert_eq!(first.outside_interval.len(), 13);
    assert_eq!(first.outside_complement_interval.len(), 9);
    assert_eq!(first.request.cylinder_interval.len(), 5);
    assert_eq!(second.outside_interval.len(), 11);
    assert_eq!(second.outside_complement_interval.len(), 11);
    assert_eq!(second.request.cylinder_interval.len(), 2);
}

struct CylinderMediatedFixture {
    brep: BRep,
    support_face: FaceHandle,
    outside_face: FaceHandle,
    cylinder_face: FaceHandle,
    support_loop: LoopHandle,
    outside_loop: LoopHandle,
    start: VertexHandle,
    end: VertexHandle,
    outside_start_half_edge: HalfEdgeHandle,
    outside_end_half_edge: HalfEdgeHandle,
    support_cylinder_half_edge: HalfEdgeHandle,
    outside_cylinder_half_edge: HalfEdgeHandle,
    cylinder_interval: Vec<HalfEdgeHandle>,
}

struct LongIntervalCylinderMediatedFixture {
    brep: BRep,
    outside_face: FaceHandle,
    cylinder_face: FaceHandle,
    outside_loop: LoopHandle,
    outside_half_edges: Vec<HalfEdgeHandle>,
    support_faces: [FaceHandle; 2],
    support_loops: [LoopHandle; 2],
    support_cylinder_half_edges: [HalfEdgeHandle; 2],
    outside_cylinder_half_edges: [HalfEdgeHandle; 2],
    cylinder_intervals: [Vec<HalfEdgeHandle>; 2],
}

impl LongIntervalCylinderMediatedFixture {
    fn first_request(&self) -> CylinderMediatedSupportContinuationRequest {
        self.request(0, 0, 12)
    }

    fn second_request(&self) -> CylinderMediatedSupportContinuationRequest {
        self.request(1, 5, 15)
    }

    fn request(
        &self,
        support_index: usize,
        start_index: usize,
        end_index: usize,
    ) -> CylinderMediatedSupportContinuationRequest {
        CylinderMediatedSupportContinuationRequest {
            support_face: self.support_faces[support_index],
            outside_face: self.outside_face,
            cylinder_face: self.cylinder_face,
            start: self.brep.half_edges[self.outside_half_edges[start_index]].origin,
            end: self.brep.half_edges[self.outside_half_edges[end_index]].origin,
            outside_loop: self.outside_loop,
            outside_start_half_edge: self.outside_half_edges[start_index],
            outside_end_half_edge: self.outside_half_edges[end_index],
            support_loop: self.support_loops[support_index],
            support_cylinder_half_edge: self.support_cylinder_half_edges[support_index],
            outside_cylinder_half_edge: self.outside_cylinder_half_edges[support_index],
            cylinder_interval: self.cylinder_intervals[support_index].clone(),
        }
    }
}

impl CylinderMediatedFixture {
    fn request(&self) -> CylinderMediatedSupportContinuationRequest {
        CylinderMediatedSupportContinuationRequest {
            support_face: self.support_face,
            outside_face: self.outside_face,
            cylinder_face: self.cylinder_face,
            start: self.start,
            end: self.end,
            outside_loop: self.outside_loop,
            outside_start_half_edge: self.outside_start_half_edge,
            outside_end_half_edge: self.outside_end_half_edge,
            support_loop: self.support_loop,
            support_cylinder_half_edge: self.support_cylinder_half_edge,
            outside_cylinder_half_edge: self.outside_cylinder_half_edge,
            cylinder_interval: self.cylinder_interval.clone(),
        }
    }
}

fn cylinder_mediated_fixture() -> CylinderMediatedFixture {
    let mut builder = BRepBuilder::new(Tolerance::default());
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);
    let vertices = [
        builder.add_vertex(Point3::new(0.0, 0.0, 0.0)),
        builder.add_vertex(Point3::new(1.0, 0.0, 0.0)),
        builder.add_vertex(Point3::new(1.0, 1.0, 0.0)),
        builder.add_vertex(Point3::new(0.0, 1.0, 0.0)),
        builder.add_vertex(Point3::new(2.0, 0.0, 0.0)),
        builder.add_vertex(Point3::new(3.0, 0.0, 0.0)),
        builder.add_vertex(Point3::new(3.0, 1.0, 0.0)),
        builder.add_vertex(Point3::new(2.0, 1.0, 0.0)),
    ];
    let [a, b, c, d, e, f, g, h] = vertices;
    let support_face = builder.add_face(shell);
    let cylinder_face = builder.add_face(shell);
    let outside_face = builder.add_face(shell);

    let (_, support_cylinder_half_edge, support_cylinder_twin) = builder.add_edge(a, b);
    let (_, support_side_1, _) = builder.add_edge(b, c);
    let (_, support_side_2, _) = builder.add_edge(c, d);
    let (_, support_side_3, _) = builder.add_edge(d, a);
    let support_loop = builder
        .add_loop(
            &[
                support_cylinder_half_edge,
                support_side_1,
                support_side_2,
                support_side_3,
            ],
            support_face,
        )
        .expect("support loop wires");

    let (_, outside_cylinder_half_edge, outside_cylinder_twin) = builder.add_edge(e, f);
    let (_, outside_end_half_edge, _) = builder.add_edge(f, g);
    let (_, outside_side_2, _) = builder.add_edge(g, h);
    let (_, outside_side_3, _) = builder.add_edge(h, e);
    let outside_loop = builder
        .add_loop(
            &[
                outside_cylinder_half_edge,
                outside_end_half_edge,
                outside_side_2,
                outside_side_3,
            ],
            outside_face,
        )
        .expect("outside loop wires");

    let (_, cylinder_mid_0, _) = builder.add_edge(b, e);
    let (_, cylinder_mid_1, _) = builder.add_edge(f, a);
    let cylinder_loop = [
        outside_cylinder_twin,
        cylinder_mid_0,
        support_cylinder_twin,
        cylinder_mid_1,
    ];
    builder
        .add_loop(&cylinder_loop, cylinder_face)
        .expect("cylinder loop wires");

    let brep = builder.build_unchecked();
    CylinderMediatedFixture {
        brep,
        support_face,
        outside_face,
        cylinder_face,
        support_loop,
        outside_loop,
        start: e,
        end: f,
        outside_start_half_edge: outside_cylinder_half_edge,
        outside_end_half_edge,
        support_cylinder_half_edge,
        outside_cylinder_half_edge,
        cylinder_interval: vec![outside_cylinder_twin, cylinder_mid_0, support_cylinder_twin],
    }
}

fn long_interval_cylinder_mediated_fixture() -> LongIntervalCylinderMediatedFixture {
    let mut builder = BRepBuilder::new(Tolerance::default());
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);
    let outside_face = builder.add_face(shell);
    let cylinder_face = builder.add_face(shell);
    let first_support_face = builder.add_face(shell);
    let second_support_face = builder.add_face(shell);

    let outside_vertices = (0..20)
        .map(|index| builder.add_vertex(Point3::new(index as f64, 0.0, 0.0)))
        .collect::<Vec<_>>();
    let mut outside_half_edges = Vec::new();
    let mut outside_twins = Vec::new();
    for index in 0..outside_vertices.len() {
        let next = (index + 1) % outside_vertices.len();
        let (_, half_edge, twin) =
            builder.add_edge(outside_vertices[index], outside_vertices[next]);
        outside_half_edges.push(half_edge);
        outside_twins.push(twin);
    }
    let outside_loop = builder
        .add_loop(&outside_half_edges, outside_face)
        .expect("outside loop wires");

    let (first_support_loop, first_support_cylinder_half_edge, first_support_cylinder_twin) =
        add_support_square(&mut builder, first_support_face, 100.0);
    let (second_support_loop, second_support_cylinder_half_edge, second_support_cylinder_twin) =
        add_support_square(&mut builder, second_support_face, 200.0);

    let cylinder_connectors = (0..8)
        .map(|index| {
            let v0 = builder.add_vertex(Point3::new(300.0 + index as f64, 0.0, 0.0));
            let v1 = builder.add_vertex(Point3::new(300.5 + index as f64, 0.0, 0.0));
            let (_, half_edge, _) = builder.add_edge(v0, v1);
            half_edge
        })
        .collect::<Vec<_>>();
    let cylinder_loop_half_edges = vec![
        outside_twins[3],
        cylinder_connectors[0],
        cylinder_connectors[1],
        cylinder_connectors[2],
        first_support_cylinder_twin,
        cylinder_connectors[3],
        outside_twins[7],
        second_support_cylinder_twin,
        cylinder_connectors[4],
        cylinder_connectors[5],
        cylinder_connectors[6],
        cylinder_connectors[7],
    ];
    builder
        .add_loop(&cylinder_loop_half_edges, cylinder_face)
        .expect("cylinder loop wires");

    let brep = builder.build_unchecked();
    LongIntervalCylinderMediatedFixture {
        brep,
        outside_face,
        cylinder_face,
        outside_loop,
        outside_half_edges: outside_half_edges.clone(),
        support_faces: [first_support_face, second_support_face],
        support_loops: [first_support_loop, second_support_loop],
        support_cylinder_half_edges: [
            first_support_cylinder_half_edge,
            second_support_cylinder_half_edge,
        ],
        outside_cylinder_half_edges: [outside_half_edges[3], outside_half_edges[7]],
        cylinder_intervals: [
            vec![
                outside_twins[3],
                cylinder_connectors[0],
                cylinder_connectors[1],
                cylinder_connectors[2],
                first_support_cylinder_twin,
            ],
            vec![outside_twins[7], second_support_cylinder_twin],
        ],
    }
}

fn add_support_square(
    builder: &mut BRepBuilder,
    support_face: FaceHandle,
    offset: f64,
) -> (LoopHandle, HalfEdgeHandle, HalfEdgeHandle) {
    let vertices = [
        builder.add_vertex(Point3::new(offset, 0.0, 0.0)),
        builder.add_vertex(Point3::new(offset + 1.0, 0.0, 0.0)),
        builder.add_vertex(Point3::new(offset + 1.0, 1.0, 0.0)),
        builder.add_vertex(Point3::new(offset, 1.0, 0.0)),
    ];
    let (_, support_cylinder_half_edge, support_cylinder_twin) =
        builder.add_edge(vertices[0], vertices[1]);
    let (_, side_1, _) = builder.add_edge(vertices[1], vertices[2]);
    let (_, side_2, _) = builder.add_edge(vertices[2], vertices[3]);
    let (_, side_3, _) = builder.add_edge(vertices[3], vertices[0]);
    let support_loop = builder
        .add_loop(
            &[support_cylinder_half_edge, side_1, side_2, side_3],
            support_face,
        )
        .expect("support loop wires");
    (
        support_loop,
        support_cylinder_half_edge,
        support_cylinder_twin,
    )
}
