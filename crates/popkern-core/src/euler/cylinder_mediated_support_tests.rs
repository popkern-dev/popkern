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
