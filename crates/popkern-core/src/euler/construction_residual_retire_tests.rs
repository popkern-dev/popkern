use crate::tolerance::Tolerance;
use crate::topology::{BRep, TopologyRole};

#[test]
fn retire_construction_residual_face_exposes_surviving_edges() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let face = brep.faces.keys().next().expect("cube has a face");
    brep.faces[face].role = TopologyRole::ConstructionResidual;

    let result = brep
        .retire_construction_residual_faces(&[face])
        .expect("residual face retires");

    assert_eq!(result.retired_faces, vec![face]);
    assert_eq!(result.exposed_edges.len(), 4);
    assert_eq!(brep.faces.len(), 5);
    assert_eq!(brep.loops.len(), 5);
    assert_eq!(brep.boundary_edges.len(), 4);
    brep.validate()
        .expect("retired residual leaves valid open shell");
}

#[test]
fn retire_construction_residual_rejects_production_face() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let face = brep.faces.keys().next().expect("cube has a face");
    let before = (brep.faces.len(), brep.loops.len(), brep.half_edges.len());

    let err = brep
        .retire_construction_residual_faces(&[face])
        .expect_err("production face is protected");

    assert!(err.to_string().contains("construction-residual"));
    assert_eq!(
        before,
        (brep.faces.len(), brep.loops.len(), brep.half_edges.len())
    );
}

#[test]
fn retire_unwired_construction_residual_edges_removes_null_loop_half_edges() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let vertex = brep.vertices.keys().next().expect("cube has a vertex");
    let edge = brep.edges.insert(crate::topology::Edge {
        half_edge: Default::default(),
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: TopologyRole::ConstructionResidual,
    });
    let half_edge = brep.half_edges.insert(crate::topology::HalfEdge {
        origin: vertex,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge,
        loop_: Default::default(),
        pcurve: None,
    });
    brep.edges[edge].half_edge = half_edge;
    brep.boundary_edges.insert(edge);

    let result = brep
        .retire_unwired_construction_residual_edges()
        .expect("unwired residual retires");

    assert_eq!(result.removed_edges, vec![edge]);
    assert_eq!(result.removed_half_edges, vec![half_edge]);
    assert!(!brep.edges.contains_key(edge));
    assert!(!brep.half_edges.contains_key(half_edge));
    brep.validate()
        .expect("unwired residual retire preserves valid production topology");
}
