//! PK-509 negative controls for `validate_all`: each targeted topology
//! corruption of a valid, geometry-bearing cube is rejected with its
//! expected structured error. The positive control (the untouched cube)
//! validates clean, so every rejection is caused by the corruption.

use crate::error::ValidationError;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).unwrap_or_else(|e| panic!("cube: {e:?}"))
}

fn variants(errors: &[ValidationError]) -> Vec<&'static str> {
    errors
        .iter()
        .map(|e| match e {
            ValidationError::DanglingReference { .. } => "DanglingReference",
            ValidationError::EulerPoincare { .. } => "EulerPoincare",
            ValidationError::BrokenTwin { .. } => "BrokenTwin",
            ValidationError::BrokenCycle { .. } => "BrokenCycle",
            ValidationError::LoopJunctionMismatch { .. } => "LoopJunctionMismatch",
            ValidationError::InconsistentVertex { .. } => "InconsistentVertex",
            ValidationError::NonManifoldEdge { .. } => "NonManifoldEdge",
            ValidationError::InconsistentOrientation { .. } => "InconsistentOrientation",
            ValidationError::OpenLoop { .. } => "OpenLoop",
            ValidationError::CurveSurfaceDeviation { .. } => "CurveSurfaceDeviation",
            ValidationError::ShellGap { .. } => "ShellGap",
            ValidationError::WindingDishonest { .. } => "WindingDishonest",
            _ => "Other",
        })
        .collect()
}

#[test]
fn the_untouched_cube_validates_clean() {
    assert!(cube().validate_all().is_empty());
}

#[test]
fn a_flipped_coedge_is_rejected() {
    let mut b = cube();
    let he = b
        .half_edges
        .keys()
        .next()
        .unwrap_or_else(|| panic!("half-edge"));
    let twin = b.half_edges[he].twin;
    b.half_edges[he].origin = b.half_edges[twin].origin;
    let found = variants(&b.validate_all());
    assert!(found.contains(&"InconsistentVertex"), "{found:?}");
}

#[test]
fn a_permuted_loop_is_rejected() {
    // One face's loop a -> b -> c -> d rewired a -> c -> b -> d: it
    // still cycles through all four and closes, but jumps at its
    // junctions (PK-943).
    let mut b = cube();
    let f = b.faces.keys().next().unwrap_or_else(|| panic!("face"));
    let a = b.loops[b.faces[f].outer_loop].half_edge;
    let hb = b.half_edges[a].next;
    let hc = b.half_edges[hb].next;
    let hd = b.half_edges[hc].next;
    b.half_edges[a].next = hc;
    b.half_edges[hc].next = hb;
    b.half_edges[hb].next = hd;
    b.half_edges[hc].prev = a;
    b.half_edges[hb].prev = hc;
    b.half_edges[hd].prev = hb;
    assert!(
        b.validate_half_edge_cycles().is_ok(),
        "the cycle still closes"
    );
    let found = variants(&b.validate_all());
    assert!(found.contains(&"LoopJunctionMismatch"), "{found:?}");
    assert!(matches!(
        b.validate(),
        Err(ValidationError::LoopJunctionMismatch { .. })
    ));
}

#[test]
fn a_dangling_edge_is_rejected() {
    let mut b = cube();
    let e = b.edges.keys().next().unwrap_or_else(|| panic!("edge"));
    let he = b.edges[e].half_edge;
    b.half_edges.remove(he);
    let found = variants(&b.validate_all());
    assert!(found.contains(&"BrokenTwin"), "{found:?}");
    assert!(found.contains(&"BrokenCycle"), "{found:?}");
}

#[test]
fn a_wrong_face_orientation_is_rejected() {
    let mut b = cube();
    let f = b.faces.keys().next().unwrap_or_else(|| panic!("face"));
    if let Some(s) = b.faces[f].surface.as_mut() {
        s.orientation = !s.orientation;
    }
    let found = variants(&b.validate_all());
    assert!(found.contains(&"WindingDishonest"), "{found:?}");
}

/// The genus-aware Euler check tolerates an odd imbalance of one, so an
/// orphan vertex passed `validate_all` until the vertex-anchored
/// reference check (PK-509).
#[test]
fn an_euler_break_by_an_orphan_vertex_is_rejected() {
    let mut b = cube();
    let v = b.vertices.keys().next().unwrap_or_else(|| panic!("vertex"));
    let position = b.vertices[v].position;
    b.vertices.insert(crate::topology::Vertex {
        position,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let found = variants(&b.validate_all());
    assert!(found.contains(&"DanglingReference"), "{found:?}");
}

/// An Euler break by a removed face is rejected too.
#[test]
fn an_euler_break_by_a_removed_face_is_rejected() {
    let mut b = cube();
    let f = b.faces.keys().next().unwrap_or_else(|| panic!("face"));
    b.faces.remove(f);
    assert!(!b.validate_all().is_empty());
}
