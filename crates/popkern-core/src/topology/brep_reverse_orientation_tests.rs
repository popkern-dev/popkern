//! Tests for whole-body orientation reversal.

use crate::geometry::curve::CurveEval;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn dishonest_count(brep: &BRep) -> usize {
    brep.audit_winding_honesty()
        .iter()
        .filter(|v| v.honest == Some(false))
        .count()
}

#[test]
fn reversal_preserves_winding_honesty() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    assert_eq!(dishonest_count(&brep), 0, "cube starts honest");
    brep.reverse_body_orientation();
    assert_eq!(
        dishonest_count(&brep),
        0,
        "reversed body must stay winding-honest (winding and flag negate together)"
    );
    assert!(
        brep.faces
            .values()
            .all(|f| f.surface.as_ref().is_some_and(|s| !s.orientation)),
        "all cube flags flip false"
    );
}

#[test]
fn reversal_is_an_involution() {
    let original = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut brep = original.clone();
    brep.reverse_body_orientation();
    brep.reverse_body_orientation();
    for (handle, he) in &original.half_edges {
        let round = &brep.half_edges[handle];
        assert_eq!(round.origin, he.origin);
        assert_eq!(round.next, he.next);
        assert_eq!(round.prev, he.prev);
        assert_eq!(
            round.mate().unwrap_or_default(),
            he.mate().unwrap_or_default()
        );
        assert_eq!(round.loop_, he.loop_);
    }
    for (vh, v) in &original.vertices {
        assert_eq!(brep.vertices[vh].outgoing, v.outgoing);
    }
    for (fh, f) in &original.faces {
        assert_eq!(
            brep.faces[fh].surface.as_ref().map(|s| s.orientation),
            f.surface.as_ref().map(|s| s.orientation)
        );
    }
}

#[test]
fn reversal_keeps_structural_validity() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    brep.reverse_body_orientation();
    // Twins stay opposite: each half-edge's origin equals its twin's
    // destination (= twin.twin.origin's neighbour via next chain).
    for (_, he) in &brep.half_edges {
        let twin = &brep.half_edges[he.mate().unwrap_or_default()];
        assert_ne!(
            he.origin, twin.origin,
            "twin halves must originate at opposite ends"
        );
    }
    // Loops still close and traverse continuously.
    for (_, loop_data) in &brep.loops {
        let start = loop_data.half_edge;
        let mut cursor = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let he = &brep.half_edges[cursor];
            let next = &brep.half_edges[he.next];
            // Continuity: this half-edge's destination is the next's origin.
            assert_eq!(
                brep.half_edges[he.mate().unwrap_or_default()].origin,
                next.origin,
                "loop traversal must stay vertex-continuous"
            );
            assert_eq!(brep.half_edges[next.prev].origin, he.origin, "prev inverse");
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
        assert_eq!(cursor, start, "loop must close");
    }
    // Vertex outgoing pointers anchor at their own vertex.
    for (vh, vertex) in &brep.vertices {
        assert_eq!(
            brep.half_edges[vertex.outgoing].origin, vh,
            "outgoing half-edge must originate at its vertex"
        );
    }
}

/// A sheet reverses like a solid: its free fins turn with the rest (the
/// mate-based reversal left them unturned, pointed their vertices'
/// outgoing at nothing and dropped the free edges' canonical fin), and
/// reversing twice restores it.
#[test]
fn reversal_turns_a_sheet_with_its_free_edges() {
    let original = BRep::open_box(Tolerance::default()).unwrap();
    let mut brep = original.clone();
    brep.reverse_body_orientation();
    brep.validate().expect("reversed open box validates");
    for (vh, vertex) in &brep.vertices {
        assert_eq!(brep.half_edges[vertex.outgoing].origin, vh);
    }
    for &edge in &brep.boundary_edges {
        let fin = brep.edges[edge].half_edge;
        assert!(brep.half_edges.contains_key(fin), "free edge keeps its fin");
        let curve = brep.edges[edge].curve.as_ref().expect("edge curve");
        let start = brep.vertices[brep.half_edges[fin].origin].position;
        assert!(
            curve.geometry.point_at(curve.t_start).distance_to(&start) < 1e-9,
            "the free fin starts at the curve's start"
        );
    }
    brep.reverse_body_orientation();
    for (h, he) in &brep.half_edges {
        let o = &original.half_edges[h];
        assert_eq!((he.origin, he.next, he.prev), (o.origin, o.next, o.prev));
    }
    for (e, edge) in &brep.edges {
        assert_eq!(edge.half_edge, original.edges[e].half_edge);
    }
}
