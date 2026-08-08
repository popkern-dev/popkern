//! Tests for whole-body orientation reversal.

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
        assert_eq!(round.twin, he.twin);
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
        let twin = &brep.half_edges[he.twin];
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
                brep.half_edges[he.twin].origin, next.origin,
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
