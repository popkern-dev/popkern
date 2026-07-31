use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle};

fn open_box_ring(brep: &BRep) -> Vec<EdgeHandle> {
    brep.boundary_edges.iter().copied().collect()
}

fn snapshot(brep: &BRep) -> String {
    serde_json::to_string(brep).expect("BRep serializes")
}

#[test]
fn open_box_cap_closes_and_validates_strictly() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let ring = open_box_ring(&brep);
    let result = brep.cap_boundary_ring(&ring).expect("cap commits");

    assert_eq!(brep.vertices.len(), 8);
    assert_eq!(brep.edges.len(), 12);
    assert_eq!(brep.half_edges.len(), 24);
    assert_eq!(brep.faces.len(), 6);
    assert_eq!(brep.loops.len(), 6);
    assert!(brep.boundary_edges.is_empty());
    assert_eq!(result.half_edges.len(), 4);

    // The strict genus-0 formula re-arms once boundary_edges empties.
    brep.validate_euler_poincare().expect("strict Euler holds");
    brep.validate().expect("full fast validation holds");
}

#[test]
fn cap_loop_is_cycle_ordered_and_twinned() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let ring = open_box_ring(&brep);
    let result = brep.cap_boundary_ring(&ring).expect("cap commits");

    for (index, &hole_side) in result.half_edges.iter().enumerate() {
        let data = &brep.half_edges[hole_side];
        assert_eq!(
            data.next,
            result.half_edges[(index + 1) % result.half_edges.len()],
            "cap loop must be committed in cycle order"
        );
        assert_eq!(data.loop_, result.loop_);
        let surviving = &brep.half_edges[data.twin];
        assert_eq!(surviving.twin, hole_side, "twin pairing must be mutual");
        assert_eq!(
            data.origin, brep.half_edges[surviving.next].origin,
            "hole side must start at the surviving side's end vertex"
        );
    }
    assert_eq!(brep.faces[result.face].outer_loop, result.loop_);
    assert!(brep.faces[result.face].surface.is_none());
}

#[test]
fn cap_accepts_any_input_order() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let mut ring = open_box_ring(&brep);
    ring.reverse();
    brep.cap_boundary_ring(&ring).expect("cap commits");
    assert!(brep.boundary_edges.is_empty());
    brep.validate_euler_poincare().expect("strict Euler holds");
}

#[test]
fn cap_declines_short_ring_untouched() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let ring: Vec<_> = open_box_ring(&brep).into_iter().take(2).collect();
    let before = snapshot(&brep);
    let error = brep.cap_boundary_ring(&ring).expect_err("2 edges decline");
    assert!(error.to_string().contains("at least 3 edges"), "{error}");
    assert_eq!(snapshot(&brep), before, "decline must leave BRep untouched");
}

#[test]
fn cap_declines_duplicate_edge_untouched() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let mut ring = open_box_ring(&brep);
    ring.push(ring[0]);
    let before = snapshot(&brep);
    let error = brep
        .cap_boundary_ring(&ring)
        .expect_err("duplicate declines");
    assert!(error.to_string().contains("listed twice"), "{error}");
    assert_eq!(snapshot(&brep), before);
}

#[test]
fn cap_declines_non_boundary_edge_untouched() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let mut ring = open_box_ring(&brep);
    let interior = brep
        .edges
        .keys()
        .find(|edge| !brep.boundary_edges.contains(edge))
        .expect("open box has interior edges");
    ring[0] = interior;
    let before = snapshot(&brep);
    let error = brep
        .cap_boundary_ring(&ring)
        .expect_err("interior edge declines");
    assert!(
        error.to_string().contains("not a registered boundary edge"),
        "{error}"
    );
    assert_eq!(snapshot(&brep), before);
}

#[test]
fn cap_declines_unclosed_ring_untouched() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let ring: Vec<_> = open_box_ring(&brep).into_iter().take(3).collect();
    let before = snapshot(&brep);
    let error = brep
        .cap_boundary_ring(&ring)
        .expect_err("3 of 4 edges decline");
    assert!(error.to_string().contains("not closed"), "{error}");
    assert_eq!(snapshot(&brep), before);
}

#[test]
fn capped_body_is_watertight_by_manifold_checks() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let ring = open_box_ring(&brep);
    brep.cap_boundary_ring(&ring).expect("cap commits");
    brep.validate_half_edge_pairing().expect("pairing holds");
    brep.validate_manifold_edges().expect("manifold holds");
    brep.validate_vertex_valence().expect("valence holds");
}
