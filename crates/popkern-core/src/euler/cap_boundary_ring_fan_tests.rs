use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle};

fn open_box_ring(brep: &BRep) -> Vec<EdgeHandle> {
    brep.boundary_edges.iter().copied().collect()
}

fn fan_pole() -> Point3 {
    // open_box removes one cube face; the centroid of the opening
    // nudged inward is a valid pole for every removal choice.
    Point3::new(0.5, 0.5, 0.5)
}

#[test]
fn open_box_fan_closes_and_validates_strictly() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let ring = open_box_ring(&brep);
    let result = brep
        .cap_boundary_ring_fan(&ring, fan_pole())
        .expect("fan commits");

    assert_eq!(brep.vertices.len(), 9);
    assert_eq!(brep.edges.len(), 16);
    assert_eq!(brep.half_edges.len(), 32);
    assert_eq!(brep.faces.len(), 9);
    assert_eq!(brep.loops.len(), 9);
    assert!(brep.boundary_edges.is_empty());
    assert_eq!(result.sector_faces.len(), 4);
    assert_eq!(result.spoke_edges.len(), 4);
    assert_eq!(result.rim_half_edges.len(), 4);

    brep.validate_euler_poincare().expect("strict Euler holds");
    brep.validate().expect("full fast validation holds");
    brep.validate_vertex_valence().expect("valence holds");
}

#[test]
fn fan_sectors_are_triangles_with_shared_spokes() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let ring = open_box_ring(&brep);
    let original_ring: std::collections::BTreeSet<_> = ring.iter().copied().collect();
    let result = brep
        .cap_boundary_ring_fan(&ring, fan_pole())
        .expect("fan commits");
    let count = result.sector_faces.len();

    assert_eq!(
        brep.vertices[result.pole_vertex].position,
        fan_pole(),
        "pole vertex sits at the requested position"
    );

    for (index, (&face, &rim)) in result
        .sector_faces
        .iter()
        .zip(&result.rim_half_edges)
        .enumerate()
    {
        assert!(
            original_ring.contains(&brep.half_edges[rim].edge),
            "sector {index} rim half-edge must consume a ring edge"
        );
        // Collect the sector loop's edges: rim k plus spokes k and k+1.
        let loop_ = brep.faces[face].outer_loop;
        let start = brep.loops[loop_].half_edge;
        let mut cursor = start;
        let mut loop_edges = std::collections::BTreeSet::new();
        let mut cycle_len = 0;
        loop {
            loop_edges.insert(brep.half_edges[cursor].edge);
            cycle_len += 1;
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
            assert!(cycle_len <= 3, "sector {index} loop exceeds a triangle");
        }
        assert_eq!(cycle_len, 3, "sector {index} loop is a triangle");
        assert!(loop_edges.contains(&brep.half_edges[rim].edge));
        assert!(
            loop_edges.contains(&result.spoke_edges[index]),
            "sector {index} shares spoke {index}"
        );
        assert!(
            loop_edges.contains(&result.spoke_edges[(index + 1) % count]),
            "sector {index} shares spoke {}",
            (index + 1) % count
        );
    }
}

#[test]
fn fan_spokes_connect_ring_vertices_to_pole() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let ring = open_box_ring(&brep);
    let result = brep
        .cap_boundary_ring_fan(&ring, fan_pole())
        .expect("fan commits");

    for (index, &spoke) in result.spoke_edges.iter().enumerate() {
        let ring_vertex = brep.half_edges[result.rim_half_edges[index]].origin;
        let he = brep.edges[spoke].half_edge;
        let origin = brep.half_edges[he].origin;
        let end = brep.half_edges[brep.half_edges[he].twin].origin;
        let endpoints = std::collections::BTreeSet::from([origin, end]);
        assert_eq!(
            endpoints,
            std::collections::BTreeSet::from([ring_vertex, result.pole_vertex]),
            "spoke {index} must connect ring vertex {index} to the pole"
        );
    }
}

#[test]
fn fan_decline_propagates_and_leaves_brep_untouched() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let ring: Vec<_> = open_box_ring(&brep).into_iter().take(3).collect();
    let before = serde_json::to_string(&brep).expect("BRep serializes");
    let error = brep
        .cap_boundary_ring_fan(&ring, fan_pole())
        .expect_err("unclosed ring declines");
    assert!(error.to_string().contains("not closed"), "{error}");
    let after = serde_json::to_string(&brep).expect("BRep serializes");
    assert_eq!(before, after, "fan decline must leave the BRep untouched");
}
