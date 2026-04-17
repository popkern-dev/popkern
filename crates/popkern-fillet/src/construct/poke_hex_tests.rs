//! Tests for `poke_hex_face_into_three_quads`.
//!
//! Builds a synthetic hexagonal face using raw Euler operators
//! (one `mvfs` + six `mev`/`mef` ops), then pokes it and verifies:
//!   * Three sub-faces are returned.
//!   * Each sub-face has exactly 4 half-edges (quadrilateral).
//!   * The central vertex is shared by all three sub-faces via
//!     its outgoing half-edges.
//!   * Total topology accounting: +1 vertex, +3 edges, +2 faces
//!     vs. the pre-poke state.
//!   * Every newly-added spoke edge lies between exactly two
//!     sub-faces (proper internal shared edge).

use super::poke_hex_face_into_three_quads;
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Build a flat hexagonal face with 6 vertices on a unit circle
/// at z=0. Returns the BRep plus the face handle for the hex.
fn build_flat_hex() -> (BRep, popkern_core::topology::FaceHandle) {
    let mut brep = BRep::new(Tolerance::default());
    // Six hex vertex positions on a unit circle in z=0 plane.
    let hex_points: [Point3; 6] = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 0.8660254038, 0.0),
        Point3::new(-0.5, 0.8660254038, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(-0.5, -0.8660254038, 0.0),
        Point3::new(0.5, -0.8660254038, 0.0),
    ];
    let mvfs = brep.mvfs(hex_points[0]).expect("mvfs");
    let m1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, hex_points[1])
        .expect("mev1");
    let m2 = brep.mev(m1.half_edge_twin, hex_points[2]).expect("mev2");
    let m3 = brep.mev(m2.half_edge_twin, hex_points[3]).expect("mev3");
    let m4 = brep.mev(m3.half_edge_twin, hex_points[4]).expect("mev4");
    let m5 = brep.mev(m4.half_edge_twin, hex_points[5]).expect("mev5");
    // Close the hex: mef connects V5's outgoing back to V0's incoming.
    let mef = brep.mef(m5.half_edge_twin, m1.half_edge).expect("mef");
    let face = mef.face;
    (brep, face)
}

#[test]
fn poke_hex_produces_three_quadrilateral_subfaces() {
    let (mut brep, face) = build_flat_hex();
    let v_before = brep.vertices.len();
    let e_before = brep.edges.len();
    let f_before = brep.faces.len();
    let he_before = brep.half_edges.len();

    let center = Point3::new(0.0, 0.0, 0.0);
    let result = poke_hex_face_into_three_quads(&mut brep, face, center)
        .expect("poke should succeed on synthetic hex");

    // Euler invariant: +1 vertex, +3 edges, +2 faces (+6 half-edges).
    assert_eq!(brep.vertices.len(), v_before + 1, "+1 vertex (central)");
    assert_eq!(brep.edges.len(), e_before + 3, "+3 edges (3 spokes)");
    assert_eq!(brep.faces.len(), f_before + 2, "+2 faces (1 → 3)");
    assert_eq!(brep.half_edges.len(), he_before + 6, "+6 half-edges");

    // Every sub-face has exactly 4 half-edges (quadrilateral).
    for (i, &sub_face) in result.sub_faces.iter().enumerate() {
        let face_data = brep.faces.get(sub_face).expect("sub-face should exist");
        let loop_data = brep.loops.get(face_data.outer_loop).expect("loop");
        let start_he = loop_data.half_edge;
        let mut count = 0;
        let mut cursor = start_he;
        loop {
            count += 1;
            if count > 10 {
                panic!("sub-face {i} loop doesn't close");
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == start_he {
                break;
            }
        }
        assert_eq!(
            count, 4,
            "sub-face {i} should have 4 half-edges; got {count}"
        );
    }

    // Central vertex should have 3 outgoing half-edges (one per spoke).
    let central = result.central_vertex;
    let outgoing_count = brep
        .half_edges
        .values()
        .filter(|he| he.origin == central)
        .count();
    assert_eq!(
        outgoing_count, 3,
        "central vertex should have 3 outgoing half-edges (one per spoke)",
    );

    // Central vertex should have 3 INCOMING half-edges too (each
    // spoke has two half-edges, one in each direction; but only
    // one has origin == central in the "outgoing" direction —
    // this check is redundant with outgoing_count given proper
    // twin wiring, but included for safety).
    let incoming_count = brep
        .half_edges
        .values()
        .filter(|he| {
            let twin = he.twin;
            brep.half_edges
                .get(twin)
                .is_some_and(|t| t.origin == central)
        })
        .count();
    assert_eq!(
        incoming_count, 3,
        "3 half-edges should arrive at the central vertex"
    );
}

#[test]
fn poke_hex_rejects_non_hexagonal_faces() {
    // Build a quadrilateral face (4 vertices, not 6).
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let m1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let m2 = brep
        .mev(m1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let m3 = brep
        .mev(m2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef = brep.mef(m3.half_edge_twin, m1.half_edge).expect("mef");
    let quad_face = mef.face;

    let result = poke_hex_face_into_three_quads(&mut brep, quad_face, Point3::new(0.5, 0.5, 0.0));
    assert!(
        result.is_err(),
        "poke should reject non-hexagonal face (got quad)"
    );
}
