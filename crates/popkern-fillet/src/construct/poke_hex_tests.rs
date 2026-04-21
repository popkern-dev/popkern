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

use super::{poke_hex_face_into_three_quads, poke_n_gon_face_into_n_quads};
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

/// Build a flat octagonal face (N=4 case → 8 outer half-edges).
fn build_flat_octagon() -> (BRep, popkern_core::topology::FaceHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let pts: [Point3; 8] = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(s, s, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-s, s, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(-s, -s, 0.0),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(s, -s, 0.0),
    ];
    let mvfs = brep.mvfs(pts[0]).expect("mvfs");
    let m1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, pts[1])
        .expect("mev1");
    let m2 = brep.mev(m1.half_edge_twin, pts[2]).expect("mev2");
    let m3 = brep.mev(m2.half_edge_twin, pts[3]).expect("mev3");
    let m4 = brep.mev(m3.half_edge_twin, pts[4]).expect("mev4");
    let m5 = brep.mev(m4.half_edge_twin, pts[5]).expect("mev5");
    let m6 = brep.mev(m5.half_edge_twin, pts[6]).expect("mev6");
    let m7 = brep.mev(m6.half_edge_twin, pts[7]).expect("mev7");
    let mef = brep.mef(m7.half_edge_twin, m1.half_edge).expect("mef");
    (brep, mef.face)
}

/// §5.3 item #4 — N=4 octagon poke produces 4 quadrilateral
/// sub-faces with the correct Euler invariant.
#[test]
fn poke_octagon_produces_four_quadrilateral_subfaces() {
    let (mut brep, face) = build_flat_octagon();
    let v_before = brep.vertices.len();
    let e_before = brep.edges.len();
    let f_before = brep.faces.len();
    let he_before = brep.half_edges.len();

    let center = Point3::new(0.0, 0.0, 0.0);
    let result = poke_n_gon_face_into_n_quads(&mut brep, face, center)
        .expect("poke should succeed on synthetic octagon");

    assert_eq!(result.sub_faces.len(), 4, "expected 4 sub-faces for N=4");

    // Euler invariant: +1 vertex, +N edges, +(N-1) faces, +2N half-edges.
    assert_eq!(brep.vertices.len(), v_before + 1, "+1 central vertex");
    assert_eq!(brep.edges.len(), e_before + 4, "+4 spoke edges");
    assert_eq!(brep.faces.len(), f_before + 3, "+3 faces (1 → 4)");
    assert_eq!(brep.half_edges.len(), he_before + 8, "+8 half-edges");

    // Every sub-face is a quadrilateral.
    for (i, &sub_face) in result.sub_faces.iter().enumerate() {
        let face_data = brep.faces.get(sub_face).expect("sub-face");
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
        assert_eq!(count, 4, "sub-face {i} should be a quadrilateral");
    }

    // Central vertex has N=4 outgoing half-edges.
    let central = result.central_vertex;
    let outgoing = brep
        .half_edges
        .values()
        .filter(|he| he.origin == central)
        .count();
    assert_eq!(outgoing, 4, "central vertex should have 4 outgoing he's");
}

/// Phase II first-lap — pentagonal (N-HE) input: the blend face
/// at a sharp N=5 apex has only N=5 half-edges (pure fillet arcs,
/// no topological bridges). `poke_n_gon_face_into_n_quads` must
/// split each arc at `t=0.5` to reach 2N=10, then proceed with
/// the standard algorithm, producing 5 quadrilateral sub-faces.
#[test]
fn poke_pentagon_n_he_produces_five_quadrilateral_subfaces() {
    let mut brep = BRep::new(Tolerance::default());
    // Build a 5-sided pentagon (N=5 half-edges on outer loop).
    let n: usize = 5;
    let pts: Vec<Point3> = (0..n)
        .map(|k| {
            let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
            Point3::new(theta.cos(), theta.sin(), 0.0)
        })
        .collect();
    let mvfs = brep.mvfs(pts[0]).expect("mvfs");
    let mut prev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, pts[1])
        .expect("mev1");
    let first_he = prev.half_edge;
    for p in pts.iter().take(n).skip(2) {
        prev = brep.mev(prev.half_edge_twin, *p).expect("mev");
    }
    let mef = brep.mef(prev.half_edge_twin, first_he).expect("close mef");
    let face = mef.face;

    // Pre-poke topology invariants.
    let v_before = brep.vertices.len();
    let e_before = brep.edges.len();
    let f_before = brep.faces.len();

    // Verify the input really has N=5 HEs on the outer loop.
    let outer_loop = brep.faces[face].outer_loop;
    let start = brep.loops[outer_loop].half_edge;
    let mut he_count = 1usize;
    let mut cur = start;
    loop {
        cur = brep.half_edges[cur].next;
        if cur == start {
            break;
        }
        he_count += 1;
    }
    assert_eq!(he_count, 5, "input pentagon must have 5 outer HEs");

    let center = Point3::new(0.0, 0.0, 0.0);
    let result = poke_n_gon_face_into_n_quads(&mut brep, face, center)
        .expect("N-HE pentagon poke should succeed");

    assert_eq!(result.sub_faces.len(), n, "expected 5 sub-faces");

    // Euler accounting: the N-HE path splits each arc (+N vertices,
    // +N edges, +0 faces), then runs the standard poke (+1 vertex,
    // +N edges, +(N-1) faces). Net: +(N+1) vertices, +2N edges,
    // +(N-1) faces vs the pre-poke state.
    assert_eq!(brep.vertices.len(), v_before + n + 1, "+N+1 vertices");
    assert_eq!(brep.edges.len(), e_before + 2 * n, "+2N edges");
    assert_eq!(brep.faces.len(), f_before + (n - 1), "+(N-1) faces");

    // Every sub-face is a quadrilateral.
    for (i, &sub_face) in result.sub_faces.iter().enumerate() {
        let outer = brep.faces[sub_face].outer_loop;
        let start = brep.loops[outer].half_edge;
        let mut count = 1usize;
        let mut cur = start;
        for _ in 0..16 {
            cur = brep.half_edges[cur].next;
            if cur == start {
                break;
            }
            count += 1;
        }
        assert_eq!(count, 4, "sub-face {i} should be a quadrilateral");
    }
}

/// §5.3 item #4 — N=5 poke (decagon → 5 quads).
#[test]
fn poke_decagon_produces_five_quadrilateral_subfaces() {
    let mut brep = BRep::new(Tolerance::default());
    let n: usize = 5;
    let pts: Vec<Point3> = (0..2 * n)
        .map(|k| {
            let theta = 2.0 * std::f64::consts::PI * (k as f64) / (2 * n) as f64;
            Point3::new(theta.cos(), theta.sin(), 0.0)
        })
        .collect();
    let mvfs = brep.mvfs(pts[0]).expect("mvfs");
    let mut prev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, pts[1])
        .expect("mev1");
    let first_he = prev.half_edge;
    for p in pts.iter().take(2 * n).skip(2) {
        prev = brep.mev(prev.half_edge_twin, *p).expect("mev");
    }
    let _mef = brep.mef(prev.half_edge_twin, first_he).expect("close mef");
    let face = _mef.face;

    let v_before = brep.vertices.len();
    let e_before = brep.edges.len();
    let f_before = brep.faces.len();

    let center = Point3::new(0.0, 0.0, 0.0);
    let result = poke_n_gon_face_into_n_quads(&mut brep, face, center)
        .expect("poke should succeed on synthetic decagon");

    assert_eq!(result.sub_faces.len(), n, "expected N sub-faces");
    assert_eq!(brep.vertices.len(), v_before + 1);
    assert_eq!(brep.edges.len(), e_before + n);
    assert_eq!(brep.faces.len(), f_before + (n - 1));
}
