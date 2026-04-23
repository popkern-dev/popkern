//! Tests for [`poke_n_gon_face_into_4n_quads`].
//!
//! Uses synthetic 2N-gon faces built via raw Euler operators
//! (`mvfs` + `mev` + closing `mef`) and verifies:
//! * `4N` sub-faces returned.
//! * Each sub-face is a quadrilateral (4 half-edges).
//! * Central vertex has `N` outgoing half-edges (one per sector's
//!   `C_star` piece — the radial midpoints sit on spokes *from* C).
//! * Euler invariant: `ΔV = +4N, ΔE = +7N, ΔF = +3N, ΔHE = +14N`
//!   (counting the outer N-way poke + N radial splits + 2N outer
//!   splits + 4N interior edges).

#![allow(clippy::needless_range_loop)]

use super::poke_n_gon_face_into_4n_quads;
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Build a flat 2N-gon face with vertices on a unit circle (z=0).
/// The 2N HE layout mimics the R8-C1 arc/bridge alternating output.
fn build_flat_2ngon(n: usize) -> (BRep, FaceHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let count = 2 * n;
    let pts: Vec<Point3> = (0..count)
        .map(|k| {
            let theta = 2.0 * std::f64::consts::PI * (k as f64) / (count as f64);
            Point3::new(theta.cos(), theta.sin(), 0.0)
        })
        .collect();
    let mvfs = brep.mvfs(pts[0]).expect("mvfs");
    let m1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, pts[1])
        .expect("mev1");
    let first_he = m1.half_edge;
    let mut prev = m1;
    for p in pts.iter().take(count).skip(2) {
        prev = brep.mev(prev.half_edge_twin, *p).expect("mev");
    }
    let mef = brep.mef(prev.half_edge_twin, first_he).expect("close mef");
    (brep, mef.face)
}

/// Count half-edges on a face's outer loop.
fn loop_he_count(brep: &BRep, face: FaceHandle) -> usize {
    let outer = brep.faces[face].outer_loop;
    let start = brep.loops[outer].half_edge;
    let mut count = 0;
    let mut cursor = start;
    let bound = brep.half_edges.len() + 1;
    for _ in 0..bound {
        count += 1;
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return count;
        }
    }
    panic!("loop did not close");
}

fn poke_and_check_basics(n: usize) {
    let (mut brep, face) = build_flat_2ngon(n);
    let v_before = brep.vertices.len();
    let e_before = brep.edges.len();
    let f_before = brep.faces.len();
    let he_before = brep.half_edges.len();

    let center = Point3::new(0.0, 0.0, 0.0);
    let result = poke_n_gon_face_into_4n_quads(&mut brep, face, center).expect("4n poke succeeds");

    // Count.
    assert_eq!(
        result.sub_faces.len(),
        4 * n,
        "expected 4N={} sub-faces at N={n}",
        4 * n
    );

    // Each sub-face is a quadrilateral.
    for (idx, &sub_face) in result.sub_faces.iter().enumerate() {
        let cnt = loop_he_count(&brep, sub_face);
        assert_eq!(cnt, 4, "sub-face {idx} (N={n}) expected 4 HEs, got {cnt}");
    }

    // Euler invariant from the composite operation:
    //  N-way poke:    +1 V, +N E, +(N-1) F, +2N HE
    //  N radial splits: +N V, +N E, 0 F, +2N HE
    //  2N outer splits: +2N V, +2N E, 0 F, +4N HE
    //  N × (mev + 3 mef): +N V, +4N E, +3N F, +8N HE
    //                                 (mev: +1 V, +1 E, +2 HE;
    //                                  mef: +1 E, +1 F, +2 HE)
    //
    // Totals: ΔV = 1 + 4N, ΔE = 8N, ΔF = 4N−1, ΔHE = 16N.
    assert_eq!(
        brep.vertices.len(),
        v_before + 1 + 4 * n,
        "ΔV mismatch at N={n}"
    );
    assert_eq!(brep.edges.len(), e_before + 8 * n, "ΔE mismatch at N={n}");
    assert_eq!(
        brep.faces.len(),
        f_before + 4 * n - 1,
        "ΔF mismatch at N={n}"
    );
    assert_eq!(
        brep.half_edges.len(),
        he_before + 16 * n,
        "ΔHE mismatch at N={n}"
    );

    // Central vertex has N outgoing half-edges (the radial spokes
    // from C out to radial midpoints).
    let central = result.central_vertex;
    let outgoing = brep
        .half_edges
        .values()
        .filter(|he| he.origin == central)
        .count();
    assert_eq!(outgoing, n, "central vertex should have {n} outgoing HEs");
}

#[test]
fn poke_4n_pentagon() {
    poke_and_check_basics(5);
}

#[test]
fn poke_4n_hexagon() {
    poke_and_check_basics(6);
}

#[test]
fn poke_4n_septagon_octagon_decagon() {
    for n in [7_usize, 8, 10] {
        poke_and_check_basics(n);
    }
}

#[test]
fn poke_4n_low_valence_n3_n4() {
    for n in [3_usize, 4] {
        poke_and_check_basics(n);
    }
}

#[test]
fn poke_4n_sub_faces_disjoint() {
    let n = 6;
    let (mut brep, face) = build_flat_2ngon(n);
    let center = Point3::new(0.0, 0.0, 0.0);
    let result = poke_n_gon_face_into_4n_quads(&mut brep, face, center).unwrap();

    // All 4N sub-faces have distinct handles.
    let mut sorted: Vec<_> = result.sub_faces.iter().collect();
    sorted.sort();
    sorted.dedup();
    assert_eq!(sorted.len(), 4 * n, "sub-faces must be distinct");
}
