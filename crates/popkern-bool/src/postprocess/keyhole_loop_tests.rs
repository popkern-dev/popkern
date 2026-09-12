//! Tests for the pinched-keyhole loop merge (PK-468).

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, Loop};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// The cylinder's lateral face in its canonical keyhole form:
/// `bot_arc → seam_up → top_arc → seam_down`. Returns
/// (face, bot_arc, seam_up, top_arc, seam_down).
fn lateral_keyhole(
    brep: &BRep,
) -> (
    popkern_core::topology::FaceHandle,
    popkern_core::topology::HalfEdgeHandle,
    popkern_core::topology::HalfEdgeHandle,
    popkern_core::topology::HalfEdgeHandle,
    popkern_core::topology::HalfEdgeHandle,
) {
    for (fh, face) in brep.faces.iter() {
        let start = brep.loops[face.outer_loop].half_edge;
        let mut members = vec![start];
        let mut cursor = brep.half_edges[start].next;
        while cursor != start {
            members.push(cursor);
            cursor = brep.half_edges[cursor].next;
        }
        if members.len() != 4 {
            continue;
        }
        // The bottom arc: its twin lives on another face, and its
        // successor (a seam half) shares its origin vertex.
        for &h in &members {
            let hd = &brep.half_edges[h];
            let twin_face = brep.loops[brep.half_edges[hd.twin].loop_].face;
            let nd = &brep.half_edges[hd.next];
            if twin_face != fh && nd.origin == hd.origin {
                let seam_up = hd.next;
                let top_arc = brep.half_edges[seam_up].next;
                let seam_down = brep.half_edges[top_arc].next;
                return (fh, h, seam_up, top_arc, seam_down);
            }
        }
    }
    panic!("cylinder lateral keyhole not found");
}

/// Inverse surgery: split the keyhole into the pinch form the walk
/// emits — outer loop = bottom arc alone, inner loop = seam-up →
/// top-arc → seam-down, sharing the seam-bottom vertex with the outer.
fn split_into_pinch_form(brep: &mut BRep) {
    let (face, bot_arc, seam_up, _top_arc, seam_down) = lateral_keyhole(brep);
    let outer = brep.faces[face].outer_loop;
    brep.half_edges[bot_arc].next = bot_arc;
    brep.half_edges[bot_arc].prev = bot_arc;
    brep.half_edges[seam_down].next = seam_up;
    brep.half_edges[seam_up].prev = seam_down;
    brep.loops[outer].half_edge = bot_arc;
    let inner = brep.loops.insert(Loop {
        half_edge: seam_up,
        face,
    });
    let mut cursor = seam_up;
    loop {
        brep.half_edges[cursor].loop_ = inner;
        cursor = brep.half_edges[cursor].next;
        if cursor == seam_up {
            break;
        }
    }
    brep.faces[face].inner_loops.push(inner);
}

#[test]
fn pinched_keyhole_merges_back_to_canonical_form() {
    let mut brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    split_into_pinch_form(&mut brep);
    let (face, bot_arc, seam_up, top_arc, seam_down) = lateral_keyhole_pinched(&brep);

    let report = merge_pinched_keyhole_loops(&mut brep);
    assert_eq!(report.merged, 1, "one pinched inner loop merged");
    assert_eq!(report.declined, 0);
    assert!(brep.faces[face].inner_loops.is_empty(), "inner loop gone");
    // Canonical wiring restored: bot_arc → seam_up → top_arc → seam_down.
    assert_eq!(brep.half_edges[bot_arc].next, seam_up);
    assert_eq!(brep.half_edges[seam_up].next, top_arc);
    assert_eq!(brep.half_edges[top_arc].next, seam_down);
    assert_eq!(brep.half_edges[seam_down].next, bot_arc);
    let outer = brep.faces[face].outer_loop;
    for h in [bot_arc, seam_up, top_arc, seam_down] {
        assert_eq!(brep.half_edges[h].loop_, outer, "loop membership updated");
    }
    assert!(brep.validate().is_ok(), "merged body validates");
}

/// Locate the pinch-form pieces (outer = single bottom arc; inner =
/// 3-half-edge seam ring).
fn lateral_keyhole_pinched(
    brep: &BRep,
) -> (
    popkern_core::topology::FaceHandle,
    popkern_core::topology::HalfEdgeHandle,
    popkern_core::topology::HalfEdgeHandle,
    popkern_core::topology::HalfEdgeHandle,
    popkern_core::topology::HalfEdgeHandle,
) {
    for (fh, face) in brep.faces.iter() {
        if face.inner_loops.len() != 1 {
            continue;
        }
        let bot_arc = brep.loops[face.outer_loop].half_edge;
        let seam_up = brep.loops[face.inner_loops[0]].half_edge;
        let top_arc = brep.half_edges[seam_up].next;
        let seam_down = brep.half_edges[top_arc].next;
        return (fh, bot_arc, seam_up, top_arc, seam_down);
    }
    panic!("pinch form not found");
}

#[test]
fn canonical_bodies_are_untouched() {
    for mut brep in [
        BRep::cylinder(1.0, 2.0, tol()).unwrap(),
        BRep::unit_cube(tol()).unwrap(),
    ] {
        let faces_before = brep.faces.len();
        let loops_before = brep.loops.len();
        let report = merge_pinched_keyhole_loops(&mut brep);
        assert_eq!(report.merged, 0, "no pinched loops in canonical bodies");
        assert_eq!(report.declined, 0);
        assert_eq!(brep.faces.len(), faces_before);
        assert_eq!(brep.loops.len(), loops_before);
    }
}
