//! Tests for region removal and the boundary-weld preconditions
//! (PK-327 re-close substrate). The weld's success path is exercised
//! end-to-end by the popkern-bool self-trim gates — a valid embedded
//! body cannot carry two coincident boundary rims to weld in
//! isolation.

use std::collections::BTreeSet;

use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn removing_one_face_matches_open_box() {
    let mut cube = BRep::unit_cube(tol()).unwrap();
    let face = cube.faces.keys().next().unwrap();
    let result = cube
        .remove_face_region_to_boundary(&BTreeSet::from([face]))
        .unwrap();
    assert_eq!(cube.faces.len(), 5);
    assert_eq!(cube.boundary_edges.len(), 4);
    assert_eq!(result.rim_edges.len(), 4);
    assert_eq!(result.removed_edges, 0);
    assert_eq!(result.removed_vertices, 0);
    assert!(cube.validate().is_ok());
}

#[test]
fn removing_two_adjacent_faces_deletes_the_shared_edge() {
    let mut cube = BRep::unit_cube(tol()).unwrap();
    // Find two faces sharing an edge.
    let faces: Vec<_> = cube.faces.keys().collect();
    let mut pair = None;
    'outer: for &fa in &faces {
        for &fb in &faces {
            if fa == fb {
                continue;
            }
            let la = cube.faces[fa].outer_loop;
            let start = cube.loops[la].half_edge;
            let mut cursor = start;
            loop {
                let he = &cube.half_edges[cursor];
                let twin = &cube.half_edges[he.twin];
                if cube.loops[twin.loop_].face == fb {
                    pair = Some((fa, fb));
                    break 'outer;
                }
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
        }
    }
    let (fa, fb) = pair.expect("adjacent pair");
    let result = cube
        .remove_face_region_to_boundary(&BTreeSet::from([fa, fb]))
        .unwrap();
    assert_eq!(cube.faces.len(), 4);
    // Two adjacent quads: 1 shared edge dies, 6 rim edges open.
    assert_eq!(result.removed_edges, 1);
    assert_eq!(cube.boundary_edges.len(), 6);
    assert_eq!(cube.edges.len(), 11);
    assert!(cube.validate().is_ok());
}

#[test]
fn removing_everything_is_a_precondition_failure() {
    let mut cube = BRep::unit_cube(tol()).unwrap();
    let all: BTreeSet<_> = cube.faces.keys().collect();
    assert!(cube.remove_face_region_to_boundary(&all).is_err());
    // Model untouched.
    assert_eq!(cube.faces.len(), 6);
    assert!(cube.validate().is_ok());
}

#[test]
fn weld_rejects_non_boundary_edges() {
    let mut cube = BRep::unit_cube(tol()).unwrap();
    let edges: Vec<_> = cube.edges.keys().collect();
    let err = cube
        .weld_boundary_edge_pair(edges[0], edges[1], &tol())
        .unwrap_err();
    assert!(
        format!("{err}").contains("not a boundary edge"),
        "unexpected: {err}"
    );
}

#[test]
fn weld_rejects_non_coincident_boundary_pairs() {
    // Open two OPPOSITE faces: the two rims are far apart and must
    // not weld.
    let mut cube = BRep::unit_cube(tol()).unwrap();
    let bottom = cube
        .faces
        .keys()
        .find(|&fh| {
            let ps = cube.face_outer_loop_positions(fh);
            !ps.is_empty() && ps.iter().all(|p| p.z.abs() < 1e-9)
        })
        .unwrap();
    let top = cube
        .faces
        .keys()
        .find(|&fh| {
            let ps = cube.face_outer_loop_positions(fh);
            !ps.is_empty() && ps.iter().all(|p| (p.z - 1.0).abs() < 1e-9)
        })
        .unwrap();
    cube.remove_face_region_to_boundary(&BTreeSet::from([bottom]))
        .unwrap();
    cube.remove_face_region_to_boundary(&BTreeSet::from([top]))
        .unwrap();
    let boundary: Vec<_> = cube.boundary_edges.iter().copied().collect();
    assert_eq!(boundary.len(), 8);
    let err = cube
        .weld_boundary_edge_pair(boundary[0], boundary[7], &tol())
        .unwrap_err();
    assert!(
        format!("{err}").contains("do not oppose") || format!("{err}").contains("disagree"),
        "unexpected: {err}"
    );
}
