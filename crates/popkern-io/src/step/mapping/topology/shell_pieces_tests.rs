//! The split of a closed shell into touching pieces, on face topology
//! alone (no geometry).

use std::collections::{BTreeMap, HashMap};

use super::touching_pieces;
use crate::step::mapping::topology::face_collector::{
    BoundInfo, FaceData, FaceInfo, OrientedEdgeInfo,
};
use crate::step::parser::types::EntityId;

/// A shell of faces, each a loop of the given edges; `ends` gives each
/// edge's vertices.
fn shell(faces: &[&[EntityId]], ends: &[(EntityId, (EntityId, EntityId))]) -> FaceData {
    let edge_ends: HashMap<_, _> = ends.iter().copied().collect();
    let mut edge_curve_ids: Vec<EntityId> = ends.iter().map(|(e, _)| *e).collect();
    edge_curve_ids.sort_unstable();
    let mut vertex_ids: Vec<EntityId> = ends.iter().flat_map(|(_, (a, b))| [*a, *b]).collect();
    vertex_ids.sort_unstable();
    vertex_ids.dedup();
    FaceData {
        vertex_ids,
        edge_curve_ids,
        edge_ends,
        faces: faces
            .iter()
            .enumerate()
            .map(|(k, edges)| FaceInfo {
                face_id: 1000 + k as EntityId,
                surface_id: 2000 + k as EntityId,
                same_sense: true,
                bounds: vec![BoundInfo {
                    is_outer: true,
                    oriented_edges: edges
                        .iter()
                        .map(|&e| OrientedEdgeInfo {
                            oriented_edge_id: 3000 + e,
                            edge_curve_id: e,
                            orientation: true,
                        })
                        .collect(),
                }],
                vertex_loops: Vec::new(),
            })
            .collect(),
    }
}

fn uses(data: &FaceData) -> BTreeMap<EntityId, usize> {
    let mut out = BTreeMap::new();
    for face in &data.faces {
        for oe in face.bounds.iter().flat_map(|b| &b.oriented_edges) {
            *out.entry(oe.edge_curve_id).or_default() += 1;
        }
    }
    out
}

/// Two "lenses" (two faces on two edges each) sharing edge 9.
#[test]
fn two_pieces_on_one_edge_split_with_each_piece_using_it_twice() {
    let ends = [(1, (1, 2)), (2, (1, 2)), (9, (1, 2))];
    let data = shell(&[&[9, 1], &[1, 9], &[9, 2], &[2, 9]], &ends);
    let pieces = touching_pieces(&data, &uses(&data))
        .ok()
        .flatten()
        .expect("split");
    assert_eq!(pieces.count(), 2);
    assert_eq!(pieces.of_face, vec![0, 0, 1, 1]);
    assert!(pieces.edges.iter().all(|e| e.contains(&9)));
}

/// The same four uses of edge 9, but the two lenses joined through
/// edge 3 into one piece: only the faces' radial order could pair them.
#[test]
fn a_piece_using_an_edge_four_times_declines() {
    let ends = [(1, (1, 2)), (2, (1, 2)), (3, (2, 3)), (9, (1, 2))];
    let data = shell(&[&[9, 1], &[1, 9, 3], &[9, 2, 3], &[2, 9]], &ends);
    let Err(unpaired) = touching_pieces(&data, &uses(&data)) else {
        panic!("a self-touching piece must decline");
    };
    assert_eq!((unpaired.edge, unpaired.uses_in_piece), (9, 4));
}

#[test]
fn pieces_meeting_at_a_vertex_split() {
    let ends = [(1, (1, 2)), (2, (2, 1)), (3, (2, 3)), (4, (3, 2))];
    let data = shell(&[&[1, 2], &[2, 1], &[3, 4], &[4, 3]], &ends);
    let pieces = touching_pieces(&data, &uses(&data))
        .ok()
        .flatten()
        .expect("split");
    assert_eq!(pieces.of_face, vec![0, 0, 1, 1]);
}

/// Pieces sharing nothing (a shell listing two apart bodies) are left
/// as one body, as before.
#[test]
fn apart_pieces_stay_whole() {
    let ends = [(1, (1, 2)), (2, (2, 1)), (3, (3, 4)), (4, (4, 3))];
    let data = shell(&[&[1, 2], &[2, 1], &[3, 4], &[4, 3]], &ends);
    assert!(matches!(touching_pieces(&data, &uses(&data)), Ok(None)));
}
