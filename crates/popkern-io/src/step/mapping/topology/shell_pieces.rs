//! A closed shell that touches itself: STEP's non-manifold bodies
//! (PK-980).
//!
//! A writer may list the faces of several solids that touch in one
//! `CLOSED_SHELL`: two blocks sharing an edge list it with four face
//! uses, blocks meeting at a corner share its vertex. The kernel's
//! half-edges pair exactly two uses per edge, so such a shell imports
//! as manifold pieces meeting at recorded contacts, the form a boolean
//! serves under the split policy: each piece its own solid with its own
//! copies of the contact edges and vertices, and every contact recorded
//! on the body (`BRep::topology_contacts`), which `validate_all`
//! excuses and the mesher keeps apart.
//!
//! The pieces are the face sets joined through edges of exactly two
//! uses. That reads the file's topology alone, so it is decided only
//! where each piece uses every contact edge exactly twice; a surface
//! touching itself along an edge (one piece using an edge four times)
//! would need the faces' radial order to pair its uses, and declines,
//! as does a contact edge some piece uses once. Pieces that share
//! nothing stay one body, as before.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, ContactEntity, TopologyContact, TopologyContactKind};

use super::face_collector::FaceData;
use crate::step::parser::types::EntityId;

/// The pieces of a shell: each face's piece, and each piece's vertices
/// and edges.
pub(super) struct ShellPieces {
    /// The piece of each face, in `FaceData::faces` order.
    pub of_face: Vec<usize>,
    /// The STEP vertices of each piece.
    pub vertices: Vec<BTreeSet<EntityId>>,
    /// The STEP edges of each piece.
    pub edges: Vec<BTreeSet<EntityId>>,
}

impl ShellPieces {
    /// One piece holding the whole shell.
    pub fn whole(face_data: &FaceData) -> Self {
        Self {
            of_face: vec![0; face_data.faces.len()],
            vertices: vec![face_data.vertex_ids.iter().copied().collect()],
            edges: vec![face_data.edge_curve_ids.iter().copied().collect()],
        }
    }

    /// The number of pieces.
    pub fn count(&self) -> usize {
        self.edges.len()
    }
}

/// Why a non-manifold shell has no split into closed pieces.
pub(super) struct Unpaired {
    /// The edge whose uses no split pairs.
    pub edge: EntityId,
    /// Its uses in the piece that holds it an odd or excess number of times.
    pub uses_in_piece: usize,
}

fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// The shell's pieces when it touches itself (module doc): `Ok(None)`
/// when it is one piece or its pieces share nothing.
pub(super) fn touching_pieces(
    face_data: &FaceData,
    uses: &BTreeMap<EntityId, usize>,
) -> Result<Option<ShellPieces>, Unpaired> {
    let n = face_data.faces.len();
    let face_edges: Vec<Vec<EntityId>> = face_data
        .faces
        .iter()
        .map(|f| {
            f.bounds
                .iter()
                .flat_map(|b| b.oriented_edges.iter().map(|oe| oe.edge_curve_id))
                .collect()
        })
        .collect();
    // Join faces through edges of exactly two uses.
    let mut parent: Vec<usize> = (0..n).collect();
    let mut first_use: BTreeMap<EntityId, usize> = BTreeMap::new();
    for (k, edges) in face_edges.iter().enumerate() {
        for &e in edges {
            if uses.get(&e) != Some(&2) {
                continue;
            }
            match first_use.get(&e) {
                Some(&j) => {
                    let (a, b) = (find(&mut parent, j), find(&mut parent, k));
                    parent[a] = b;
                }
                None => {
                    first_use.insert(e, k);
                }
            }
        }
    }
    let mut index: BTreeMap<usize, usize> = BTreeMap::new();
    let of_face: Vec<usize> = (0..n)
        .map(|k| {
            let root = find(&mut parent, k);
            let next = index.len();
            *index.entry(root).or_insert(next)
        })
        .collect();
    let count = index.len();
    let mut edges = vec![BTreeSet::new(); count];
    let mut vertices = vec![BTreeSet::new(); count];
    let mut piece_uses: BTreeMap<(EntityId, usize), usize> = BTreeMap::new();
    for (k, list) in face_edges.iter().enumerate() {
        let p = of_face[k];
        for &e in list {
            edges[p].insert(e);
            *piece_uses.entry((e, p)).or_default() += 1;
            if let Some(&(a, b)) = face_data.edge_ends.get(&e) {
                vertices[p].insert(a);
                vertices[p].insert(b);
            }
        }
    }
    if let Some((&(edge, _), &uses_in_piece)) = piece_uses.iter().find(|(_, u)| **u != 2) {
        return Err(Unpaired {
            edge,
            uses_in_piece,
        });
    }
    let touching = (0..count)
        .any(|p| (p + 1..count).any(|q| vertices[p].intersection(&vertices[q]).next().is_some()));
    if count < 2 || !touching {
        return Ok(None);
    }
    Ok(Some(ShellPieces {
        of_face,
        vertices,
        edges,
    }))
}

/// Each of `ids`' copies over `count` pieces, in piece order.
pub(super) fn copies<H>(
    ids: &[EntityId],
    count: usize,
    copy: impl Fn(&(EntityId, usize)) -> Option<H>,
) -> BTreeMap<EntityId, Vec<H>> {
    ids.iter()
        .map(|&id| (id, (0..count).filter_map(|p| copy(&(id, p))).collect()))
        .collect()
}

/// Record the contacts between a shell's pieces: every STEP vertex and
/// edge copied into more than one piece, each copy named, as the split
/// policy records them.
pub(super) fn record_contacts(
    brep: &mut BRep,
    vertex_copies: &BTreeMap<EntityId, Vec<popkern_core::topology::VertexHandle>>,
    edge_copies: &BTreeMap<EntityId, Vec<popkern_core::topology::EdgeHandle>>,
) {
    let mut records = Vec::new();
    for copies in vertex_copies.values().filter(|c| c.len() > 1) {
        let Some(location) = brep.vertices.get(copies[0]).map(|v| v.position) else {
            continue;
        };
        records.push(TopologyContact {
            kind: TopologyContactKind::Vertex,
            location,
            entities: copies.iter().map(|&v| ContactEntity::Vertex(v)).collect(),
        });
    }
    for copies in edge_copies.values().filter(|c| c.len() > 1) {
        let Some(edge) = brep.edges.get(copies[0]) else {
            continue;
        };
        let location = match edge.curve.as_ref() {
            Some(c) => c.geometry.point_at(0.5 * (c.t_start + c.t_end)),
            None => {
                let Some((a, b)) = brep.half_edges.get(edge.half_edge).and_then(|h| {
                    Some((
                        h.origin,
                        h.mate().and_then(|m| brep.half_edges.get(m))?.origin,
                    ))
                }) else {
                    continue;
                };
                let (a, b) = (brep.vertices[a].position, brep.vertices[b].position);
                a + (b - a) * 0.5
            }
        };
        records.push(TopologyContact {
            kind: TopologyContactKind::Edge,
            location,
            entities: copies.iter().map(|&e| ContactEntity::Edge(e)).collect(),
        });
    }
    brep.topology_contacts.extend(records);
}

#[cfg(test)]
#[path = "shell_pieces_tests.rs"]
mod tests;
