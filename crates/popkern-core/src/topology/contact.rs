//! Recorded zero-thickness contacts of a body (PK-964).
//!
//! A boolean whose operands touch along an edge or at a point makes a
//! result that is non-manifold there. Under the opt-in split policy it
//! is served as manifold pieces that meet at the contact — two copies
//! of the contact edge or vertex, or two shells touching — and every
//! such contact is recorded here. The record is what makes the split
//! honest: `validate_all` accepts faces meeting only at a recorded
//! contact and still rejects every unrecorded touch, and consumers
//! (mesh welding, export) can tell a deliberate contact from a defect.
//!
//! Records name entities, so they hold only while those entities live:
//! an operator that rebuilds the body drops them, and the touch then
//! reads as the defect it is unless the operator records it again.

use serde::{Deserialize, Serialize};

use crate::geometry::Point3;
use crate::topology::{EdgeHandle, FaceHandle, VertexHandle};

/// What touches at a contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TopologyContactKind {
    /// Two pieces meet along a curve (two coincident edges).
    Edge,
    /// Two pieces meet at a point (coincident vertices, or a vertex,
    /// edge or face of one piece on another).
    Vertex,
}

/// An entity taking part in a contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ContactEntity {
    /// A vertex at the contact.
    Vertex(VertexHandle),
    /// An edge along (or through) the contact.
    Edge(EdgeHandle),
    /// A face the contact lies on.
    Face(FaceHandle),
}

/// One recorded zero-thickness contact.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TopologyContact {
    /// Edge or point contact.
    pub kind: TopologyContactKind,
    /// A point of the contact (an edge contact's midpoint).
    pub location: Point3,
    /// The entities that meet there, from every piece involved.
    pub entities: Vec<ContactEntity>,
}

impl crate::topology::BRep {
    /// Faces incident to `entity`: a vertex's faces around it, an
    /// edge's two faces, a face itself.
    pub fn contact_entity_faces(&self, entity: ContactEntity) -> Vec<FaceHandle> {
        let face_of = |he: crate::topology::HalfEdgeHandle| {
            let l = self.half_edges.get(he)?.loop_;
            Some(self.loops.get(l)?.face)
        };
        let mut faces: Vec<FaceHandle> = match entity {
            ContactEntity::Face(f) => vec![f],
            ContactEntity::Edge(e) => self
                .edges
                .get(e)
                .map(|edge| {
                    let he = edge.half_edge;
                    let twin = self.half_edges.get(he).map(|h| h.twin);
                    [Some(he), twin]
                        .into_iter()
                        .flatten()
                        .filter_map(face_of)
                        .collect()
                })
                .unwrap_or_default(),
            ContactEntity::Vertex(v) => self
                .half_edges
                .iter()
                .filter(|(_, he)| he.origin == v)
                .filter_map(|(h, _)| face_of(h))
                .collect(),
        };
        faces.sort();
        faces.dedup();
        faces
    }

    /// Face pairs that meet only at recorded contacts: every pair of
    /// distinct faces incident to one contact's entities. The
    /// self-intersection check treats them as adjacent.
    pub fn contact_face_pairs(&self) -> std::collections::BTreeSet<(FaceHandle, FaceHandle)> {
        let mut pairs = std::collections::BTreeSet::new();
        for contact in &self.topology_contacts {
            let mut faces: Vec<FaceHandle> = contact
                .entities
                .iter()
                .flat_map(|&e| self.contact_entity_faces(e))
                .collect();
            faces.sort();
            faces.dedup();
            for (i, &a) in faces.iter().enumerate() {
                for &b in &faces[i + 1..] {
                    pairs.insert((a, b));
                }
            }
        }
        pairs
    }

    /// [`Self::contact_face_pairs`] as index pairs `(i < j)` into
    /// `face_keys`: faces meeting only at a recorded contact are
    /// adjacent there by design, while an unrecorded touch is not.
    pub fn contact_index_pairs(&self, face_keys: &[FaceHandle]) -> Vec<(usize, usize)> {
        if self.topology_contacts.is_empty() {
            return Vec::new();
        }
        let index = |f: FaceHandle| face_keys.iter().position(|&k| k == f);
        self.contact_face_pairs()
            .into_iter()
            .filter_map(|(fa, fb)| {
                let (i, j) = (index(fa)?, index(fb)?);
                Some((i.min(j), i.max(j)))
            })
            .collect()
    }
}
