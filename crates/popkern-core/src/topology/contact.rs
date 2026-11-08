//! A place where a body touches itself (PK-964, PK-1110, PK-1115).
//!
//! A general body holds its contacts in its topology: an edge whose
//! ring has more than two face uses, a vertex with more than one use.
//! [`crate::topology::BRep::contacts`] lists them as [`TopologyContact`]s, named
//! by the shared entity. Nothing is stored: manifold pieces that touch
//! without sharing are separate bodies.

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
