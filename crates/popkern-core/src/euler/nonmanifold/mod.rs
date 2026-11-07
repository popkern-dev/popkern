//! Non-manifold Euler operators (PK-1105, design §6): acorns, wire
//! edges, faces spanning wire cycles, faces splitting regions, and the
//! glue operators that merge coincident vertices, edges and faces into
//! one general body.
//!
//! An error leaves the body as it was, and each operator has an
//! inverse: the make, kill and edge/vertex glue operators check every
//! precondition before they change anything (sewing glues many edges,
//! so no whole-body copy), while `unglue_edge`, `glue_faces` and
//! `unglue_face` rewire a working copy, because whether a shell falls
//! into pieces is read after the rewiring. The
//! operators work on a body of any kind and leave it valid when it was
//! (a wire standing on a face, a lamina in a solid, two solids glued
//! along an edge are general bodies, not errors). Counts move by the
//! delta each documents over [`NonManifoldCounts`], which the tests
//! assert alongside each inverse's round trip.

mod acorn;
mod cells;
mod coincident;
mod contacts;
mod contacts_face;
mod error;
mod face_region;
mod glue_edge;
mod glue_face;
mod glue_vertex;
mod imprint_edge;
mod shells;
mod touch;
mod touch_detect;
mod vertex_loop;
mod wire_edge;
mod wire_face;

pub use acorn::AcornResult;
pub use coincident::CoincidentGlue;
pub use contacts::{ContactGlue, ContactGlueError};
pub use error::EulerError;
pub use glue_edge::GlueEdgeResult;
pub use glue_face::GlueFaceResult;
pub use glue_vertex::GlueVertexResult;
pub use imprint_edge::ImprintedEdge;
pub use touch::TouchGlue;
pub use vertex_loop::VertexLoopResult;
pub use wire_edge::{WireEdgeKill, WireEdgeResult};
pub use wire_face::WireFaceResult;

use crate::topology::BRep;

/// Every count a non-manifold operator moves.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct NonManifoldCounts {
    /// Vertices.
    pub vertices: usize,
    /// Edges.
    pub edges: usize,
    /// Fins (half-edges), face uses and wire uses alike.
    pub fins: usize,
    /// Loops.
    pub loops: usize,
    /// Faces.
    pub faces: usize,
    /// Shells of every kind (face sides, wires, acorns).
    pub shells: usize,
    /// Regions, the outside included.
    pub regions: usize,
    /// Wire edges held by shells.
    pub wire_edges: usize,
    /// Acorn shells.
    pub acorns: usize,
}

impl NonManifoldCounts {
    /// The counts of `brep`.
    pub fn of(brep: &BRep) -> Self {
        Self {
            vertices: brep.vertices.len(),
            edges: brep.edges.len(),
            fins: brep.half_edges.len(),
            loops: brep.loops.len(),
            faces: brep.faces.len(),
            shells: brep.shells.len(),
            regions: brep.regions.len(),
            wire_edges: brep.shells.values().map(|s| s.wire_edges.len()).sum(),
            acorns: brep.shells.values().filter(|s| s.acorn.is_some()).count(),
        }
    }

    /// `after - self`, field by field.
    pub fn delta(self, after: Self) -> [isize; 9] {
        let d = |a: usize, b: usize| b as isize - a as isize;
        [
            d(self.vertices, after.vertices),
            d(self.edges, after.edges),
            d(self.fins, after.fins),
            d(self.loops, after.loops),
            d(self.faces, after.faces),
            d(self.shells, after.shells),
            d(self.regions, after.regions),
            d(self.wire_edges, after.wire_edges),
            d(self.acorns, after.acorns),
        ]
    }
}

#[cfg(test)]
#[path = "nonmanifold_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "wire_face_tests.rs"]
mod wire_face_tests;

#[cfg(test)]
#[path = "glue_vertex_tests.rs"]
mod glue_vertex_tests;

#[cfg(test)]
#[path = "glue_edge_tests.rs"]
mod glue_edge_tests;

#[cfg(test)]
#[path = "glue_face_tests.rs"]
mod glue_face_tests;

#[cfg(test)]
#[path = "face_region_tests.rs"]
mod face_region_tests;

#[cfg(test)]
#[path = "contacts_tests.rs"]
mod contacts_tests;

#[cfg(test)]
#[path = "cells_tests.rs"]
mod cells_tests;
