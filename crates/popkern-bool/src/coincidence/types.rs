//! Data types for the coincidence engine.
//!
//! The [`CoincidenceMap`] stores all detected coincidences between two bodies
//! and is consumed during classification and assembly to handle ON-face rules.

use popkern_core::geometry::Point3;
use popkern_core::topology::{EdgeHandle, FaceHandle, VertexHandle};

/// Orientation relationship between two coincident faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoincidenceOrientation {
    /// Face normals point in the same direction.
    Same,
    /// Face normals point in opposite directions.
    Opposite,
}

/// A detected coincidence between two faces.
#[derive(Debug, Clone)]
pub struct FaceOverlap {
    /// Face from body A.
    pub face_a: FaceHandle,
    /// Face from body B.
    pub face_b: FaceHandle,
    /// Normal orientation relationship.
    pub orientation: CoincidenceOrientation,
    /// Whether face A entirely contains face B (or vice versa).
    pub full_overlap: bool,
    /// Approximate center of the overlap region (for diagnostics).
    pub overlap_center: Point3,
}

/// A detected overlap between two edges.
#[derive(Debug, Clone)]
pub struct EdgeOverlap {
    /// Edge from body A.
    pub edge_a: EdgeHandle,
    /// Edge from body B.
    pub edge_b: EdgeHandle,
    /// Whether the edge directions are reversed relative to each other.
    pub reversed: bool,
}

/// A pair of vertices that should be merged.
#[derive(Debug, Clone)]
pub struct VertexMerge {
    /// Vertex from body A.
    pub vertex_a: VertexHandle,
    /// Vertex from body B.
    pub vertex_b: VertexHandle,
    /// Distance between the two vertices.
    pub distance: f64,
}

/// Complete coincidence map between two bodies.
///
/// Computed during the intersection phase and consumed during classification
/// and assembly. Stores all vertex merges, edge overlaps, and face overlaps
/// detected between the two input bodies.
#[derive(Debug, Clone, Default)]
pub struct CoincidenceMap {
    /// Vertex pairs to merge.
    pub vertex_merges: Vec<VertexMerge>,
    /// Overlapping edge pairs.
    pub edge_overlaps: Vec<EdgeOverlap>,
    /// Overlapping face pairs.
    pub face_overlaps: Vec<FaceOverlap>,
}

impl CoincidenceMap {
    /// Returns true if no coincidences were detected.
    pub fn is_empty(&self) -> bool {
        self.vertex_merges.is_empty()
            && self.edge_overlaps.is_empty()
            && self.face_overlaps.is_empty()
    }

    /// Total number of coincidences across all levels.
    pub fn total(&self) -> usize {
        self.vertex_merges.len() + self.edge_overlaps.len() + self.face_overlaps.len()
    }

    /// Look up the coincidence orientation for a face pair, if any.
    pub fn face_orientation(
        &self,
        face_a: FaceHandle,
        face_b: FaceHandle,
    ) -> Option<CoincidenceOrientation> {
        self.face_overlaps
            .iter()
            .find(|fo| {
                (fo.face_a == face_a && fo.face_b == face_b)
                    || (fo.face_a == face_b && fo.face_b == face_a)
            })
            .map(|fo| fo.orientation)
    }
}
