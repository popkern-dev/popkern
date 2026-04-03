//! Result types for the pairwise intersection hierarchy.
//!
//! The hierarchy computes intersections bottom-up: VV → VE → EE → VF → EF,
//! producing a [`HierarchyResult`] that is consumed by the FF intersection
//! phase and by the imprinting/merge pipeline.

use popkern_core::geometry::Point3;
use popkern_core::topology::{EdgeHandle, FaceHandle, VertexHandle};

/// A detected vertex-vertex coincidence between two bodies.
#[derive(Debug, Clone)]
pub struct VertexVertexHit {
    /// Vertex from body A.
    pub vertex_a: VertexHandle,
    /// Vertex from body B.
    pub vertex_b: VertexHandle,
    /// Distance between the two vertices.
    pub distance: f64,
}

/// A detected vertex-on-edge contact.
#[derive(Debug, Clone)]
pub struct VertexEdgeHit {
    /// The vertex (from one body).
    pub vertex: VertexHandle,
    /// Which body the vertex belongs to.
    pub vertex_body: Operand,
    /// The edge it lies on (from the other body).
    pub edge: EdgeHandle,
    /// Parameter t on the edge curve at the contact point.
    pub t_edge: f64,
    /// Distance from vertex to the closest point on the edge.
    pub distance: f64,
}

/// A detected edge-edge intersection point.
#[derive(Debug, Clone)]
pub struct EdgeEdgeHit {
    /// Edge from body A.
    pub edge_a: EdgeHandle,
    /// Edge from body B.
    pub edge_b: EdgeHandle,
    /// Parameter on edge A's curve.
    pub t_a: f64,
    /// Parameter on edge B's curve.
    pub t_b: f64,
    /// The 3D intersection point.
    pub point: Point3,
}

/// A detected vertex-on-face contact.
#[derive(Debug, Clone)]
pub struct VertexFaceHit {
    /// The vertex.
    pub vertex: VertexHandle,
    /// Which body the vertex belongs to.
    pub vertex_body: Operand,
    /// The face it lies on.
    pub face: FaceHandle,
    /// UV parameters on the face surface.
    pub uv: (f64, f64),
    /// Distance from vertex to the surface.
    pub distance: f64,
}

/// A detected edge-face intersection (pierce point).
#[derive(Debug, Clone)]
pub struct EdgeFaceHit {
    /// The edge.
    pub edge: EdgeHandle,
    /// Which body the edge belongs to.
    pub edge_body: Operand,
    /// The face it pierces.
    pub face: FaceHandle,
    /// The 3D intersection point.
    pub point: Point3,
    /// Parameter t on the edge curve.
    pub t_edge: f64,
    /// UV parameters on the face surface.
    pub uv: (f64, f64),
}

/// Identifies which operand body an entity belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operand {
    /// The first (left-hand) operand.
    A,
    /// The second (right-hand) operand.
    B,
}

/// Complete result of the sub-FF intersection hierarchy.
///
/// Computed bottom-up before FF intersection. Lower-dimensional hits
/// constrain and inform the FF phase (e.g. VE hits provide starting
/// points for intersection curve tracing, EF hits provide edge split
/// points for imprinting).
#[derive(Debug, Clone, Default)]
pub struct HierarchyResult {
    /// Vertex-vertex coincidences.
    pub vv_hits: Vec<VertexVertexHit>,
    /// Vertex-on-edge contacts.
    pub ve_hits: Vec<VertexEdgeHit>,
    /// Edge-edge intersection points.
    pub ee_hits: Vec<EdgeEdgeHit>,
    /// Vertex-on-face contacts.
    pub vf_hits: Vec<VertexFaceHit>,
    /// Edge-face pierce points.
    pub ef_hits: Vec<EdgeFaceHit>,
}

impl HierarchyResult {
    /// Total number of sub-FF intersection hits found.
    pub fn total_hits(&self) -> usize {
        self.vv_hits.len()
            + self.ve_hits.len()
            + self.ee_hits.len()
            + self.vf_hits.len()
            + self.ef_hits.len()
    }

    /// Returns true if no sub-FF intersections were found.
    pub fn is_empty(&self) -> bool {
        self.total_hits() == 0
    }
}
