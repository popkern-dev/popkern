//! Target topology spec for Phase 6c-A2 apex construction.
//!
//! `TargetTopology` is a **declarative** description of the desired
//! final apex region: vertices to ensure exist, edges with their
//! curve geometry, and the cap face spec.  The builder
//! ([`construct_apex`](super::construct_apex)) is responsible for
//! choosing the Euler op sequence to realise it.

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{FaceHandle, VertexHandle};

/// Declarative spec of the desired final apex region.
///
/// The planner emits this; the builder consumes it and applies
/// Euler ops to realise it.  `TargetTopology` is INTROSPECTABLE
/// (UI can render it as a preview before construction).
#[derive(Debug, Clone)]
pub struct TargetTopology {
    /// Vertices that should exist at these positions after
    /// construction.  Existing vertices within tolerance are
    /// reused; missing vertices are created via `mev` (or the
    /// bulk equivalent).
    pub apex_vertices: Vec<Point3>,
    /// Edges that should exist between vertex pairs with these
    /// curve geometries.
    pub apex_edges: Vec<TargetEdge>,
    /// Cap face spec — single Sphere or N Gregory sub-patches.
    pub cap_face: TargetCapFace,
    /// Updates to existing non-cap face boundaries (cube/pyramid
    /// face boundaries truncated at the apex region).  Each entry
    /// describes one face whose outer loop must be modified to
    /// match the new apex topology.
    pub face_boundary_updates: Vec<FaceBoundaryUpdate>,
}

/// One target edge in the final apex topology.
///
/// Specifies the edge's endpoints and the curve3d that should be
/// stored on it.  When [`construct_apex`](super::construct_apex)
/// realises this, it creates the edge via `mef` (or finds an
/// existing edge between the endpoints) and stores the curve
/// directly (no separate wire step).
#[derive(Debug, Clone)]
pub struct TargetEdge {
    /// Index into `TargetTopology::apex_vertices` for the edge's
    /// origin vertex (in BRep direction sense — primary HE origin).
    pub origin_vertex_idx: usize,
    /// Index into `TargetTopology::apex_vertices` for the edge's
    /// dest vertex.
    pub dest_vertex_idx: usize,
    /// The curve3d to store on this edge.  Endpoints in BRep
    /// direction: `curve.point_at(curve.parameter_range().0)`
    /// must equal the origin vertex's position;
    /// `curve.point_at(curve.parameter_range().1)` the dest's.
    pub curve: NurbsCurve,
    /// Edge role classification — used by the builder to pick
    /// face memberships and pcurve construction.
    pub role: TargetEdgeRole,
}

/// Edge role classification for `TargetEdge`.
///
/// The builder uses this to decide which surface(s) the edge
/// should be associated with (cap face, cylinder face, plane face).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TargetEdgeRole {
    /// Cap-vs-Cyl edge — typically the FilletEnd arc terminating
    /// on the fillet's cylinder.  Identifies which fillet.
    CapVsCyl {
        /// Fillet index this edge sits on.
        fillet_idx: usize,
    },
    /// Cap-vs-Plane edge — typically a FaceEdge bridge connecting
    /// adjacent FilletEnd terminations on a shared face.
    /// Identifies the shared face index.
    CapVsPlane {
        /// Index into the face_boundary_updates / shared-face list.
        shared_face_idx: usize,
    },
    /// Inner cap-vs-cap spoke — only relevant for Gregory caps with
    /// 3-quad poke decomposition; not used for Sphere caps.
    InnerSpoke,
}

/// Cap face spec for the target topology.
///
/// For Sphere caps: a single face with Sphere surface assigned.
/// For Gregory caps: an N-sided face that the cap-construction
/// pipeline will subsequently decompose into N sub-patches via
/// 3-quad poke.
#[derive(Debug, Clone)]
pub enum TargetCapFace {
    /// Single Sphere cap face.  Outer loop walks the apex_vertices
    /// in cap-loop order, all edges are CapVsCyl.
    Sphere {
        /// Cap-sphere centre.
        center: Point3,
        /// Cap-sphere radius.
        radius: f64,
        /// Outer-loop edge indices into `TargetTopology::apex_edges`,
        /// in cap-loop walking order.
        outer_loop_edge_indices: Vec<usize>,
    },
    /// N-sided Gregory cap.  Outer loop has 2N edges alternating
    /// CapVsCyl and CapVsPlane.
    Gregory {
        /// Outer-loop edge indices into `TargetTopology::apex_edges`,
        /// in cap-loop walking order (2N entries).
        outer_loop_edge_indices: Vec<usize>,
    },
}

/// Boundary update for an existing non-cap face.
///
/// When the apex topology is built, existing face boundaries (e.g.
/// cube face, fillet cylinder face) need their outer-loop chains
/// updated to connect to the new apex vertices and replace any
/// stale apex-side curves.
#[derive(Debug, Clone)]
pub struct FaceBoundaryUpdate {
    /// The face whose boundary needs updating.
    pub face: FaceHandle,
    /// The vertex on this face's outer loop that should be replaced
    /// (if any).  Typically the rolling-ball footprint vertex that
    /// no longer exists.
    pub remove_vertex: Option<VertexHandle>,
    /// New apex vertex indices (in `TargetTopology::apex_vertices`)
    /// that the face's outer loop should walk through after the
    /// update, in walk order.
    pub insert_vertex_indices: Vec<usize>,
    /// Apex-side base edge curves to update (cyl-vs-Plane apex
    /// generators that previously ended at `remove_vertex`).
    /// Each entry is (edge handle, new curve).
    pub update_edge_curves: Vec<(popkern_core::topology::EdgeHandle, NurbsCurve)>,
}
