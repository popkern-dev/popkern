//! Why a non-manifold Euler operator declined.

use crate::topology::{
    EdgeHandle, FaceHandle, HalfEdgeHandle, RadialOrderError, RegionHandle, ShellHandle,
    VertexHandle,
};

/// Why a non-manifold Euler operator declined. The body is unchanged.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum EulerError {
    /// A vertex is not in the body.
    #[error("vertex {0:?} is not in the body")]
    UnknownVertex(VertexHandle),
    /// An edge is not in the body.
    #[error("edge {0:?} is not in the body")]
    UnknownEdge(EdgeHandle),
    /// A region is not in the body.
    #[error("region {0:?} is not in the body")]
    UnknownRegion(RegionHandle),
    /// The vertex is not an acorn (an isolated vertex held by a shell).
    #[error("vertex {vertex:?} is not an acorn")]
    NotAnAcorn {
        /// The vertex.
        vertex: VertexHandle,
    },
    /// The edge is not a wire edge held by a shell.
    #[error("edge {edge:?} is not a held wire edge")]
    NotAWireEdge {
        /// The edge.
        edge: EdgeHandle,
    },
    /// A wire edge from a vertex to itself.
    #[error("a wire edge from vertex {vertex:?} to itself")]
    DegenerateWireEdge {
        /// The vertex.
        vertex: VertexHandle,
    },
    /// The vertex lies in no shell of the region: it is held by another
    /// region's shell, or only faces whose sides face other regions
    /// reach it.
    #[error("vertex {vertex:?} does not lie in region {region:?}")]
    VertexNotInRegion {
        /// The vertex.
        vertex: VertexHandle,
        /// The region asked for.
        region: RegionHandle,
    },
    /// The wire would join two shells that both hold face sides: a
    /// shell's face sides are one surface connected through its edges,
    /// so faces joined only by a wire stay separate shells.
    #[error("a wire edge would join face shells {0:?} and {1:?}")]
    WireJoinsFaceShells(ShellHandle, ShellHandle),
    /// The two vertices to glue are one.
    #[error("vertex {vertex:?} cannot be glued to itself")]
    SameVertex {
        /// The vertex.
        vertex: VertexHandle,
    },
    /// The entities to glue lie further apart than the tolerance.
    #[error("not coincident: {distance:e} apart, tolerance {tolerance:e}")]
    NotCoincident {
        /// The measured distance (the largest, for curves).
        distance: f64,
        /// The tolerance it was judged against.
        tolerance: f64,
    },
    /// The two edges to glue are one.
    #[error("edge {edge:?} cannot be glued to itself")]
    SameEdge {
        /// The edge.
        edge: EdgeHandle,
    },
    /// The edge has no face use (glue faces' edges; a wire edge lying
    /// along a face edge goes with `kill_wire_edge`).
    #[error("edge {edge:?} has no face use")]
    NotAFaceEdge {
        /// The edge.
        edge: EdgeHandle,
    },
    /// The edges to glue do not end at the same two vertices.
    #[error(
        "edges {keep:?} and {gone:?} do not share their end vertices (glue the vertices first)"
    )]
    EdgeEndsDiffer {
        /// The edge kept.
        keep: EdgeHandle,
        /// The edge glued away.
        gone: EdgeHandle,
    },
    /// The merged ring's radial order is undecided.
    #[error("the glued ring's radial order is undecided: {0}")]
    RadialOrder(Box<RadialOrderError>),
    /// A wedge of the glued ring would lie between two regions: the
    /// faces cross there, which gluing does not decide.
    #[error("gluing edge {edge:?} puts the wedge between fins {fins:?} in regions {regions:?}")]
    RegionConflict {
        /// The edge kept.
        edge: EdgeHandle,
        /// The wedge's fins, in ring order.
        fins: (HalfEdgeHandle, HalfEdgeHandle),
        /// The regions of its two sides.
        regions: (RegionHandle, RegionHandle),
    },
    /// The fins are not a proper consecutive run of the edge's ring.
    #[error("the fins are not a consecutive run of edge {edge:?}'s ring")]
    NotARingRun {
        /// The edge.
        edge: EdgeHandle,
    },
    /// Ungluing would leave a bounded region's shell in pieces.
    #[error("ungluing would split the boundary of region {region:?}")]
    UnglueSplitsRegion {
        /// The region.
        region: RegionHandle,
    },
    /// The faces' boundaries do not pair fin for fin, each running the
    /// other way between the same vertices along a coincident edge.
    #[error("faces {keep:?} and {gone:?} do not share their boundary")]
    FacesDoNotMatch {
        /// The face kept.
        keep: FaceHandle,
        /// The face glued away.
        gone: FaceHandle,
    },
    /// No side of the one face lies in the region of a side of the
    /// other: nothing lies between them to close.
    #[error("faces {keep:?} and {gone:?} do not face each other across one region")]
    FacesNotFacing {
        /// The face kept.
        keep: FaceHandle,
        /// The face glued away.
        gone: FaceHandle,
    },
    /// The face does not lie between two solids in different shells (a
    /// lamina, a sheet's face, a solid's boundary face), or a fin's sense
    /// along its edge is undecidable.
    #[error("face {face:?} is not an internal face between two solids")]
    NotAnInternalFace {
        /// The face.
        face: FaceHandle,
    },
    /// The shell does not hold both sides of a closed surface alone.
    #[error("shell {shell:?} does not hold a closed surface on both sides")]
    NotAClosedSurface {
        /// The shell.
        shell: ShellHandle,
    },
    /// The region is not bounded by one surface with another region's
    /// shell across every face of it.
    #[error("region {region:?} is not a face region")]
    NotAFaceRegion {
        /// The region.
        region: RegionHandle,
    },
    /// The general body has no manifold pieces to split into.
    #[error("not splittable into manifold pieces: {reason}")]
    NotSplittable {
        /// What stands in the way.
        reason: &'static str,
    },
    /// The fins given are not a closed cycle of distinct wire edges
    /// held by one shell.
    #[error("not a wire cycle: {reason}")]
    NotAWireCycle {
        /// What the cycle breaks.
        reason: &'static str,
    },
    /// The face is not a lamina spanning a wire cycle: it has holes, its
    /// sides lie in different shells, or an edge of its loop has another
    /// use.
    #[error("face {face:?} is not a lamina on free edges")]
    NotALaminaFace {
        /// The face.
        face: FaceHandle,
    },
}
