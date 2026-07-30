//! Boundary-edge and retired-span records for the planned
//! repeated-apex partition (split from the plan module, 500-line cap).

use crate::apex_construction::repeated_apex_loop_plan::RepeatedApexArcKind;
use popkern_core::topology::{EdgeHandle, HalfEdgeHandle, VertexHandle};

/// One piece of the partition's intended cap boundary.
pub(in crate::apex_construction) enum PartitionBoundaryEdge {
    /// FilletEnd seam between a sector's two T-vertices.
    FilletEnd {
        /// Sector index.
        fillet_index: usize,
        /// Sector start vertex T_a[k].
        t_a: VertexHandle,
        /// Sector end vertex T_b[k].
        t_b: VertexHandle,
    },
    /// Bridged corner FaceEdge (absent for welded shared corners).
    FaceEdge {
        /// Corner index.
        fillet_index: usize,
        /// The bridge edge.
        edge: EdgeHandle,
        /// Cap-side half-edge of the bridge.
        cap_half_edge: HalfEdgeHandle,
        /// Support-side half-edge of the bridge.
        neighbor_half_edge: HalfEdgeHandle,
    },
}

/// A retired arc of the pre-normalization cap loop.
pub(in crate::apex_construction) struct RetiredApexSpan {
    /// Arc index in the retiring loop walk.
    pub(in crate::apex_construction) arc_index: usize,
    /// Arc classification.
    pub(in crate::apex_construction) kind: RepeatedApexArcKind,
    /// The arc's half-edges in walk order.
    pub(in crate::apex_construction) half_edges: Vec<HalfEdgeHandle>,
}

impl PartitionBoundaryEdge {
    /// One-line diagnostic form for partition trace summaries.
    pub(in crate::apex_construction) fn summary(&self) -> String {
        match self {
            PartitionBoundaryEdge::FilletEnd {
                fillet_index,
                t_a,
                t_b,
            } => format!("FilletEnd[{fillet_index}]({t_a:?}->{t_b:?})"),
            PartitionBoundaryEdge::FaceEdge {
                fillet_index,
                edge,
                cap_half_edge,
                neighbor_half_edge,
            } => format!(
                "FaceEdge[{fillet_index}]({edge:?}, cap={cap_half_edge:?}, \
                 neighbor={neighbor_half_edge:?})"
            ),
        }
    }
}
