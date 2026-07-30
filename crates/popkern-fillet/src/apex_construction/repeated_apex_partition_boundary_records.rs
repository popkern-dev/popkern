//! Boundary-edge and retired-span records for the planned
//! repeated-apex partition (split from the plan module, 500-line cap).

use crate::apex_construction::repeated_apex_loop_plan::RepeatedApexArcKind;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};

use crate::apex_construction::error::ApexBuildError;

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

/// Retired-arc ownership for a WELDED shared corner (PK-219): the
/// retiring cap loop passes `V -> C -> V` through the single welded
/// corner vertex, so the arc has exactly two half-edge legs and no
/// bridge edge.
#[derive(Clone)]
pub(in crate::apex_construction) struct SharedCornerOwnership {
    /// Retiring leg `V -> C` (plays the FaceEdge arc's `apex_to_t_b`
    /// role for the corner's sector).
    pub(in crate::apex_construction) apex_to_corner: HalfEdgeHandle,
    /// Retiring leg `C -> V` (plays `t_a_next_to_apex` for sector
    /// `corner_index + 1`).
    pub(in crate::apex_construction) corner_to_apex: HalfEdgeHandle,
}

impl SharedCornerOwnership {
    /// Build from a classified SharedCorner retired span.
    pub(in crate::apex_construction) fn from_retired_span(
        brep: &BRep,
        span: &RetiredApexSpan,
        corner_index: usize,
    ) -> Result<Self, ApexBuildError> {
        let [apex_to_corner, corner_to_apex] = span.half_edges.as_slice() else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition expected welded SharedCorner[{corner_index}] \
                     arc to have 2 half-edges (V->C, C->V), got {}",
                    span.half_edges.len()
                ),
            });
        };
        for half_edge in [apex_to_corner, corner_to_apex] {
            if !brep.half_edges.contains_key(*half_edge) {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition welded SharedCorner[{corner_index}] \
                         retired leg {half_edge:?} is stale"
                    ),
                });
            }
        }
        Ok(Self {
            apex_to_corner: *apex_to_corner,
            corner_to_apex: *corner_to_apex,
        })
    }
}

/// Per-corner retired-arc ownership: a bridged corner owns a
/// FaceEdge triple; a welded corner owns its two retiring legs.
#[derive(Clone)]
pub(in crate::apex_construction) enum CornerArcOwnership {
    /// Bridged corner (historical 2N contract).
    Bridged(super::super::partition_ownership::FaceEdgeOwnership),
    /// Welded shared corner (PK-219).
    Welded(SharedCornerOwnership),
}

impl CornerArcOwnership {
    /// Retiring leg that plays `apex_to_t_b` for this corner's sector.
    pub(in crate::apex_construction) fn apex_to_t_b(&self) -> HalfEdgeHandle {
        match self {
            CornerArcOwnership::Bridged(ownership) => ownership.apex_to_t_b,
            CornerArcOwnership::Welded(ownership) => ownership.apex_to_corner,
        }
    }

    /// Retiring leg that plays `t_a_next_to_apex` for the NEXT sector.
    pub(in crate::apex_construction) fn t_a_next_to_apex(&self) -> HalfEdgeHandle {
        match self {
            CornerArcOwnership::Bridged(ownership) => ownership.t_a_next_to_apex,
            CornerArcOwnership::Welded(ownership) => ownership.corner_to_apex,
        }
    }

    /// Cap-cycle half-edge contributed by this corner (bridged only).
    pub(in crate::apex_construction) fn cap_cycle_half_edge(&self) -> Option<HalfEdgeHandle> {
        match self {
            CornerArcOwnership::Bridged(ownership) => Some(ownership.face_edge_cap_half_edge),
            CornerArcOwnership::Welded(_) => None,
        }
    }
}
