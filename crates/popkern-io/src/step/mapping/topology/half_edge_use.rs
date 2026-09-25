//! Which of an edge's two half-edges a STEP oriented edge uses.

use popkern_core::topology::brep::BRep;
use popkern_core::topology::{HalfEdgeHandle, LoopHandle};

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::EntityId;

/// The half-edge an oriented edge takes: the one its orientation
/// prefers, or its twin when the preferred one is already wired into a
/// loop (one face use per half-edge), with a warning. Returns the
/// half-edge and whether it runs along the edge.
pub(super) fn choose_step_half_edge_use(
    brep: &BRep,
    he_fwd: HalfEdgeHandle,
    he_rev: HalfEdgeHandle,
    oriented_forward: bool,
    edge_curve_id: EntityId,
    log: &mut DiagnosticLog,
) -> (HalfEdgeHandle, bool) {
    let preferred = if oriented_forward { he_fwd } else { he_rev };
    let alternate = if oriented_forward { he_rev } else { he_fwd };
    if !half_edge_is_wired(brep, preferred) {
        return (preferred, preferred == he_fwd);
    }
    if !half_edge_is_wired(brep, alternate) {
        log.warn(
            Phase::MapTopology,
            Some(edge_curve_id),
            "STEP oriented edge preferred an already wired half-edge; using available twin to preserve one face-use per half-edge",
        );
        return (alternate, alternate == he_fwd);
    }
    (preferred, preferred == he_fwd)
}

fn half_edge_is_wired(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    brep.half_edges
        .get(half_edge)
        .is_some_and(|data| data.loop_ != LoopHandle::default())
}
