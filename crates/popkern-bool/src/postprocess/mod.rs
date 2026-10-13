//! Post-processing for boolean results: redundant edge removal and face simplification.

pub mod collinear_vertex;
pub(crate) mod curveless_edge;
pub mod exact_rims;
pub mod face_simplify;
pub(crate) mod junction_reconcile;
pub mod keyhole_loop;
pub mod loop_repair;
pub(crate) mod merge_rechart;
pub mod parallel_weld;
pub mod point_loop_retire;
pub mod quality;
pub(crate) mod redundant_chain;
pub mod redundant_edge;
pub mod self_twin_stub;
pub mod slit_retire;
pub mod slit_seam;
pub mod sliver_removal;
pub mod spur_edge;

#[cfg(test)]
#[path = "postprocess_tests.rs"]
mod tests;

/// Closed-edge sense at the producer exit (PK-478's repair arm), after
/// every pcurve author and the orientation reconcile its verdicts read:
/// a bored tube's start-cap hole pcurve ran with its outer loop (the
/// cap read outer + hole, PK-819) and a cone-through-cube union's apex
/// rim served canonical-backward once the late slit sweep had run. The
/// repair reverses only a chain-convicted pcurve or swaps a backward
/// canonical side; it never acts on an indeterminate verdict.
pub(crate) fn repair_closed_edge_sense(result: &mut popkern_core::topology::BRep) {
    let repaired = result.repair_canonical_forward();
    if repaired > 0 && popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "[boolean] post-processing: {repaired} closed edge(s) repaired to canonical-forward"
        );
    }
}
