//! Non-manifold sewing (PK-1109, S10): each pair glued into one ring
//! in radial order (`BRep::glue_edges`), so a third face sewn onto an
//! edge already joined makes a ring of three instead of declining.

use super::brep_sewing::{
    BrepSewingDeclineReason, BrepSewingEdgePair, BrepSewingOptions, edge_segment,
    reversed_endpoint_gap,
};
use crate::topology::{BRep, EdgeHandle};

/// Glue `pair` (its edges resolved through earlier glues) into one ring
/// on a working copy, committed when every step succeeds. Returns the
/// endpoint gap.
pub(super) fn glue_pair(
    brep: &mut BRep,
    pair: &BrepSewingEdgePair,
    options: BrepSewingOptions,
    alias: &mut std::collections::BTreeMap<EdgeHandle, EdgeHandle>,
) -> Result<f64, BrepSewingDeclineReason> {
    let resolve = |alias: &std::collections::BTreeMap<EdgeHandle, EdgeHandle>, mut e| {
        while let Some(&next) = alias.get(&e) {
            e = next;
        }
        e
    };
    let (keep, gone) = (resolve(alias, pair.first), resolve(alias, pair.second));
    if keep == gone {
        return Err(BrepSewingDeclineReason::SameEdge);
    }
    for e in [keep, gone] {
        if !brep.edges.contains_key(e) {
            return Err(BrepSewingDeclineReason::MissingEdge);
        }
        if brep.radial_fins(e).is_empty() {
            return Err(BrepSewingDeclineReason::NotBoundaryEdge);
        }
    }
    let (Some(a), Some(b)) = (edge_segment(brep, keep), edge_segment(brep, gone)) else {
        return Err(BrepSewingDeclineReason::MissingEdge);
    };
    // Either way round: a third face may run the edge either way.
    let straight = a.0.distance_to(&b.0).max(a.1.distance_to(&b.1));
    let gap = reversed_endpoint_gap(a, b).min(straight);
    if gap > options.tolerance {
        return Err(BrepSewingDeclineReason::EndpointGapExceeded {
            gap,
            tolerance: options.tolerance,
        });
    }
    let mut work = brep.clone();
    let mut vertices = Default::default();
    let glued = work
        .glue_edge_ends(keep, gone, &mut vertices)
        .and_then(|()| work.glue_edges(keep, gone).map(|_| ()));
    if let Err(e) = glued {
        return Err(BrepSewingDeclineReason::GlueDeclined {
            reason: e.to_string(),
        });
    }
    *brep = work;
    alias.insert(gone, keep);
    Ok(gap)
}

#[cfg(test)]
#[path = "brep_sewing_glue_tests.rs"]
mod tests;
