//! Redundant-edge dissolution for the serve postprocess (split from
//! `operation_postprocess.rs` under the file-size cap).

use popkern_core::{tolerance::Tolerance, topology::BRep};

/// Safely remove redundant edges, validating after each removal.
///
/// Finds redundant edges and removes them one at a time via KEF.
/// After each removal, validates the topology. If validation fails,
/// stops immediately (leaving remaining redundant edges in place).
/// Returns how many edits it made (KEFs applied, chains dissolved) —
/// zero when the body is unchanged.
pub(crate) fn remove_redundant_edges_safe(result: &mut BRep, tolerance: &Tolerance) -> usize {
    use crate::postprocess::redundant_edge::find_redundant_edges;

    let report = find_redundant_edges(result, tolerance);

    if report.redundant_edges.is_empty() {
        return dissolve_redundant_chains(result, tolerance);
    }

    let (mut removed, mut edits) = (0, 0usize);
    for eh in &report.redundant_edges {
        // Re-check against the CURRENT topology: the report is
        // computed once, and earlier KEFs in this loop change both
        // which edges two faces share (a stale single-share verdict
        // dissolves one edge of a multi-edge boundary and leaves the
        // survivor self-twinned) and how far a merged face spans the
        // chart (a merge that was span-safe at report time can wrap
        // a period by the time its edge comes up).
        if !crate::postprocess::redundant_edge::edge_redundant_now(result, *eh, tolerance) {
            continue;
        }
        if !crate::clip::chart_cell::edge_dissolve_keeps_cell_resolvable(result, *eh, tolerance) {
            continue;
        }
        let he = match result.edges.get(*eh) {
            Some(e) => e.half_edge,
            None => continue,
        };

        // Snapshot vertex/edge/face counts before KEF
        let pre_v = result.vertices.len();
        let pre_e = result.edges.len();
        let pre_f = result.faces.len();

        let sides = crate::postprocess::merge_rechart::merge_sides(result, he);
        if let Ok(kef) = result.kef(he) {
            edits += 1;
            crate::postprocess::merge_rechart::rechart_inherited(result, &sides, kef.removed_face);
            // Validate after removal — if topology is corrupted, stop
            if result.validate().is_err() {
                // Revert is impossible, but at least stop making it worse
                if cfg!(debug_assertions) {
                    eprintln!(
                        "[boolean] KEF corrupted topology after removing edge, stopping. \
                         Before: V={pre_v} E={pre_e} F={pre_f}"
                    );
                }
                break;
            }
            removed += 1;
        }
    }

    if removed > 0 && cfg!(debug_assertions) {
        eprintln!("[boolean] post-processing: removed {removed} redundant edge(s)");
    }
    edits + dissolve_redundant_chains(result, tolerance)
}

/// Merge co-surfacial face pairs whose shared boundary is a multi-edge
/// PATH (PK-379): staged per pair, kept only when the body still
/// validates. Single-edge pairs were handled above; pairs whose
/// shared edges are not one path stay split (a ring would need an
/// inner loop, not a merge).
fn dissolve_redundant_chains(result: &mut BRep, tolerance: &Tolerance) -> usize {
    use crate::postprocess::redundant_chain::{dissolve_chain, shared_coplanar_chain};
    let mut merged = 0usize;
    loop {
        let mut progressed = false;
        let edges: Vec<_> = result.edges.keys().collect();
        for eh in edges {
            let Some(chain) = shared_coplanar_chain(result, eh, tolerance) else {
                continue;
            };
            if chain.len() < 2 {
                continue;
            }
            if !crate::clip::chart_cell::edge_dissolve_keeps_cell_resolvable(
                result, chain[0], tolerance,
            ) {
                continue;
            }
            let mut staged = result.clone();
            if dissolve_chain(&mut staged, &chain).is_ok() && staged.validate().is_ok() {
                *result = staged;
                merged += 1;
                progressed = true;
                break;
            }
        }
        if !progressed {
            break;
        }
    }
    if merged > 0 && cfg!(debug_assertions) {
        eprintln!("[boolean] post-processing: dissolved {merged} redundant shared chain(s)");
    }
    merged
}
