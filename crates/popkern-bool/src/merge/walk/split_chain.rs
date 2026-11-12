//! The loops of a face split from its outer loop by records: those
//! reached from the outer through shared edges (PK-1173).
//!
//! The size election files a loop as a hole when it is under half the
//! largest loop's area and shares no edge with the elected outer. A
//! face cut by two records near one of its edges is three strips: the
//! outer shares the nearer record with the middle strip, which shares
//! the other with the strip at the edge, a loop sharing nothing with
//! the outer itself and filed as its hole (a torus less a slab 0.001
//! to 0.011 off its meridian edges served a strip 2e-4 wide in its
//! chart as an inner loop winding with its outer, and declined). Every
//! strip of such a chain is a region of the face. A genuine hole shares
//! edges only with the interior piece it bounds, which shares none
//! with the outer's chain, so the chain never reaches it.

use std::collections::HashSet;

use super::super::graph::IntersectionGraph;

/// For each of `loops`, whether it is reached from the outer loop
/// (`outer_hes`) through half-edges whose twins lie in the outer or in
/// a loop already reached.
pub(super) fn split_from_outer(
    graph: &IntersectionGraph,
    outer_hes: &HashSet<usize>,
    loops: &[Vec<usize>],
) -> Vec<bool> {
    let mut chain = outer_hes.clone();
    let mut reached = vec![false; loops.len()];
    loop {
        let mut grew = false;
        for (k, lp) in loops.iter().enumerate() {
            if reached[k] {
                continue;
            }
            if lp
                .iter()
                .any(|&h| chain.contains(&graph.half_edges[h].twin))
            {
                reached[k] = true;
                chain.extend(lp.iter().copied());
                grew = true;
            }
        }
        if !grew {
            return reached;
        }
    }
}
