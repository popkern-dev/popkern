//! Dof copies of grid nodes around split cells: each node gets one copy
//! per group of entries around it that touch through face neighbours
//! within the node's eight cells.

use std::collections::BTreeSet;
use std::ops::Range;

use crate::grid::Grid;

use super::{FragmentError, Fragments, NODE_COPIES, find};

/// Sets `out.corner_copy` for the entries of every node of the cells in
/// `need` (`entry_range` gives a cell's entries), from `out.pair_touch`.
pub(super) fn assign(
    grid: &Grid,
    need: &BTreeSet<usize>,
    entry_range: &dyn Fn(usize) -> Range<usize>,
    out: &mut Fragments,
) -> Result<(), FragmentError> {
    let mut nodes: BTreeSet<usize> = BTreeSet::new();
    for &c in need {
        if let Some(ijk) = grid.cell_ijk(c) {
            nodes.extend(grid.cell_nodes(ijk).unwrap_or_default());
        }
    }
    for node in nodes {
        let Some(nijk) = grid.node_ijk(node) else {
            continue;
        };
        let mut star: Vec<(usize, [usize; 3])> = Vec::new();
        for o in crate::grid::CELL_CORNERS {
            let c = [0, 1, 2].map(|a| nijk[a].checked_sub(o[a]));
            if let [Some(i), Some(j), Some(k)] = c
                && let Some(cell) = grid.cell_index([i, j, k])
            {
                for e in entry_range(cell) {
                    star.push((e, [i, j, k]));
                }
            }
        }
        let mut parent: Vec<usize> = (0..star.len()).collect();
        for x in 0..star.len() {
            for y in x + 1..star.len() {
                let (ex, cx) = star[x];
                let (ey, cy) = star[y];
                let diff: usize = (0..3).map(|a| cx[a].abs_diff(cy[a])).sum();
                if diff != 1 {
                    continue;
                }
                let key = (ex.min(ey), ex.max(ey));
                let touch = out
                    .pair_touch
                    .binary_search_by_key(&key, |p| p.0)
                    .map_or(true, |i| out.pair_touch[i].1);
                if touch {
                    let (rx, ry) = (find(&mut parent, x), find(&mut parent, y));
                    if rx != ry {
                        parent[rx.max(ry)] = rx.min(ry);
                    }
                }
            }
        }
        // Copies in order of each group's first entry.
        let mut roots: Vec<usize> = Vec::new();
        for x in 0..star.len() {
            let r = find(&mut parent, x);
            if !roots.contains(&r) {
                roots.push(r);
            }
        }
        if roots.len() > NODE_COPIES {
            return Err(FragmentError::TooMany(node));
        }
        for (x, &(e, cijk)) in star.iter().enumerate() {
            let r = find(&mut parent, x);
            let copy = roots.iter().position(|&q| q == r).unwrap_or(0) as u8;
            let corner = crate::grid::CELL_CORNERS
                .iter()
                .position(|o| (0..3).all(|a| cijk[a] + o[a] == nijk[a]))
                .unwrap_or(0);
            out.corner_copy[e][corner] = copy;
        }
    }
    Ok(())
}
