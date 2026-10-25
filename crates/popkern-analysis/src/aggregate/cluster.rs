//! Clusters of extended cells: the free dofs that cells whose nodes are
//! all extended (no dof of their own) are interpolated from, joined
//! where such cells share dofs, with those cells and their roots.
//!
//! On a cluster the field is a combination of its dofs alone, so a
//! bound between quadratic forms over its cells (the Nitsche flux on
//! its fixed pieces against its energy, PH-067) holds for any field.
//! A fin thinner than the lattice spacing along its whole length is one
//! cluster: each layer of fin cells extends from its own root, and the
//! cells across two layers' nodes join them.

use std::collections::BTreeMap;

use crate::grid::Grid;
use crate::moments::SolidMoments;
use crate::moments::fragments::NODE_COPIES;

use super::{Aggregation, Row};

fn find(p: &mut [usize], mut x: usize) -> usize {
    while p[x] != x {
        p[x] = p[p[x]];
        x = p[x];
    }
    x
}

/// A node's value as (position in a cluster's dofs, weight) terms.
pub type NodeTerms = Vec<(usize, f64)>;

/// One cluster.
#[derive(Debug, Clone, PartialEq)]
pub struct Cluster {
    /// Its free dofs, ascending.
    pub dofs: Vec<u32>,
    /// Its cells: the entry and, per node (`CELL_CORNERS` order), the
    /// node's value as (position in `dofs`, weight) terms.
    pub cells: Vec<(usize, [NodeTerms; 8])>,
}

impl Cluster {
    /// Physical position of each of its dofs' nodes.
    pub fn positions(&self, grid: &Grid, aggregation: &Aggregation) -> Vec<[f64; 3]> {
        let e = &aggregation.extension;
        let mut key_of: BTreeMap<u32, usize> = BTreeMap::new();
        for (key, row) in e.nodes.iter().zip(&e.rows) {
            if let Row::Free(d) = row {
                key_of.insert(*d, *key);
            }
        }
        self.dofs
            .iter()
            .map(|d| {
                let key = key_of.get(d).copied().unwrap_or(0);
                grid.node_position(grid.node_ijk(key / NODE_COPIES).unwrap_or_default())
            })
            .collect()
    }
}

/// The clusters of `aggregation`, in order of their least dof.
pub fn clusters(aggregation: &Aggregation, grid: &Grid, moments: &SolidMoments) -> Vec<Cluster> {
    let e = &aggregation.extension;
    let row = |key: usize| e.position(key).and_then(|p| e.rows.get(p));
    // Each node's value as (free dof, weight) terms.
    let terms = |key: usize| -> Vec<(u32, f64)> {
        match row(key) {
            Some(Row::Free(d)) => vec![(*d, 1.0)],
            Some(Row::Extended(w)) => w.to_vec(),
            None => Vec::new(),
        }
    };
    // Cells with every node extended, and their roots.
    let mut cells: Vec<usize> = Vec::new();
    for entry in 0..moments.cells.len() {
        let Some(keys) = moments.node_keys(grid, entry) else {
            continue;
        };
        if keys
            .iter()
            .all(|&k| matches!(row(k), Some(Row::Extended(_))))
        {
            cells.push(entry);
            if let Some(&r) = aggregation.root.get(entry) {
                cells.push(r);
            }
        }
    }
    cells.sort_unstable();
    cells.dedup();
    let mut parent: Vec<usize> = (0..e.free).collect();
    let mut node_terms: Vec<[Vec<(u32, f64)>; 8]> = Vec::with_capacity(cells.len());
    for &entry in &cells {
        let keys = moments.node_keys(grid, entry).unwrap_or_default();
        let t: [Vec<(u32, f64)>; 8] = std::array::from_fn(|a| terms(keys[a]));
        let mut first: Option<usize> = None;
        for &(d, _) in t.iter().flatten() {
            let d = d as usize;
            if d >= parent.len() {
                continue;
            }
            match first {
                None => first = Some(d),
                Some(f) => {
                    let (rf, rd) = (find(&mut parent, f), find(&mut parent, d));
                    parent[rf.max(rd)] = rf.min(rd);
                }
            }
        }
        node_terms.push(t);
    }
    // Clusters by the root of their cells' dofs.
    let mut by_root: BTreeMap<usize, (Vec<u32>, Vec<usize>)> = BTreeMap::new();
    for (i, t) in node_terms.iter().enumerate() {
        let Some(&(d, _)) = t.iter().flatten().next() else {
            continue;
        };
        let r = find(&mut parent, d as usize);
        let slot = by_root.entry(r).or_default();
        slot.0.extend(t.iter().flatten().map(|x| x.0));
        slot.1.push(i);
    }
    by_root
        .into_values()
        .map(|(mut dofs, members)| {
            dofs.sort_unstable();
            dofs.dedup();
            let local = |d: u32| dofs.binary_search(&d).unwrap_or(0);
            let cells = members
                .iter()
                .map(|&i| {
                    let t = &node_terms[i];
                    let nodes =
                        std::array::from_fn(|a| t[a].iter().map(|&(d, w)| (local(d), w)).collect());
                    (cells[i], nodes)
                })
                .collect();
            Cluster { dofs, cells }
        })
        .collect()
}

/// Σ over `cluster`'s cells of Wᵀ M W with M = `per_cell(entry)`, an
/// (8·`block`)² node-block matrix (row 8-node index · block + component)
/// and W the cell's node weights, as a dense (dofs·block)² matrix
/// (row-major). Cells `per_cell` gives `None` for add nothing.
pub fn gather(
    cluster: &Cluster,
    block: usize,
    per_cell: &mut dyn FnMut(usize) -> Option<Vec<f64>>,
) -> Vec<f64> {
    let n = cluster.dofs.len() * block;
    let m = 8 * block;
    let mut out = vec![0.0; n * n];
    for (entry, nodes) in &cluster.cells {
        let Some(k) = per_cell(*entry) else {
            continue;
        };
        if k.len() != m * m {
            continue;
        }
        for a in 0..8 {
            for b in 0..8 {
                for &(da, wa) in &nodes[a] {
                    for &(db, wb) in &nodes[b] {
                        for i in 0..block {
                            for j in 0..block {
                                out[(da * block + i) * n + db * block + j] +=
                                    wa * wb * k[(a * block + i) * m + b * block + j];
                            }
                        }
                    }
                }
            }
        }
    }
    out
}
