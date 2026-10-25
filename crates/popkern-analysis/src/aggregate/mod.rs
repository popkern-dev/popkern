//! Aggregated FEM (AgFEM): the discrete extension that keeps cut-cell
//! systems well conditioned (map §1.3, notes 06).
//!
//! A cell is large when its volume fraction η ≥ γ. Every small cell is
//! assigned to the nearest large cell of the same solid by a
//! multi-source breadth-first search over face neighbours, layer by
//! layer: among the touching neighbours one layer nearer, it takes the
//! root whose material centroid is nearest its own (lowest entry on
//! ties). At a re-entrant corner this keeps a sliver of one leg on that
//! leg's large cell rather than on a cell across the corner that it
//! touches only through the sliver's thin face. A node all
//! of whose active cells are small is ill-posed: its value is the
//! trilinear extrapolation of its root cell's field, so it carries no
//! dof of its own. The reduced system is Eᵀ K E.

pub mod cluster;

use crate::element::poly::{eval, shape};
use crate::grid::Grid;
use crate::model::SolidId;
use crate::moments::SolidMoments;
use crate::moments::fragments::NODE_COPIES;
use crate::report::Warning;

/// Aggregate BFS distance (in cells) beyond which a thin feature is
/// reported.
pub const THIN_AGGREGATE_DISTANCE: usize = 3;

/// Face-neighbour offsets in the fixed search order.
const NEIGHBOURS: [(usize, i64); 6] = [(0, -1), (0, 1), (1, -1), (1, 1), (2, -1), (2, 1)];

/// How one active node maps to dofs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Row {
    /// A free dof.
    Free(u32),
    /// Extrapolated from eight free dofs of its root cell.
    Extended([(u32, f64); 8]),
}

/// The extension operator E from free dofs to active nodes.
#[derive(Debug, Clone, PartialEq)]
pub struct Extension {
    /// Active node keys (grid node · `NODE_COPIES` + copy), ascending.
    pub nodes: Vec<usize>,
    /// The row of each active node, parallel to `nodes`.
    pub rows: Vec<Row>,
    /// Number of free dofs.
    pub free: usize,
}

impl Extension {
    /// Position of node key `node` among the active nodes.
    pub fn position(&self, node: usize) -> Option<usize> {
        self.nodes.binary_search(&node).ok()
    }

    /// E as a sparse matrix: rows are active-node positions, columns
    /// free dofs.
    pub fn to_csr(&self) -> Result<crate::linalg::sparse::Csr, crate::linalg::sparse::SparseError> {
        let mut triplets = Vec::new();
        for (i, row) in self.rows.iter().enumerate() {
            match row {
                Row::Free(d) => triplets.push((i, *d as usize, 1.0)),
                Row::Extended(w) => triplets.extend(w.iter().map(|(d, c)| (i, *d as usize, *c))),
            }
        }
        crate::linalg::sparse::Csr::from_triplets(self.nodes.len(), self.free, &triplets)
    }

    /// Finest-lattice coordinates of each free dof.
    pub fn free_positions(&self, grid: &Grid) -> Vec<[usize; 3]> {
        self.free_keys()
            .map(|k| grid.node_ijk(k / NODE_COPIES).unwrap_or_default())
            .collect()
    }

    /// Dof copy of each free dof (0 except across a slit).
    pub fn free_copies(&self) -> Vec<u8> {
        self.free_keys().map(|k| (k % NODE_COPIES) as u8).collect()
    }

    fn free_keys(&self) -> impl Iterator<Item = usize> + '_ {
        self.nodes
            .iter()
            .zip(&self.rows)
            .filter(|(_, r)| matches!(r, Row::Free(_)))
            .map(|(n, _)| *n)
    }

    /// E x: values at every active node from free values `x`.
    pub fn expand(&self, x: &[f64]) -> Vec<f64> {
        self.rows
            .iter()
            .map(|row| match row {
                Row::Free(d) => x[*d as usize],
                Row::Extended(w) => w.iter().map(|(d, c)| c * x[*d as usize]).sum(),
            })
            .collect()
    }

    /// Eᵀ y: free values from nodal values `y` (parallel to `nodes`),
    /// accumulated in node order.
    pub fn restrict(&self, y: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.free];
        for (row, v) in self.rows.iter().zip(y) {
            match row {
                Row::Free(d) => out[*d as usize] += v,
                Row::Extended(w) => {
                    for (d, c) in w {
                        out[*d as usize] += c * v;
                    }
                }
            }
        }
        out
    }
}

/// The aggregation of one solid's active cells.
#[derive(Debug, Clone, PartialEq)]
pub struct Aggregation {
    /// The solid.
    pub solid: SolidId,
    /// γ actually used (lowered when the solid has no large cell).
    pub gamma: f64,
    /// Root entry of each active entry, parallel to
    /// `SolidMoments::cells`.
    pub root: Vec<usize>,
    /// BFS distance of each active cell from its root.
    pub distance: Vec<usize>,
    /// Number of small (aggregated) cells.
    pub cells_aggregated: usize,
    /// Small cells no face path reached, assigned to the nearest large
    /// cell centre instead.
    pub unreached: usize,
    /// The extension operator.
    pub extension: Extension,
    /// Warnings raised.
    pub warnings: Vec<Warning>,
}

impl Aggregation {
    /// Entry `entry`'s field as its root's: W with value at its node a =
    /// Σ_b W[a][b] · value at the root's node b (`CELL_CORNERS` order),
    /// when every node of `entry` is extended from its root's free dofs
    /// alone (the field there is the root's trilinear). `None` for a
    /// root, or an entry with a free node or one extended from another
    /// root.
    pub fn root_weights(
        &self,
        grid: &Grid,
        moments: &SolidMoments,
        entry: usize,
    ) -> Option<[[f64; 8]; 8]> {
        let r = *self.root.get(entry)?;
        if r == entry {
            return None;
        }
        let e = &self.extension;
        let row = |key: usize| e.position(key).and_then(|p| e.rows.get(p));
        let root_dofs: Vec<u32> = moments
            .node_keys(grid, r)?
            .iter()
            .map(|&k| match row(k) {
                Some(Row::Free(d)) => Some(*d),
                _ => None,
            })
            .collect::<Option<_>>()?;
        let mut w = [[0.0; 8]; 8];
        for (a, &key) in moments.node_keys(grid, entry)?.iter().enumerate() {
            let Some(Row::Extended(weights)) = row(key) else {
                return None;
            };
            for &(d, c) in weights {
                let b = root_dofs.iter().position(|&x| x == d)?;
                w[a][b] += c;
            }
        }
        Some(w)
    }
}

/// Why aggregation failed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AggregationError {
    /// The solid has no active cell.
    #[error("solid {0:?} has no active cell")]
    Empty(SolidId),
    /// A root cell has an ill-posed node (an internal invariant).
    #[error("root cell {0} has an ill-posed node")]
    RootNotFree(usize),
}

fn neighbour(grid: &Grid, cell: usize, axis: usize, step: i64) -> Option<usize> {
    let mut ijk = grid.cell_ijk(cell)?;
    let v = ijk[axis] as i64 + step;
    if v < 0 {
        return None;
    }
    ijk[axis] = v as usize;
    grid.cell_index(ijk)
}

/// Aggregates `moments` with threshold `gamma`. `surface_area` (of the
/// solid) gives the wall estimate 2V/A when no cell is large.
pub fn aggregate(
    grid: &Grid,
    moments: &SolidMoments,
    gamma: f64,
    surface_area: f64,
) -> Result<Aggregation, AggregationError> {
    let cells = &moments.cells;
    if cells.is_empty() {
        return Err(AggregationError::Empty(moments.solid));
    }
    let h = grid.cell_size();
    let eta: Vec<f64> = cells.iter().map(|c| c.volume_fraction(h)).collect();
    let mut warnings = Vec::new();
    let mut gamma_used = gamma;
    if !eta.iter().any(|&e| e >= gamma) {
        gamma_used = eta.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        warnings.push(Warning::UnderResolved {
            solid: moments.solid.0,
            wall: 2.0 * moments.volume / surface_area,
            cell: h,
        });
    }
    let large: Vec<bool> = eta.iter().map(|&e| e >= gamma_used).collect();

    // Multi-source BFS from the large cells, one layer at a time.
    let centroid: Vec<[f64; 3]> = cells
        .iter()
        .map(|c| {
            let lo = grid.cell_min(grid.cell_ijk(c.cell).unwrap_or_default());
            let local = if c.m[0] > 0.0 {
                [c.m[1] / c.m[0], c.m[3] / c.m[0], c.m[9] / c.m[0]]
            } else {
                [0.0; 3]
            };
            [0, 1, 2].map(|a| lo[a] + 0.5 * h * (local[a] + 1.0))
        })
        .collect();
    let gap = |i: usize, j: usize| -> f64 {
        (0..3)
            .map(|a| (centroid[i][a] - centroid[j][a]) * (centroid[i][a] - centroid[j][a]))
            .sum()
    };
    let mut root: Vec<Option<usize>> = vec![None; cells.len()];
    let mut distance = vec![0usize; cells.len()];
    let mut layer: Vec<usize> = Vec::new();
    for i in 0..cells.len() {
        if large[i] {
            root[i] = Some(i);
            layer.push(i);
        }
    }
    let mut depth = 0;
    while !layer.is_empty() {
        depth += 1;
        // Best (gap, root) per newly reached entry.
        let mut best: std::collections::BTreeMap<usize, (f64, usize)> =
            std::collections::BTreeMap::new();
        for &i in &layer {
            let Some(r) = root[i] else {
                continue;
            };
            for (axis, step) in NEIGHBOURS {
                let Some(d) = neighbour(grid, cells[i].cell, axis, step) else {
                    continue;
                };
                for j in moments.entries(d) {
                    if root[j].is_some() || !moments.touch(i, j) {
                        continue;
                    }
                    let g = gap(j, r);
                    let better = best
                        .get(&j)
                        .is_none_or(|&(bg, br)| g < bg || (g == bg && r < br));
                    if better {
                        best.insert(j, (g, r));
                    }
                }
            }
        }
        layer = Vec::with_capacity(best.len());
        for (j, (_, r)) in best {
            root[j] = Some(r);
            distance[j] = depth;
            layer.push(j);
        }
    }
    // Small cells with no face path to a large cell: nearest large cell
    // centre (lowest index on ties).
    let centre = |cell: usize| {
        let ijk = grid.cell_ijk(cell).unwrap_or_default();
        ijk.map(|v| v as f64)
    };
    let mut unreached = 0;
    for i in 0..cells.len() {
        if root[i].is_some() {
            continue;
        }
        unreached += 1;
        let c = centre(cells[i].cell);
        let best = cells
            .iter()
            .enumerate()
            .zip(&large)
            .filter(|(_, l)| **l)
            .map(|((e, cell), _)| {
                let d = centre(cell.cell);
                let r2: f64 = (0..3).map(|a| (d[a] - c[a]) * (d[a] - c[a])).sum();
                (r2, e)
            })
            .fold(None, |acc: Option<(f64, usize)>, x| match acc {
                Some(a) if a.0 <= x.0 => Some(a),
                _ => Some(x),
            });
        root[i] = best.map(|b| b.1);
        distance[i] = usize::MAX;
    }
    let root: Vec<usize> = root
        .into_iter()
        .enumerate()
        .map(|(e, r)| r.unwrap_or(e))
        .collect();
    let max_distance = distance
        .iter()
        .copied()
        .filter(|&d| d != usize::MAX)
        .max()
        .unwrap_or(0);
    if max_distance > THIN_AGGREGATE_DISTANCE {
        warnings.push(Warning::ThinAggregates {
            solid: moments.solid.0,
            max_distance,
        });
    }
    let extension = extension(grid, moments, &large, &root)?;
    Ok(Aggregation {
        solid: moments.solid,
        gamma: gamma_used,
        cells_aggregated: large.iter().filter(|l| !**l).count(),
        root,
        distance,
        unreached,
        extension,
        warnings,
    })
}

fn extension(
    grid: &Grid,
    moments: &SolidMoments,
    large: &[bool],
    root: &[usize],
) -> Result<Extension, AggregationError> {
    let cells = &moments.cells;
    // (node, well-posed, owner cell position): cells in ascending order,
    // so the first cell touching a node is its lowest-index cell.
    let mut touches: Vec<(usize, bool, usize)> = Vec::new();
    for (i, &is_large) in large.iter().enumerate().take(cells.len()) {
        for node in moments.node_keys(grid, i).unwrap_or_default() {
            touches.push((node, is_large, i));
        }
    }
    touches.sort_by_key(|t| (t.0, t.2));
    let mut nodes = Vec::new();
    let mut well_posed = Vec::new();
    let mut owner = Vec::new();
    for t in &touches {
        if nodes.last() == Some(&t.0) {
            if let Some(w) = well_posed.last_mut() {
                *w |= t.1;
            }
        } else {
            nodes.push(t.0);
            well_posed.push(t.1);
            owner.push(t.2);
        }
    }
    let mut free_of = vec![u32::MAX; nodes.len()];
    let mut free = 0u32;
    for (slot, &w) in free_of.iter_mut().zip(&well_posed) {
        if w {
            *slot = free;
            free += 1;
        }
    }
    let h = grid.cell_size();
    let mut rows = Vec::with_capacity(nodes.len());
    for (p, &node) in nodes.iter().enumerate() {
        if well_posed[p] {
            rows.push(Row::Free(free_of[p]));
            continue;
        }
        let r = root[owner[p]];
        let r_cell = cells[r].cell;
        let r_ijk = grid.cell_ijk(r_cell).unwrap_or_default();
        let r_lo = grid.cell_min(r_ijk);
        let x = grid.node_position(grid.node_ijk(node / NODE_COPIES).unwrap_or_default());
        let xi = [0, 1, 2].map(|a| 2.0 * (x[a] - r_lo[a]) / h - 1.0);
        let r_nodes = moments.node_keys(grid, r).unwrap_or_default();
        let mut w = [(0u32, 0.0); 8];
        for (a, slot) in w.iter_mut().enumerate() {
            let q = nodes
                .binary_search(&r_nodes[a])
                .map_err(|_| AggregationError::RootNotFree(r_cell))?;
            if free_of[q] == u32::MAX {
                return Err(AggregationError::RootNotFree(r_cell));
            }
            *slot = (free_of[q], eval(&shape(a), xi));
        }
        rows.push(Row::Extended(w));
    }
    Ok(Extension {
        nodes,
        rows,
        free: free as usize,
    })
}

#[cfg(test)]
#[path = "aggregate_tests.rs"]
mod aggregate_tests;
