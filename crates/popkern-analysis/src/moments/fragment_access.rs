//! Lookups on fragmented moments: entries of a cell, the entry holding
//! a piece or a point, node keys and face contact; a sliver's material.

use std::ops::Range;

use crate::grid::Grid;
use crate::moments::SolidMoments;

use super::{NODE_COPIES, label_at, sliver_at};

impl SolidMoments {
    /// Entries (indices into `cells`) of grid cell `cell`.
    pub fn entries(&self, cell: usize) -> Range<usize> {
        let lo = self.cells.partition_point(|e| e.cell < cell);
        let hi = self.cells.partition_point(|e| e.cell <= cell);
        lo..hi
    }

    /// Whether entry `e` is one fragment of a split cell.
    pub fn is_fragment(&self, e: usize) -> bool {
        self.cells
            .get(e)
            .is_some_and(|c| self.entries(c.cell).len() > 1)
    }

    /// The entry holding clipped piece `index` (in cell `cell`).
    pub fn entry_of_piece(&self, index: usize, cell: usize) -> Option<usize> {
        let f = self
            .fragments
            .piece_fragment
            .get(index)
            .copied()
            .unwrap_or(0);
        self.entries(cell).find(|&e| self.cells[e].fragment == f)
    }

    /// The entry of `cell` holding point `p` (stepped half a lattice
    /// spacing along `inward` first, for a point on the surface).
    pub fn entry_at(
        &self,
        grid: &Grid,
        cell: usize,
        p: [f64; 3],
        inward: Option<[f64; 3]>,
    ) -> Option<usize> {
        let range = self.entries(cell);
        if range.len() <= 1 {
            return (!range.is_empty()).then_some(range.start);
        }
        let split = &self.fragments.split;
        let k = split.binary_search_by_key(&cell, |s| s.cell).ok()?;
        let ijk = grid.cell_ijk(cell)?;
        let depth = self.fragments.depth;
        let s = grid.cell_size() / (1usize << depth) as f64;
        let q = inward.map_or(p, |n| [0, 1, 2].map(|a| p[a] + 0.5 * s * n[a]));
        // A point on a sliver: the step inward can leave material
        // thinner than it, so the point itself is tried first.
        let lattice_count = split[k].labels.iter().copied().max().unwrap_or(0);
        let tag = match sliver_at(grid, &split[k].slivers, p) {
            Some(s) => lattice_count.saturating_add(s as u8 + 1),
            None => label_at(grid, ijk, depth, &split[k].labels, &split[k].slivers, q),
        };
        range
            .clone()
            .find(|&e| self.cells[e].fragment + 1 == tag)
            .or(Some(range.start))
    }

    /// Node keys (node · NODE_COPIES + copy) of entry `e`'s corners in
    /// `CELL_CORNERS` order.
    pub fn node_keys(&self, grid: &Grid, e: usize) -> Option<[usize; 8]> {
        let c = self.cells.get(e)?;
        let nodes = grid.cell_nodes(grid.cell_ijk(c.cell)?)?;
        let copies = self.fragments.corner_copy.get(e).copied().unwrap_or([0; 8]);
        Some(std::array::from_fn(|a| {
            nodes[a] * NODE_COPIES + copies[a] as usize
        }))
    }

    /// Whether face-adjacent entries `a` and `b` share material across
    /// their common face.
    pub fn touch(&self, a: usize, b: usize) -> bool {
        let key = (a.min(b), a.max(b));
        let p = &self.fragments.pair_touch;
        p.binary_search_by_key(&key, |x| x.0)
            .map_or(true, |i| p[i].1)
    }
}

impl super::Sliver {
    /// Whether `p` lies in one of its material boxes, to `margin`.
    pub(crate) fn holds(&self, p: [f64; 3], margin: f64) -> bool {
        self.boxes
            .iter()
            .any(|(lo, hi)| (0..3).all(|a| p[a] >= lo[a] - margin && p[a] <= hi[a] + margin))
    }

    /// Its material boxes with a face on the plane x_axis = `plane` (to
    /// `tol`).
    pub(crate) fn boxes_on(
        &self,
        axis: usize,
        plane: f64,
        tol: f64,
    ) -> impl Iterator<Item = &([f64; 3], [f64; 3])> {
        self.boxes.iter().filter(move |(lo, hi)| {
            (lo[axis] - plane).abs() <= tol || (hi[axis] - plane).abs() <= tol
        })
    }
}
