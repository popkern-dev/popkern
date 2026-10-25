//! Lexicographic cell and node numbering (x fastest, then y, then z).

use super::Grid;

/// Offsets of a cell's eight corner nodes, in the fixed local order
/// used by every element: x fastest, then y, then z.
pub const CELL_CORNERS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

fn linear(dims: [usize; 3], ijk: [usize; 3]) -> Option<usize> {
    if (0..3).any(|a| ijk[a] >= dims[a]) {
        return None;
    }
    Some(ijk[0] + dims[0] * (ijk[1] + dims[1] * ijk[2]))
}

fn unlinear(dims: [usize; 3], index: usize) -> Option<[usize; 3]> {
    if index >= dims[0] * dims[1] * dims[2] {
        return None;
    }
    let i = index % dims[0];
    let rest = index / dims[0];
    Some([i, rest % dims[1], rest / dims[1]])
}

impl Grid {
    /// Linear index of cell `ijk`, or `None` outside the grid.
    pub fn cell_index(&self, ijk: [usize; 3]) -> Option<usize> {
        linear(self.dims(), ijk)
    }

    /// Cell coordinates of linear index `index`, or `None` past the end.
    pub fn cell_ijk(&self, index: usize) -> Option<[usize; 3]> {
        unlinear(self.dims(), index)
    }

    /// Linear index of node `ijk`, or `None` outside the grid.
    pub fn node_index(&self, ijk: [usize; 3]) -> Option<usize> {
        linear(self.node_dims(), ijk)
    }

    /// Node coordinates of linear index `index`, or `None` past the end.
    pub fn node_ijk(&self, index: usize) -> Option<[usize; 3]> {
        unlinear(self.node_dims(), index)
    }

    /// Linear node indices of cell `ijk`'s corners in [`CELL_CORNERS`]
    /// order, or `None` outside the grid.
    pub fn cell_nodes(&self, ijk: [usize; 3]) -> Option<[usize; 8]> {
        self.cell_index(ijk)?;
        let nd = self.node_dims();
        let base = ijk[0] + nd[0] * (ijk[1] + nd[1] * ijk[2]);
        Some(CELL_CORNERS.map(|[dx, dy, dz]| base + dx + nd[0] * (dy + nd[1] * dz)))
    }

    /// Every cell, in linear-index order.
    pub fn cells(&self) -> CellIter {
        CellIter {
            dims: self.dims(),
            next: 0,
            count: self.cell_count(),
        }
    }
}

/// Iterator over a grid's cells in lexicographic order (x fastest).
#[derive(Debug, Clone)]
pub struct CellIter {
    dims: [usize; 3],
    next: usize,
    count: usize,
}

impl Iterator for CellIter {
    type Item = [usize; 3];

    fn next(&mut self) -> Option<[usize; 3]> {
        if self.next >= self.count {
            return None;
        }
        let ijk = unlinear(self.dims, self.next);
        self.next += 1;
        ijk
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let left = self.count - self.next;
        (left, Some(left))
    }
}

impl ExactSizeIterator for CellIter {}
