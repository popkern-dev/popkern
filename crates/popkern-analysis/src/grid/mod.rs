//! The uniform background grid.
//!
//! Cubic cells of size h cover the model's bounding box with one full
//! padding cell on every side, so every boundary cell has a complete
//! neighbourhood. Cells and nodes are numbered lexicographically with x
//! fastest, and every iteration over them follows that order.

use popkern_math::PortableMath;

use crate::options::AnalysisOptions;

mod index;

pub use index::{CELL_CORNERS, CellIter};

/// Largest node count a grid may have, so node indices fit in `u32`.
pub const MAX_NODES: usize = u32::MAX as usize;

/// `n < MAX_NODES`, compared in `u64` so it reads the same on 32-bit
/// targets, where `MAX_NODES` is `usize::MAX`.
fn below_max_nodes(n: usize) -> bool {
    (n as u64) < MAX_NODES as u64
}

/// Growth factor applied to a derived cell size until the budget fits.
const DERIVE_GROWTH: f64 = 1.01;

/// Most cells added to one axis to absorb rounding in extent / h.
const ROUNDING_STEPS: usize = 4;

/// Bound on derive iterations; 1.01^2000 exceeds any rounding excess.
const DERIVE_MAX_STEPS: usize = 2000;

/// A uniform grid of cubic cells.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Grid {
    origin: [f64; 3],
    cell_size: f64,
    dims: [usize; 3],
}

/// Why a grid cannot be built.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum GridError {
    /// A bounding-box coordinate is NaN or infinite.
    #[error("bounding box is not finite")]
    NonFiniteBounds,
    /// The bounding box has zero or negative extent on an axis.
    #[error("bounding box has no extent on axis {axis}")]
    ZeroExtent {
        /// Axis index (0 = x).
        axis: usize,
    },
    /// The cell size is not a positive finite number.
    #[error("cell size {0} is not positive and finite")]
    BadCellSize(f64),
    /// The grid would have more nodes than [`MAX_NODES`].
    #[error("grid of {dims:?} cells exceeds the node index range")]
    TooLarge {
        /// Requested cells per axis (saturated).
        dims: [usize; 3],
    },
    /// No cell size within the search keeps the cells within budget.
    #[error("no cell size fits {budget} cells")]
    BudgetUnreachable {
        /// The cell budget.
        budget: usize,
    },
    /// The cell budget is zero.
    #[error("cell budget is zero")]
    NoBudget,
}

fn check_bounds(min: [f64; 3], max: [f64; 3]) -> Result<[f64; 3], GridError> {
    if min.iter().chain(&max).any(|c| !c.is_finite()) {
        return Err(GridError::NonFiniteBounds);
    }
    let mut extent = [0.0; 3];
    for axis in 0..3 {
        extent[axis] = max[axis] - min[axis];
        if extent[axis] <= 0.0 {
            return Err(GridError::ZeroExtent { axis });
        }
    }
    Ok(extent)
}

/// Unpadded cells per axis needed to cover `extent` with cells of `h`.
fn cover(extent: [f64; 3], h: f64) -> [usize; 3] {
    // `as` saturates for huge values; TooLarge catches that later.
    extent.map(|e| ((e / h).ceil() as usize).max(1))
}

impl Grid {
    /// The grid of cells of size `cell_size` covering the box
    /// `[min, max]` with one padding cell on every side.
    pub fn covering(min: [f64; 3], max: [f64; 3], cell_size: f64) -> Result<Grid, GridError> {
        let extent = check_bounds(min, max)?;
        if !(cell_size.is_finite() && cell_size > 0.0) {
            return Err(GridError::BadCellSize(cell_size));
        }
        let covered = cover(extent, cell_size);
        let too_large = GridError::TooLarge {
            dims: covered.map(|n| n.saturating_add(2)),
        };
        let mut dims = [0; 3];
        let mut origin = [0.0; 3];
        for axis in 0..3 {
            let mut n = covered[axis];
            if !below_max_nodes(n) {
                return Err(too_large);
            }
            origin[axis] = min[axis] - cell_size;
            // Rounding in e / h can leave the last cell short of max;
            // add cells until the far padding cell is whole. Below
            // MAX_NODES cells one or two steps suffice.
            let mut steps = 0;
            while origin[axis] + cell_size * (n as f64 + 1.0) < max[axis] {
                n += 1;
                steps += 1;
                if steps > ROUNDING_STEPS {
                    return Err(too_large);
                }
            }
            dims[axis] = n + 2;
        }
        let nodes = dims
            .iter()
            .try_fold(1usize, |acc, &d| acc.checked_mul(d.checked_add(1)?));
        match nodes {
            Some(n) if (n as u64) <= MAX_NODES as u64 => Ok(Grid {
                origin,
                cell_size,
                dims,
            }),
            _ => Err(GridError::TooLarge { dims }),
        }
    }

    /// The grid for a model with bounding box `[min, max]` under
    /// `options`: cells of `target_cell_size` when set, otherwise the
    /// smallest cell size (searched in steps of 1%) whose unpadded
    /// cells covering the box number at most `max_active_cells`. Active
    /// cells are a subset of those, so the budget holds before any
    /// rasterising.
    pub fn for_model(
        min: [f64; 3],
        max: [f64; 3],
        options: &AnalysisOptions,
    ) -> Result<Grid, GridError> {
        let extent = check_bounds(min, max)?;
        // The grid origin moves by `grid_shift` independently of the
        // geometry; without it, translating a part moves the grid with
        // it and the cut pattern never changes.
        let min = [0, 1, 2].map(|a| min[a] - options.grid_shift[a]);
        if let Some(h) = options.target_cell_size {
            return Grid::covering(min, max, h);
        }
        let budget = options.max_active_cells;
        if budget == 0 {
            return Err(GridError::NoBudget);
        }
        let volume = extent[0] * extent[1] * extent[2];
        let mut h = (volume / budget as f64).pcbrt();
        for _ in 0..DERIVE_MAX_STEPS {
            // Too large is not final here: a bigger h may still fit.
            if let Ok(grid) = Grid::covering(min, max, h)
                && grid.interior_cell_count() <= budget
            {
                return Ok(grid);
            }
            h *= DERIVE_GROWTH;
        }
        Err(GridError::BudgetUnreachable { budget })
    }

    /// Position of node (0, 0, 0), the grid's minimum corner.
    pub fn origin(&self) -> [f64; 3] {
        self.origin
    }

    /// Edge length of every cell.
    pub fn cell_size(&self) -> f64 {
        self.cell_size
    }

    /// Cells per axis, padding included.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Nodes per axis (cells + 1).
    pub fn node_dims(&self) -> [usize; 3] {
        self.dims.map(|d| d + 1)
    }

    /// Total number of cells.
    pub fn cell_count(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Number of cells inside the padding layer.
    pub fn interior_cell_count(&self) -> usize {
        self.dims.iter().map(|d| d - 2).product()
    }

    /// Total number of nodes.
    pub fn node_count(&self) -> usize {
        let [a, b, c] = self.node_dims();
        a * b * c
    }

    /// Minimum corner of cell `ijk`.
    pub fn cell_min(&self, ijk: [usize; 3]) -> [f64; 3] {
        let h = self.cell_size;
        [0, 1, 2].map(|a| self.origin[a] + h * ijk[a] as f64)
    }

    /// Position of node `ijk`.
    pub fn node_position(&self, ijk: [usize; 3]) -> [f64; 3] {
        self.cell_min(ijk)
    }

    /// The cell containing point `p`, if `p` is inside the grid. A
    /// point on an interior cell boundary belongs to the upper cell.
    pub fn cell_of_point(&self, p: [f64; 3]) -> Option<[usize; 3]> {
        let mut ijk = [0; 3];
        for a in 0..3 {
            let t = (p[a] - self.origin[a]) / self.cell_size;
            if t.is_nan() || t < 0.0 || t > self.dims[a] as f64 {
                return None;
            }
            ijk[a] = (t.floor() as usize).min(self.dims[a] - 1);
        }
        Some(ijk)
    }
}

#[cfg(test)]
#[path = "grid_tests.rs"]
mod grid_tests;
