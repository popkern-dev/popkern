//! Normal extrapolation at flat surfaces (PH-069).
//!
//! Trilinear interpolation of nodal σ* at a surface point leans on the
//! boundary nodes, whose values are extrapolated from patches that move
//! with the grid: a wall's surface stress varies 8% over grid shifts at
//! 3 cells. Extrapolating along the inward normal from σ* one and two
//! cells inside, σ(x) = 2σ*(x + h n) − σ*(x + 2h n), reads the field
//! where the patches are interior: 3.3% (PH-036's measurement). It
//! loses peaks where the surface is concave within that reach (a notch,
//! bore or fillet, a re-entrant corner), so it applies only where the
//! point lies behind or on the plane of every surface piece within two
//! cells: a flat wall, its opposite face, a convex edge. A tessellated
//! curve (more distinct planes than a cap) keeps interpolation too.
//! Raisers keep interpolation, and their peaks are re-solved by
//! submodels.

use std::collections::BTreeMap;

use crate::clip::ClippedSolid;
use crate::grid::Grid;

/// Surface planes within this many cells (Chebyshev) of a point's
/// cell are checked.
const REACH_CELLS: usize = 2;

/// More distinct planes than this within reach is a curved surface
/// (tessellated), left to interpolation.
const MAX_PLANES: usize = 24;

/// A surface plane: unit outward normal n and offset d = n·x on it.
type Plane = ([f64; 3], f64);

/// The surface planes within [`REACH_CELLS`] of each cell holding surface pieces,
/// `None` past [`MAX_PLANES`].
pub(super) struct FlatCells {
    planes: BTreeMap<usize, Option<Vec<Plane>>>,
    tol: f64,
}

impl FlatCells {
    /// Collects the planes around every cell holding pieces of
    /// `clipped`.
    pub fn new(grid: &Grid, clipped: &ClippedSolid) -> FlatCells {
        let h = grid.cell_size();
        let key = |x: f64| (x / (1.0e-9 * h)).round() as i64;
        // Distinct planes by cell, in piece order.
        let mut by_cell: BTreeMap<usize, Vec<Plane>> = BTreeMap::new();
        for p in &clipped.pieces {
            let v = p.polygon.vertices();
            let Some(a) = v.first() else {
                continue;
            };
            let n = p.normal;
            let plane = (n, n[0] * a[0] + n[1] * a[1] + n[2] * a[2]);
            let planes = by_cell.entry(p.cell).or_default();
            let same = |q: &Plane| {
                (0..3).all(|k| (q.0[k] * 1e9).round() == (n[k] * 1e9).round())
                    && key(q.1) == key(plane.1)
            };
            if !planes.iter().any(same) {
                planes.push(plane);
            }
        }
        let dims = grid.dims();
        let r = REACH_CELLS as i64;
        let mut planes = BTreeMap::new();
        for &cell in by_cell.keys() {
            let Some(ijk) = grid.cell_ijk(cell) else {
                continue;
            };
            let mut near: Vec<Plane> = Vec::new();
            let mut over = false;
            'search: for dz in -r..=r {
                for dy in -r..=r {
                    for dx in -r..=r {
                        let q = [ijk[0] as i64 + dx, ijk[1] as i64 + dy, ijk[2] as i64 + dz];
                        if (0..3).any(|a| q[a] < 0 || q[a] >= dims[a] as i64) {
                            continue;
                        }
                        let Some(c) = grid.cell_index(q.map(|v| v as usize)) else {
                            continue;
                        };
                        for plane in by_cell.get(&c).into_iter().flatten() {
                            let same = |x: &Plane| {
                                (0..3).all(|k| (x.0[k] * 1e9).round() == (plane.0[k] * 1e9).round())
                                    && key(x.1) == key(plane.1)
                            };
                            if !near.iter().any(same) {
                                near.push(*plane);
                                if near.len() > MAX_PLANES {
                                    over = true;
                                    break 'search;
                                }
                            }
                        }
                    }
                }
            }
            planes.insert(cell, (!over).then_some(near));
        }
        FlatCells {
            planes,
            tol: 1.0e-6 * h,
        }
    }

    /// Whether surface point `x` in `cell` lies behind or on every
    /// surface plane within reach: the surface around it is flat or
    /// convex (a wall, its opposite face, a convex edge), with no
    /// re-entrant corner or concave curve to lose a peak at.
    pub fn is_flat(&self, cell: usize, x: [f64; 3]) -> bool {
        matches!(self.planes.get(&cell), Some(Some(planes))
            if planes.iter().all(|(n, d)| n[0] * x[0] + n[1] * x[1] + n[2] * x[2] - d <= self.tol))
    }
}

#[cfg(test)]
#[path = "normal_tests.rs"]
mod tests;
