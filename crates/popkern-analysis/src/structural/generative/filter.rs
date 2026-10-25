//! The density filter (Bruns & Tortorelli; top88 eqn 9) over a solve's
//! entries: x̃_e = Σ_j w_ej x_j / S_e with w_ej = H_ej v_j,
//! H_ej = max(0, r − |c_e − c_j|) and S_e = Σ_j w_ej, centres c and
//! volumes v from the moments. Entries of different solids never mix.
//! The sensitivities go back by the chain rule (top88 eqn 10), the
//! transpose of the weights.

use crate::grid::Grid;
use crate::moments::SolidMoments;
use crate::recover::cell::centroid;

/// One entry of the design: its solid, its index in that solid, its
/// centre (mm) and volume (mm³).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesignEntry {
    /// Solid index.
    pub solid: usize,
    /// Entry index within the solid's moments.
    pub entry: usize,
    /// Centre, mm.
    pub centre: [f64; 3],
    /// Volume, mm³.
    pub volume: f64,
    /// The lower corner of its grid cell, mm.
    pub corner: [f64; 3],
}

/// Every active entry of every solid, in (solid, entry) order.
pub fn design_entries(grid: &Grid, solids: &[&SolidMoments]) -> Vec<DesignEntry> {
    let h = grid.cell_size();
    let mut out = Vec::new();
    for (s, moments) in solids.iter().enumerate() {
        for (e, cell) in moments.cells.iter().enumerate() {
            let Some(ijk) = grid.cell_ijk(cell.cell) else {
                continue;
            };
            let lo = grid.cell_min(ijk);
            let xi = centroid(&cell.m);
            out.push(DesignEntry {
                solid: s,
                entry: e,
                centre: [0, 1, 2].map(|a| lo[a] + 0.5 * h * (xi[a] + 1.0)),
                volume: cell.m[0].max(0.0),
                corner: lo,
            });
        }
    }
    out
}

/// The filter's weights, one row of (neighbour, w) per entry.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityFilter {
    starts: Vec<usize>,
    neighbours: Vec<(u32, f64)>,
    sums: Vec<f64>,
}

impl DensityFilter {
    /// The filter of radius `radius` (mm) over `entries` (from
    /// [`design_entries`] on the same solids). A radius below the
    /// spacing between centres leaves every entry its own.
    pub fn new(
        grid: &Grid,
        solids: &[&SolidMoments],
        entries: &[DesignEntry],
        radius: f64,
    ) -> DensityFilter {
        let h = grid.cell_size();
        let reach = (radius / h).ceil().max(0.0) as i64;
        let dims = grid.dims();
        // Global index of each solid's first entry.
        let mut first = Vec::with_capacity(solids.len());
        let mut count = 0;
        for m in solids {
            first.push(count);
            count += m.cells.len();
        }
        let mut starts = Vec::with_capacity(entries.len() + 1);
        let mut neighbours = Vec::new();
        let mut sums = Vec::with_capacity(entries.len());
        for (me, d) in entries.iter().enumerate() {
            starts.push(neighbours.len());
            let moments = solids[d.solid];
            let ijk = moments.cells[d.entry].cell;
            let Some(ijk) = grid.cell_ijk(ijk) else {
                sums.push(0.0);
                continue;
            };
            let mut sum = 0.0;
            for dz in -reach..=reach {
                for dy in -reach..=reach {
                    for dx in -reach..=reach {
                        let q = [ijk[0] as i64 + dx, ijk[1] as i64 + dy, ijk[2] as i64 + dz];
                        if (0..3).any(|a| q[a] < 0 || q[a] >= dims[a] as i64) {
                            continue;
                        }
                        let Some(c) = grid.cell_index(q.map(|v| v as usize)) else {
                            continue;
                        };
                        for j in moments.entries(c) {
                            let g = first[d.solid] + j;
                            let other = &entries[g];
                            let dist = (0..3)
                                .map(|a| {
                                    let t = d.centre[a] - other.centre[a];
                                    t * t
                                })
                                .sum::<f64>()
                                .sqrt();
                            // An entry always weighs itself, even when the
                            // radius is below its own size.
                            let hw = if g == me {
                                (radius - dist).max(f64::MIN_POSITIVE)
                            } else {
                                radius - dist
                            };
                            if hw > 0.0 {
                                let w = hw * other.volume;
                                neighbours.push((g as u32, w));
                                sum += w;
                            }
                        }
                    }
                }
            }
            sums.push(sum);
        }
        starts.push(neighbours.len());
        DensityFilter {
            starts,
            neighbours,
            sums,
        }
    }

    fn row(&self, e: usize) -> &[(u32, f64)] {
        &self.neighbours[self.starts[e]..self.starts[e + 1]]
    }

    /// x̃ = filter(x).
    pub fn apply(&self, x: &[f64]) -> Vec<f64> {
        (0..self.sums.len())
            .map(|e| {
                let s = self.sums[e];
                if s <= 0.0 {
                    return x[e];
                }
                self.row(e)
                    .iter()
                    .map(|&(j, w)| w * x[j as usize])
                    .sum::<f64>()
                    / s
            })
            .collect()
    }

    /// ∂ψ/∂x from ∂ψ/∂x̃ (the chain rule through [`Self::apply`]): the
    /// transpose, each row e scattering g_e w_ej / S_e to its
    /// neighbours j in row order. (Reading H back from row j as
    /// w_je / v_e lost every entry of zero volume, the cells that only
    /// touch the part, whose sensitivity the print filter makes
    /// non-zero.)
    pub fn back(&self, g: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.sums.len()];
        for (e, &s) in self.sums.iter().enumerate() {
            if s <= 0.0 {
                out[e] += g[e];
                continue;
            }
            let ge = g[e] / s;
            for &(j, w) in self.row(e) {
                out[j as usize] += w * ge;
            }
        }
        out
    }

    /// Stored neighbour weights (for the memory report).
    pub fn weights(&self) -> usize {
        self.neighbours.len()
    }
}
