//! The polynomial surrogate of a dual cross-jet bi-nonic along one
//! patch boundary (PK-541).
//!
//! On a boundary through an activated corner, the squared-weight blend
//! `Q = (wv²·P_u + wu²·P_v) / (wu² + wv²)` is exactly that boundary's
//! claim through second derivatives. For example, on `u = 0`
//! (`wu = 0`):
//! - `Q = P_u`;
//! - `∂u Q`, `∂v Q`, `∂uv Q` and `∂vv Q` all carry a factor `wu`, so
//!   they vanish;
//! - `∂uu Q = 2Δ/wv²` survives, but only multiplied by a basis factor
//!   `B_i(u)` with `i ≥ 1`, which vanishes on `u = 0`.
//!
//! A corner not on the boundary contributes through basis factors that
//! vanish to order ≥ 7 there, with its blend smooth away from its
//! corner, so its claims do not reach the boundary 2-jet.
//!
//! The patch's position, first and second partials along the boundary
//! are therefore those of the plain bi-nonic whose blended cells take
//! the boundary's claims. That holds at every point of the open
//! boundary, and for the limits at its ends taken along the boundary.
//!
//! At an activated corner point itself, the twist `S_uv` is
//! direction-dependent. The cell `(1,1)` enters it through
//! `B₁′(0)² = 81`, and `Q₁₁` tends to `P_u` along one boundary and to
//! `P_v` along the other. The evaluator's corner convention (the pair
//! midpoint) is neither limit. A seam's quantity at its end is the
//! limit along the seam, which is the surrogate's value.
//!
//! The rational exact form cannot give these enclosures: its weight
//! `Π (wu² + wv²)` vanishes at every activated corner, and its
//! derivative bounds blow up approaching one.

use crate::geometry::gregory::RectangularGregoryBiNonic;
use crate::geometry::gregory::rectangular_binonic_dual_jet::{
    DUAL_JET_CORNER_POSITIONS, DualJetCell, RectangularGregoryBiNonicDualJet,
    dual_jet_cell_position,
};

/// One of the four boundaries of a `[0, 1]²` patch chart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PatchBoundary {
    /// The `u = 0` boundary.
    U0,
    /// The `u = 1` boundary.
    U1,
    /// The `v = 0` boundary.
    V0,
    /// The `v = 1` boundary.
    V1,
}

impl PatchBoundary {
    /// The boundary a chart line `u = value` (`along_v` true) or
    /// `v = value` lies on, when `value` is exactly 0 or 1.
    pub fn of_chart_line(along_v: bool, value: f64) -> Option<Self> {
        let side = if value == 0.0 {
            0
        } else if value == 1.0 {
            1
        } else {
            return None;
        };
        Some(match (along_v, side) {
            (true, 0) => Self::U0,
            (true, _) => Self::U1,
            (false, 0) => Self::V0,
            (false, _) => Self::V1,
        })
    }
}

impl RectangularGregoryBiNonicDualJet {
    /// The plain bi-nonic whose position, first and second partials
    /// equal this patch's everywhere on `boundary` (see the module
    /// docs). Each activated corner on the boundary contributes the
    /// claims of that boundary. Corners off it contribute their
    /// u-boundary claims, which the boundary 2-jet does not read.
    pub fn boundary_surrogate(&self, boundary: PatchBoundary) -> RectangularGregoryBiNonic {
        let mut grid = self.grid;
        for (corner, block) in self.corners.iter().enumerate() {
            let Some(block) = block else {
                continue;
            };
            let (_, jc) = DUAL_JET_CORNER_POSITIONS[corner];
            let v_boundary_claim = match boundary {
                PatchBoundary::V0 => jc == 0,
                PatchBoundary::V1 => jc == 9,
                PatchBoundary::U0 | PatchBoundary::U1 => false,
            };
            for (cell, (p_u, p_v)) in [
                (DualJetCell::Cell11, block.cell_11),
                (DualJetCell::Cell12, block.cell_12),
                (DualJetCell::Cell21, block.cell_21),
                (DualJetCell::Cell22, block.cell_22),
            ] {
                let (i, j) = dual_jet_cell_position(corner, cell);
                grid[i][j] = if v_boundary_claim { p_v } else { p_u };
            }
        }
        RectangularGregoryBiNonic { grid }
    }
}

#[cfg(test)]
#[path = "rectangular_binonic_dual_jet_boundary_tests.rs"]
mod tests;
