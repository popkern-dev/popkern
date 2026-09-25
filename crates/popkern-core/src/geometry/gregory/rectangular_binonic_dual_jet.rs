//! Bidegree-(9, 9) Bezier patch with dual cross-jet rational corner
//! cells (PK-519).
//!
//! Identical to [`RectangularGregoryBiNonic`] except that, at each
//! *activated* corner, the four interior cells that both boundary
//! 2-jets read — corner-local `(1,1)` (the twist), `(1,2)`, `(2,1)`,
//! and `(2,2)` — are evaluation-time rational blends of per-boundary
//! claim pairs with **squared** quadrant weights:
//!
//! ```text
//! Q(u,v) = (wv²·P_u + wu²·P_v) / (wu² + wv²)
//! ```
//!
//! where `wu` vanishes on the u-constant boundary through the corner
//! and `wv` on the v-constant one. This is the G2 extension of the
//! Chiyokura-Kimura dual-twist blend
//! (`rectangular_binonic_dual_twist.rs`): with linear weights the
//! off-boundary derivative of the blend is `(P_u − P_v)/w` on the
//! boundary — bounded but nonzero, leaking the other boundary's claim
//! into second derivatives. Squared weights kill that first
//! derivative *on* each boundary, so every activated cell's
//! contribution to a boundary's position, first, and second
//! cross-derivative jets is exactly that boundary's own claim.
//!
//! Why this exists: the KP20 setback fan's corner piece carries two
//! rim sides meeting at the setback corner `c`, on two *different*
//! fillet strips that agree only to second order at `c`. A plain
//! tensor cell near the corner is demanded by side A's depth-2
//! (curvature) row and side B's depth-0/1 (position/tangency) rows at
//! once, and no single value serves both (measured floor ≈ 1 in rim
//! transverse curvature against the 5e-2 gate — see
//! `docs/decisions/2026-09-23-setback-corner-dual-jet-substrate.md`).
//!
//! The twist cell is dual as well: with a shared twist the two rim
//! boundaries' cross-tangent DERIVATIVES along the rim (the
//! geodesic-torsion / mixed-second-fundamental-form term) cannot both
//! be served — measured as a 4.6e-2 diagonal-direction curvature jump
//! on the flagship after the depth-2 cells alone were dualized.

use crate::geometry::gregory::RectangularGregoryBiNonic;
use crate::geometry::gregory::bernstein_binonic::{
    bernstein_binonic_all, bernstein_binonic_deriv_all,
};
use crate::geometry::gregory::bernstein_binonic_second_deriv_all;
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Per-corner dual cross-jet claim block. Each pair is `(P_u, P_v)`:
/// `P_u` is realized on the u-constant boundary through the corner
/// (`u = 0` or `u = 1`), `P_v` on the v-constant one — the dual-twist
/// pair convention.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DualJetCorner {
    /// Corner-local `(1,1)` — the twist cell, read by BOTH
    /// boundaries' d1 (and d2) jets. Dual claims let each boundary
    /// carry its own cross-tangent derivative along the rim (the
    /// geodesic-torsion term the shared twist cannot serve for two
    /// different strips).
    pub cell_11: (Point3, Point3),
    /// Corner-local `(1,2)`: depth 1 from the u-constant boundary
    /// (its tangency band) and depth 2 from the v-constant boundary
    /// (its curvature row). `P_u` therefore carries the u-boundary's
    /// jet-reconciled d1 value; `P_v` is a free d2-row coefficient
    /// for the v-boundary.
    pub cell_12: (Point3, Point3),
    /// Corner-local `(2,1)` — the mirror of [`Self::cell_12`].
    pub cell_21: (Point3, Point3),
    /// Corner-local `(2,2)` — both boundaries' curvature rows cross
    /// here; each claim is free for its own boundary's d2 row.
    pub cell_22: (Point3, Point3),
}

impl DualJetCorner {
    /// Claim block whose pairs both take the given plain cell values
    /// — evaluation is bit-identical to the plain patch until a
    /// claim is written separately.
    pub fn from_plain_cells(
        cell_11: Point3,
        cell_12: Point3,
        cell_21: Point3,
        cell_22: Point3,
    ) -> Self {
        Self {
            cell_11: (cell_11, cell_11),
            cell_12: (cell_12, cell_12),
            cell_21: (cell_21, cell_21),
            cell_22: (cell_22, cell_22),
        }
    }
}

/// Bi-nonic patch with dual cross-jet rational corner cells.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RectangularGregoryBiNonicDualJet {
    /// 10×10 control-point grid, row-major (`grid[i][j]`, `i`
    /// indexes `u`). At an activated corner the four corner-local
    /// `(1,1)`/`(1,2)`/`(2,1)`/`(2,2)` entries are placeholders —
    /// evaluation replaces them with the rational blend from
    /// [`corners`](Self::corners).
    pub grid: [[Point3; 10]; 10],
    /// Per-corner claim blocks, CCW from the parameter origin:
    /// `[(0,0), (0,9), (9,9), (9,0)]` in grid `(i, j)` positions.
    /// `None` = plain corner (cells evaluate from the grid).
    pub corners: [Option<DualJetCorner>; 4],
}

/// Which of a corner block's three cells a grid index addresses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DualJetCell {
    /// Corner-local `(1,1)` — the twist cell.
    Cell11,
    /// Corner-local `(1,2)`.
    Cell12,
    /// Corner-local `(2,1)`.
    Cell21,
    /// Corner-local `(2,2)`.
    Cell22,
}

/// Grid `(i, j)` position of each corner, CCW from the origin.
pub const DUAL_JET_CORNER_POSITIONS: [(usize, usize); 4] = [(0, 0), (0, 9), (9, 9), (9, 0)];

/// Map a grid index to the corner and cell role it addresses when
/// that corner is activated.
#[inline]
pub fn dual_jet_cell_of(i: usize, j: usize) -> Option<(usize, DualJetCell)> {
    match (i, j) {
        (1, 1) => Some((0, DualJetCell::Cell11)),
        (1, 8) => Some((1, DualJetCell::Cell11)),
        (8, 8) => Some((2, DualJetCell::Cell11)),
        (8, 1) => Some((3, DualJetCell::Cell11)),
        (1, 2) => Some((0, DualJetCell::Cell12)),
        (2, 1) => Some((0, DualJetCell::Cell21)),
        (2, 2) => Some((0, DualJetCell::Cell22)),
        (1, 7) => Some((1, DualJetCell::Cell12)),
        (2, 8) => Some((1, DualJetCell::Cell21)),
        (2, 7) => Some((1, DualJetCell::Cell22)),
        (8, 7) => Some((2, DualJetCell::Cell12)),
        (7, 8) => Some((2, DualJetCell::Cell21)),
        (7, 7) => Some((2, DualJetCell::Cell22)),
        (8, 2) => Some((3, DualJetCell::Cell12)),
        (7, 1) => Some((3, DualJetCell::Cell21)),
        (7, 2) => Some((3, DualJetCell::Cell22)),
        _ => None,
    }
}

/// Grid `(i, j)` position of a corner's cell.
#[inline]
pub fn dual_jet_cell_position(corner: usize, cell: DualJetCell) -> (usize, usize) {
    let (ic, jc) = DUAL_JET_CORNER_POSITIONS[corner];
    let (a, b) = match cell {
        DualJetCell::Cell11 => (1, 1),
        DualJetCell::Cell12 => (1, 2),
        DualJetCell::Cell21 => (2, 1),
        DualJetCell::Cell22 => (2, 2),
    };
    let i = if ic == 0 { a } else { 9 - a };
    let j = if jc == 0 { b } else { 9 - b };
    (i, j)
}

/// Quadrant weights for a corner: `wu` vanishes on the u-constant
/// boundary through it, `wv` on the v-constant one. `du`/`dv` are
/// the parameter-derivative signs of `wu`/`wv`.
#[inline]
fn corner_weights(corner: usize, u: f64, v: f64) -> (f64, f64, f64, f64) {
    let (ic, jc) = DUAL_JET_CORNER_POSITIONS[corner];
    let (wu, du) = if ic == 0 { (u, 1.0) } else { (1.0 - u, -1.0) };
    let (wv, dv) = if jc == 0 { (v, 1.0) } else { (1.0 - v, -1.0) };
    (wu, du, wv, dv)
}

#[inline]
fn claim_pair(block: &DualJetCorner, cell: DualJetCell) -> (Point3, Point3) {
    match cell {
        DualJetCell::Cell11 => block.cell_11,
        DualJetCell::Cell12 => block.cell_12,
        DualJetCell::Cell21 => block.cell_21,
        DualJetCell::Cell22 => block.cell_22,
    }
}

/// Blend jet at one cell: value plus first and second parameter
/// derivatives of `Q(u,v) = (wv²·P_u + wu²·P_v) / (wu² + wv²)`.
///
/// With `s = wu²`, `t = wv²`, `D = s + t`, `Δ = P_v − P_u`:
///
/// ```text
/// Q      = P_u + (s/D)·Δ
/// ∂u Q   =  2·wu·du · t·Δ / D²          (zero on both boundaries)
/// ∂v Q   = −2·wv·dv · s·Δ / D²
/// ∂uu Q  =  2·t·Δ·(D − 4s) / D³
/// ∂uv Q  =  4·wu·du·wv·dv·Δ·(s − t) / D³
/// ∂vv Q  = −2·s·Δ·(D − 4t) / D³
/// ```
///
/// Every derivative carries `Δ`, so equal claims short-circuit to
/// the plain cell (bit-identical to the un-promoted patch). At the
/// corner (`D` below the guard) the value takes the pair midpoint
/// and derivatives are zero. Every blended cell's Bernstein factor
/// vanishes at least quadratically there, so the position and first
/// partials are unaffected. The twist `S_uv` reads the `(1,1)` cell
/// through `B₁′(0)² = 81`, however, so at the corner point it takes
/// the midpoint claim. The true twist is direction-dependent there:
/// its limit is `P_u` along the u-constant boundary and `P_v` along
/// the v-constant one (see `boundary_surrogate`).
struct BlendJet {
    q: Point3,
    du: Vec3,
    dv: Vec3,
    duu: Vec3,
    duv: Vec3,
    dvv: Vec3,
}

#[inline]
fn blend_jet(block: &DualJetCorner, cell: DualJetCell, corner: usize, u: f64, v: f64) -> BlendJet {
    let (p_u, p_v) = claim_pair(block, cell);
    let zero = Vec3::new(0.0, 0.0, 0.0);
    let delta = p_v - p_u;
    if delta.length_squared() == 0.0 {
        return BlendJet {
            q: p_u,
            du: zero,
            dv: zero,
            duu: zero,
            duv: zero,
            dvv: zero,
        };
    }
    let (wu, du_sign, wv, dv_sign) = corner_weights(corner, u, v);
    let (s, t) = (wu * wu, wv * wv);
    let d = s + t;
    if d < 1e-30 {
        return BlendJet {
            q: Point3::new(
                0.5 * (p_u.x + p_v.x),
                0.5 * (p_u.y + p_v.y),
                0.5 * (p_u.z + p_v.z),
            ),
            du: zero,
            dv: zero,
            duu: zero,
            duv: zero,
            dvv: zero,
        };
    }
    let d2 = d * d;
    let d3 = d2 * d;
    BlendJet {
        q: Point3::new(
            p_u.x + (s / d) * delta.x,
            p_u.y + (s / d) * delta.y,
            p_u.z + (s / d) * delta.z,
        ),
        du: delta * (2.0 * wu * du_sign * t / d2),
        dv: delta * (-2.0 * wv * dv_sign * s / d2),
        duu: delta * (2.0 * t * (d - 4.0 * s) / d3),
        duv: delta * (4.0 * wu * du_sign * wv * dv_sign * (s - t) / d3),
        dvv: delta * (-2.0 * s * (d - 4.0 * t) / d3),
    }
}

impl RectangularGregoryBiNonicDualJet {
    /// Lift a plain bi-nonic with no activated corners — evaluation
    /// is bit-identical to the plain patch.
    pub fn from_plain(plain: &RectangularGregoryBiNonic) -> Self {
        Self {
            grid: plain.grid,
            corners: [None; 4],
        }
    }

    /// Activate a corner: its claim pairs initialize to the current
    /// grid cell values, so evaluation stays bit-identical until a
    /// claim is written separately.
    pub fn activate_corner(&mut self, corner: usize) {
        let cell = |c: DualJetCell| {
            let (i, j) = dual_jet_cell_position(corner, c);
            self.grid[i][j]
        };
        self.corners[corner] = Some(DualJetCorner::from_plain_cells(
            cell(DualJetCell::Cell11),
            cell(DualJetCell::Cell12),
            cell(DualJetCell::Cell21),
            cell(DualJetCell::Cell22),
        ));
    }

    /// The active claim block addressed by grid index `(i, j)`, when
    /// that cell is a blended placeholder.
    #[inline]
    fn active_cell(&self, i: usize, j: usize) -> Option<(usize, DualJetCell, &DualJetCorner)> {
        let (corner, cell) = dual_jet_cell_of(i, j)?;
        self.corners[corner]
            .as_ref()
            .map(|block| (corner, cell, block))
    }

    /// Evaluate the patch at `(u, v) ∈ [0, 1]²`.
    ///
    /// Accumulation order matches the plain bi-nonic's factored
    /// inner-`v` sum, so an un-activated (or equal-claims) patch
    /// evaluates BIT-identically to [`RectangularGregoryBiNonic`].
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_all(v);
        let mut out = Point3::new(0.0, 0.0, 0.0);
        for (i, &wu) in bu.iter().enumerate() {
            if wu == 0.0 {
                continue;
            }
            let mut inner = Point3::new(0.0, 0.0, 0.0);
            for (j, &w) in bv.iter().enumerate() {
                let p = match self.active_cell(i, j) {
                    Some((corner, cell, block)) => blend_jet(block, cell, corner, u, v).q,
                    None => self.grid[i][j],
                };
                inner.x += p.x * w;
                inner.y += p.y * w;
                inner.z += p.z * w;
            }
            out.x += inner.x * wu;
            out.y += inner.y * wu;
            out.z += inner.z * wu;
        }
        out
    }

    /// `∂S/∂u`: product rule over the tensor basis plus the rational
    /// blend's own `u`-derivative.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_all(u);
        let bu_prime = bernstein_binonic_deriv_all(u);
        let bv = bernstein_binonic_all(v);
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for (i, (&bu_i, &bu_prime_i)) in bu.iter().zip(bu_prime.iter()).enumerate() {
            for (j, &bv_j) in bv.iter().enumerate() {
                let (p, dp_du) = match self.active_cell(i, j) {
                    Some((corner, cell, block)) => {
                        let jet = blend_jet(block, cell, corner, u, v);
                        (jet.q, jet.du)
                    }
                    None => (self.grid[i][j], Vec3::new(0.0, 0.0, 0.0)),
                };
                let w_prime = bu_prime_i * bv_j;
                let w = bu_i * bv_j;
                out = out
                    + Vec3::new(
                        w_prime * p.x + w * dp_du.x,
                        w_prime * p.y + w * dp_du.y,
                        w_prime * p.z + w * dp_du.z,
                    );
            }
        }
        out
    }

    /// `∂S/∂v`: symmetric to [`partial_u`](Self::partial_u).
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_all(v);
        let bv_prime = bernstein_binonic_deriv_all(v);
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for (i, &bu_i) in bu.iter().enumerate() {
            for (j, (&bv_j, &bv_prime_j)) in bv.iter().zip(bv_prime.iter()).enumerate() {
                let (p, dp_dv) = match self.active_cell(i, j) {
                    Some((corner, cell, block)) => {
                        let jet = blend_jet(block, cell, corner, u, v);
                        (jet.q, jet.dv)
                    }
                    None => (self.grid[i][j], Vec3::new(0.0, 0.0, 0.0)),
                };
                let w_prime = bu_i * bv_prime_j;
                let w = bu_i * bv_j;
                out = out
                    + Vec3::new(
                        w_prime * p.x + w * dp_dv.x,
                        w_prime * p.y + w * dp_dv.y,
                        w_prime * p.z + w * dp_dv.z,
                    );
            }
        }
        out
    }

    /// `∂²S/∂u²` — fully analytic: `Σ (Q·Bu″Bv + 2·∂uQ·Bu′Bv +
    /// ∂uuQ·BuBv)`. Analytic seconds keep the boundary 2-jets exact
    /// (an FD scheme would sample the blend's interior transition
    /// zone and read a false curvature near activated corners).
    pub fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        self.second_partial(u, v, SecondKind::Uu)
    }

    /// `∂²S/∂u∂v` — fully analytic.
    pub fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        self.second_partial(u, v, SecondKind::Uv)
    }

    /// `∂²S/∂v²` — fully analytic.
    pub fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        self.second_partial(u, v, SecondKind::Vv)
    }

    fn second_partial(&self, u: f64, v: f64, kind: SecondKind) -> Vec3 {
        let bu = bernstein_binonic_all(u);
        let bu1 = bernstein_binonic_deriv_all(u);
        let bu2 = bernstein_binonic_second_deriv_all(u);
        let bv = bernstein_binonic_all(v);
        let bv1 = bernstein_binonic_deriv_all(v);
        let bv2 = bernstein_binonic_second_deriv_all(v);
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..10 {
            for j in 0..10 {
                // Basis weight of the cell VALUE for this derivative.
                let w2 = match kind {
                    SecondKind::Uu => bu2[i] * bv[j],
                    SecondKind::Uv => bu1[i] * bv1[j],
                    SecondKind::Vv => bu[i] * bv2[j],
                };
                match self.active_cell(i, j) {
                    Some((corner, cell, block)) => {
                        let jet = blend_jet(block, cell, corner, u, v);
                        // Product rule: value·basis″ + cross terms +
                        // blend-second·basis.
                        let (cross, second) = match kind {
                            SecondKind::Uu => {
                                (jet.du * (2.0 * bu1[i] * bv[j]), jet.duu * (bu[i] * bv[j]))
                            }
                            SecondKind::Uv => (
                                jet.du * (bu[i] * bv1[j]) + jet.dv * (bu1[i] * bv[j]),
                                jet.duv * (bu[i] * bv[j]),
                            ),
                            SecondKind::Vv => {
                                (jet.dv * (2.0 * bu[i] * bv1[j]), jet.dvv * (bu[i] * bv[j]))
                            }
                        };
                        out = out
                            + Vec3::new(jet.q.x * w2, jet.q.y * w2, jet.q.z * w2)
                            + cross
                            + second;
                    }
                    None => {
                        let p = self.grid[i][j];
                        out = out + Vec3::new(p.x * w2, p.y * w2, p.z * w2);
                    }
                }
            }
        }
        out
    }

    /// Parameter range, `[0,1] × [0,1]`.
    pub fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, 1.0), (0.0, 1.0))
    }
}

#[derive(Clone, Copy)]
enum SecondKind {
    Uu,
    Uv,
    Vv,
}

#[cfg(test)]
#[path = "rectangular_binonic_dual_jet_tests.rs"]
mod tests;
