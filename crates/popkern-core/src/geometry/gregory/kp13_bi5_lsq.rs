//! Phase I.11e.3 — Least-squares system for the `g△` boundary fit.
//!
//! Minimises KP13 §3 eq (4):
//!
//! ```text
//!   min Σ_k Σ_{t ∈ {m_k, v_k}} ‖H_t(g△^k ∘ χ^k) − H_t(b^k)‖²
//! ```
//!
//! # Linearity + decoupling
//!
//! [`apply_free_params`] is **linear** in its inputs `(params,
//! cc_limit)` — every step (star pin, Q propagation, per-sector
//! writes, shared-edge sync, circulant, eq-3 back-fills) is an
//! affine map. The composite jet at each unit-square corner is
//! therefore also linear in the free params.
//!
//! Because `g△` is `Point3`-valued and the tri-5 basis is applied
//! component-wise, the three output axes decouple: each is its own
//! `(12n) × (6n+5)` LSQ with a shared design matrix `A` and three
//! distinct right-hand sides.
//!
//! # Bilinear χ scope
//!
//! This module uses [`BilinearChi`]. For the bilinear map the χ-
//! domain `(u, v)` coincides with the sector-wedge natural-simplex
//! `(u, v)` of `g△^k`, so the composite jet at corner `t` is simply
//! `g△^k`'s natural-simplex jet evaluated at the corner's `(u, v)`
//! — no chain-rule conversion. Full CC χ would introduce a
//! Jacobian-inversion step that remains additive work.
//!
//! # System layout (per axis)
//!
//! * Unknowns: `6n + 5` scalars in the order below (see
//!   [`flat_index`]):
//!   * `0..5`    — Q's α, β, γ, δ, ε
//!   * `5..5+n`  — m_k[0..n]
//!   * `5+n..5+2n` — g_30[0..n]
//!   * `5+2n..5+3n` — g_40[0..n]
//!   * `5+3n..5+4n` — g_122[0..n]
//!   * `5+4n..5+5n` — g_023[0..n]
//!   * `5+5n..5+6n` — g_032[0..n]
//! * Residual rows: `6 × 2 × n = 12n` per axis, in the order below
//!   (see [`row_index`]):
//!   * For each sector `k ∈ 0..n`:
//!     * At corner `Mk`: 6 components (value, ∂u, ∂v, ∂²u, ∂²v, ∂²uv)
//!     * At corner `Vk`: 6 components
//!
//! So `A` is `(12n) × (6n+5)` and the same matrix applies to all
//! three axes (only the target vectors differ).

#![allow(dead_code)]

use super::kp13_bi5_boundary::TensorBorder;
use super::kp13_bi5_chi::{BilinearChi, ChiCorner};
use super::kp13_bi5_free_params::{ApplyError, Bi5FreeParams, apply_free_params};
use super::kp13_bi5_gtri5::Gtri5;
use super::kp13_bi5_jet::{TaylorJet2, outer_hessian_at_natural};
use crate::geometry::Point3;

/// Number of free-parameter slots per axis for valence `n`.
#[inline]
pub const fn num_free_params(n: usize) -> usize {
    6 * n + 5
}

/// Number of residual rows per axis.
#[inline]
pub const fn num_rows(n: usize) -> usize {
    12 * n
}

/// Identify a scalar slot in the per-axis flat parameter vector.
///
/// Layout: `[α, β, γ, δ, ε, m_0, ..., m_{n-1}, g30_0, ..., g30_{n-1},
/// g40_0, ..., g40_{n-1}, g122_0, ..., g122_{n-1}, g023_0, ...,
/// g023_{n-1}, g032_0, ..., g032_{n-1}]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParamSlot {
    QAlpha,
    QBeta,
    QGamma,
    QDelta,
    QEpsilon,
    MidEdge(usize),
    G30(usize),
    G40(usize),
    G122(usize),
    G023(usize),
    G032(usize),
}

impl ParamSlot {
    /// Decode a flat index into the corresponding `ParamSlot`.
    pub fn from_flat(flat: usize, n: usize) -> Self {
        if flat < 5 {
            match flat {
                0 => Self::QAlpha,
                1 => Self::QBeta,
                2 => Self::QGamma,
                3 => Self::QDelta,
                4 => Self::QEpsilon,
                _ => unreachable!(),
            }
        } else {
            let off = flat - 5;
            let k = off % n;
            let block = off / n;
            match block {
                0 => Self::MidEdge(k),
                1 => Self::G30(k),
                2 => Self::G40(k),
                3 => Self::G122(k),
                4 => Self::G023(k),
                5 => Self::G032(k),
                _ => panic!(
                    "ParamSlot::from_flat: flat={flat} > num_free_params(n={n})={}",
                    num_free_params(n)
                ),
            }
        }
    }

    /// Convert to flat index.
    pub fn to_flat(self, n: usize) -> usize {
        match self {
            Self::QAlpha => 0,
            Self::QBeta => 1,
            Self::QGamma => 2,
            Self::QDelta => 3,
            Self::QEpsilon => 4,
            Self::MidEdge(k) => 5 + k,
            Self::G30(k) => 5 + n + k,
            Self::G40(k) => 5 + 2 * n + k,
            Self::G122(k) => 5 + 3 * n + k,
            Self::G023(k) => 5 + 4 * n + k,
            Self::G032(k) => 5 + 5 * n + k,
        }
    }
}

/// Identify a scalar residual row in the per-axis flat vector.
/// Layout: for each `k ∈ 0..n`, 6 components at `Mk` followed by 6
/// at `Vk`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowKind {
    Value,
    Du,
    Dv,
    Duu,
    Dvv,
    Duv,
}

impl RowKind {
    pub fn from_offset(off: usize) -> Self {
        match off {
            0 => Self::Value,
            1 => Self::Du,
            2 => Self::Dv,
            3 => Self::Duu,
            4 => Self::Dvv,
            5 => Self::Duv,
            _ => panic!("RowKind offset {off} > 5"),
        }
    }

    pub fn extract_axis(&self, jet: &TaylorJet2, axis: usize) -> f64 {
        let p = match self {
            Self::Value => jet.value,
            Self::Du => jet.d_u,
            Self::Dv => jet.d_v,
            Self::Duu => jet.d_uu,
            Self::Dvv => jet.d_vv,
            Self::Duv => jet.d_uv,
        };
        match axis {
            0 => p.x,
            1 => p.y,
            2 => p.z,
            _ => panic!("axis {axis} out of range"),
        }
    }
}

/// Flat row index for `(sector k, corner c, component kind)`.
#[inline]
pub fn row_index(n: usize, k: usize, corner: ChiCorner, kind: RowKind) -> usize {
    assert!(k < n, "row_index: k={k} >= n={n}");
    let per_sector = 12;
    let corner_off = match corner {
        ChiCorner::Mk => 0,
        ChiCorner::Vk => 6,
        _ => panic!("row_index: LSQ only uses Mk and Vk, got {corner:?}"),
    };
    let kind_off = kind as usize;
    k * per_sector + corner_off + kind_off
}

/// Extract the composite jet of `g△^k ∘ χ_bil^k` at the given corner
/// of the χ-domain unit square.
///
/// Bilinear identification: the χ-domain `(u, v)` equals the wedge
/// natural-simplex `(u, v)`, so we just evaluate `g△^k`'s
/// natural-simplex jet at the corner's `(u, v)`.
///
/// The `outer_hessian_at_natural` helper returns partials in the
/// natural-simplex parameter; the field names `d_x, d_y, ...` are
/// labelled for the physical-coords chain-rule pipeline but under
/// bilinear χ they **are** the composite's `(u, v)` partials.
pub fn composite_jet_bilinear(gtri5: &Gtri5, k: usize, corner: ChiCorner) -> TaylorJet2 {
    let (u, v) = corner.uv();
    let oh = outer_hessian_at_natural(gtri5.sector(k), u, v);
    TaylorJet2 {
        value: oh.value,
        d_u: oh.d_x,
        d_v: oh.d_y,
        d_uu: oh.d_xx,
        d_vv: oh.d_yy,
        d_uv: oh.d_xy,
    }
}

/// The assembled LSQ system. Design matrix `a` has shape
/// `num_rows(n) × num_free_params(n)`; the three axes share it.
/// `b` holds the per-axis target residuals.
#[derive(Debug, Clone)]
pub struct LsqSystem {
    pub n: usize,
    /// `a[row][col]` with `row ∈ 0..12n`, `col ∈ 0..6n+5`.
    pub a: Vec<Vec<f64>>,
    /// `b[axis][row]`, axis ∈ {0, 1, 2} for x, y, z.
    pub b: [Vec<f64>; 3],
}

/// Build the LSQ normal system for the g△ boundary fit.
///
/// `target`: tensor-border `b^k` jets per sector.
/// `cc_limit`: Catmull-Clark limit point for the star (pinned).
///
/// Uses [`BilinearChi`] — χ-domain `(u, v)` ≡ natural-simplex `(u, v)`.
pub fn build_lsq_system(target: &TensorBorder, cc_limit: Point3) -> Result<LsqSystem, ApplyError> {
    let n = target.n();
    let rows = num_rows(n);
    let cols = num_free_params(n);

    // "Constant offset" c: composite jets when params = zero but
    // star is pinned to cc_limit.
    let zero_params = Bi5FreeParams::zero(n);
    let g_const = apply_free_params(&zero_params, cc_limit)?;
    let const_jets_mk: Vec<TaylorJet2> = (0..n)
        .map(|k| composite_jet_bilinear(&g_const, k, ChiCorner::Mk))
        .collect();
    let const_jets_vk: Vec<TaylorJet2> = (0..n)
        .map(|k| composite_jet_bilinear(&g_const, k, ChiCorner::Vk))
        .collect();

    // Target RHS per axis: b[axis][row] = target_jet - const_jet.
    let mut b: [Vec<f64>; 3] = [vec![0.0; rows], vec![0.0; rows], vec![0.0; rows]];
    for (axis, b_axis) in b.iter_mut().enumerate() {
        for k in 0..n {
            let tgt_mk = target.jet(k, ChiCorner::Mk);
            let tgt_vk = target.jet(k, ChiCorner::Vk);
            for kind in [
                RowKind::Value,
                RowKind::Du,
                RowKind::Dv,
                RowKind::Duu,
                RowKind::Dvv,
                RowKind::Duv,
            ] {
                let row_mk = row_index(n, k, ChiCorner::Mk, kind);
                let row_vk = row_index(n, k, ChiCorner::Vk, kind);
                b_axis[row_mk] =
                    kind.extract_axis(&tgt_mk, axis) - kind.extract_axis(&const_jets_mk[k], axis);
                b_axis[row_vk] =
                    kind.extract_axis(&tgt_vk, axis) - kind.extract_axis(&const_jets_vk[k], axis);
            }
        }
    }

    // Design matrix: column j = composite jets with axis-x param j
    // set to 1, cc_limit = 0. Axis choice is arbitrary — same A for
    // all axes by linearity/decoupling.
    let zero_star = Point3::new(0.0, 0.0, 0.0);
    let mut a = vec![vec![0.0; cols]; rows];
    #[allow(clippy::needless_range_loop)]
    for col in 0..cols {
        let slot = ParamSlot::from_flat(col, n);
        let params = unit_params(slot, 0, n);
        let g = apply_free_params(&params, zero_star)?;
        for k in 0..n {
            let j_mk = composite_jet_bilinear(&g, k, ChiCorner::Mk);
            let j_vk = composite_jet_bilinear(&g, k, ChiCorner::Vk);
            for kind in [
                RowKind::Value,
                RowKind::Du,
                RowKind::Dv,
                RowKind::Duu,
                RowKind::Dvv,
                RowKind::Duv,
            ] {
                let row_mk = row_index(n, k, ChiCorner::Mk, kind);
                let row_vk = row_index(n, k, ChiCorner::Vk, kind);
                a[row_mk][col] = kind.extract_axis(&j_mk, 0);
                a[row_vk][col] = kind.extract_axis(&j_vk, 0);
            }
        }
    }

    Ok(LsqSystem { n, a, b })
}

/// Construct a `Bi5FreeParams` with one scalar slot set to `1.0` in
/// the given axis, and all other scalars set to `0.0`.
fn unit_params(slot: ParamSlot, axis: usize, n: usize) -> Bi5FreeParams {
    let mut p = Bi5FreeParams::zero(n);
    let unit = match axis {
        0 => Point3::new(1.0, 0.0, 0.0),
        1 => Point3::new(0.0, 1.0, 0.0),
        2 => Point3::new(0.0, 0.0, 1.0),
        _ => panic!("axis {axis} out of range"),
    };
    match slot {
        ParamSlot::QAlpha => p.q.alpha = unit,
        ParamSlot::QBeta => p.q.beta = unit,
        ParamSlot::QGamma => p.q.gamma = unit,
        ParamSlot::QDelta => p.q.delta = unit,
        ParamSlot::QEpsilon => p.q.epsilon = unit,
        ParamSlot::MidEdge(k) => p.m_k[k] = unit,
        ParamSlot::G30(k) => p.g_30[k] = unit,
        ParamSlot::G40(k) => p.g_40[k] = unit,
        ParamSlot::G122(k) => p.g_122[k] = unit,
        ParamSlot::G023(k) => p.g_023[k] = unit,
        ParamSlot::G032(k) => p.g_032[k] = unit,
    }
    p
}

// Keep BilinearChi referenced so rust-analyzer notes it as intentional.
#[cfg(test)]
#[path = "kp13_bi5_lsq_tests.rs"]
mod tests;

// Suppress the unused import warning on `BilinearChi` in release
// builds — we document it in the module docs but don't need it
// until I.11e.4 wires the CC variant.
#[allow(unused)]
fn _bilinear_chi_ref() -> BilinearChi {
    BilinearChi::new(5)
}
