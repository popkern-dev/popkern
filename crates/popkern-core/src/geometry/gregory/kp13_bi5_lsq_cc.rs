//! KP13 §3 eq (4) LSQ system — **full Catmull-Clark χ** variant.
//!
//! Companion to [`super::kp13_bi5_lsq`], which minimises the same
//! functional using a bilinear approximation to `χ^k`. This module
//! threads the real Catmull-Clark characteristic map through the
//! chain rule so the composite jet `(g△^k ∘ χ^k)(t)` at each corner
//! `t` is the paper-faithful CC value — including the Jacobian of
//! `χ^k` and its curvature at `t`.
//!
//! # The missing piece over the bilinear module
//!
//! The bilinear χ satisfies `χ_bil^k(u, v) := u·m_k + v·m_{k+1}`,
//! so the χ-domain `(u, v)` **is** the simplex barycentric in the
//! sector wedge `△_k`: no chain-rule conversion is needed and the
//! inner jet has identity Jacobian + vanishing Hessian.
//!
//! Full Catmull-Clark χ is not linear in `(u, v)`: its first
//! partials `∂χ/∂u, ∂χ/∂v` differ from `(m_k, m_{k+1})` (they are
//! `λ m_k, λ m_{k+1}` at the star, for instance) and its second
//! partials are non-zero. Threading it correctly requires:
//!
//! 1. Convert the 2-D-ambient CC-χ jet at the corner to a
//!    **simplex-coord** jet by applying `M⁻¹` component-wise,
//!    where `M = [m_k | m_{k+1}]` is the sector basis matrix.
//! 2. Find the simplex coords of `χ^k(corner)` (for the chain-
//!    rule outer evaluation point).
//! 3. Evaluate `g△^k`'s `OuterHessian` at those simplex coords
//!    via [`outer_hessian_at_natural`].
//! 4. Apply the second-order chain rule in simplex space.
//!
//! Because `(x, y) = u_s · m_k + v_s · m_{k+1}`, the map
//! `(x, y) ↦ (u_s, v_s)` is **linear** (constant Jacobian `M⁻¹`),
//! so transforming the χ-jet's partials by `M⁻¹` moves each of
//! `d_u, d_v, d_uu, d_uv, d_vv` into simplex space without
//! introducing a new curvature term. That's why a single `M⁻¹`
//! multiply per partial suffices.
//!
//! # System layout
//!
//! Identical to [`super::kp13_bi5_lsq::LsqSystem`] — `(12n) × (6n+5)`
//! design matrix, shared across the three output axes. Only the
//! composite-jet evaluator changes.

#![allow(dead_code)]

use super::kp13_bi5_boundary::TensorBorder;
use super::kp13_bi5_chi::{CatmullClarkChi, ChiCorner, ChiJet};
use super::kp13_bi5_free_params::{ApplyError, Bi5FreeParams, apply_free_params};
use super::kp13_bi5_gtri5::{Gtri5, sector_corner};
use super::kp13_bi5_jet::{TaylorJet2, compose_chain_rule, outer_hessian_at_natural};
use super::kp13_bi5_lsq::{
    LsqError, LsqSystem, ParamSlot, RowKind, num_free_params, num_rows, row_index,
};
use crate::geometry::Point3;

/// Inverse of the sector basis matrix `M = [m_k | m_{k+1}]`, as a
/// flat `[[row0], [row1]]` 2×2 matrix.
///
/// Using `m_k = (cos θ_k, sin θ_k)` with `θ_k = 2πk/n`, the
/// determinant simplifies to `sin(2π/n)`, so `M⁻¹` always exists
/// for `n ≥ 3` (singular only in the degenerate valence-2 case,
/// which is out of scope for the KP13 solver).
pub fn inverse_sector_basis(n: usize, k: usize) -> [[f64; 2]; 2] {
    let (mkx, mky) = sector_corner(n, k);
    let (mk1x, mk1y) = sector_corner(n, k + 1);
    let det = mkx * mk1y - mk1x * mky;
    assert!(
        det.abs() > 1e-12,
        "inverse_sector_basis: degenerate M at n={n}, k={k} (det={det})"
    );
    let inv = 1.0 / det;
    [[inv * mk1y, -inv * mk1x], [-inv * mky, inv * mkx]]
}

/// Apply a 2×2 matrix to a 2-vector.
#[inline]
fn mat_vec(m: &[[f64; 2]; 2], v: (f64, f64)) -> (f64, f64) {
    (m[0][0] * v.0 + m[0][1] * v.1, m[1][0] * v.0 + m[1][1] * v.1)
}

/// Transform a 2-D-ambient [`ChiJet`] into a simplex-coord
/// [`ChiJet`] by applying `M⁻¹` to each component.
///
/// Because the ambient-to-simplex map is linear, each partial of
/// the χ-jet transforms independently — no cross-coupling
/// curvature terms arise. The returned jet has `pos`, `d_u`,
/// `d_v`, `d_uu`, `d_uv`, `d_vv` all expressed in simplex coords
/// `(u_s, v_s)`.
pub fn ambient_to_simplex_jet(jet: &ChiJet, n: usize, k: usize) -> ChiJet {
    let m_inv = inverse_sector_basis(n, k);
    ChiJet {
        pos: mat_vec(&m_inv, jet.pos),
        d_u: mat_vec(&m_inv, jet.d_u),
        d_v: mat_vec(&m_inv, jet.d_v),
        d_uu: mat_vec(&m_inv, jet.d_uu),
        d_uv: mat_vec(&m_inv, jet.d_uv),
        d_vv: mat_vec(&m_inv, jet.d_vv),
    }
}

/// Composite jet of `(g△^k ∘ χ^k)(t)` at the given corner of the
/// χ-domain unit square, using the full Catmull-Clark χ.
///
/// Computes:
///
/// 1. `χ_ambient := cc_chi.jet(k, corner)`          [2-D ambient]
/// 2. `χ_simplex := ambient_to_simplex_jet(χ_ambient)` [simplex]
/// 3. `outer := outer_hessian_at_natural(sector_k, χ_simplex.pos)`
/// 4. `compose_chain_rule(outer, χ_simplex)`
pub fn composite_jet_cc(
    gtri5: &Gtri5,
    k: usize,
    corner: ChiCorner,
    cc_chi: &CatmullClarkChi,
) -> TaylorJet2 {
    let n = cc_chi.n;
    let jet_amb = cc_chi.jet(k, corner);
    let jet_simplex = ambient_to_simplex_jet(&jet_amb, n, k);
    let (u_s, v_s) = jet_simplex.pos;
    let outer = outer_hessian_at_natural(gtri5.sector(k), u_s, v_s);
    compose_chain_rule(&outer, &jet_simplex)
}

/// Build the LSQ normal system for the `g△` boundary fit using
/// full Catmull-Clark χ. Same shape as
/// [`super::kp13_bi5_lsq::build_lsq_system`]; only the composite
/// jet evaluator differs.
pub fn build_lsq_system_cc(
    target: &TensorBorder,
    cc_limit: Point3,
) -> Result<LsqSystem, ApplyError> {
    let n = target.n();
    let cc_chi = CatmullClarkChi::new(n);
    let rows = num_rows(n);
    let cols = num_free_params(n);

    let zero_params = Bi5FreeParams::zero(n);
    let g_const = apply_free_params(&zero_params, cc_limit)?;
    let const_jets_mk: Vec<TaylorJet2> = (0..n)
        .map(|k| composite_jet_cc(&g_const, k, ChiCorner::Mk, &cc_chi))
        .collect();
    let const_jets_vk: Vec<TaylorJet2> = (0..n)
        .map(|k| composite_jet_cc(&g_const, k, ChiCorner::Vk, &cc_chi))
        .collect();

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

    let zero_star = Point3::new(0.0, 0.0, 0.0);
    let mut a = vec![vec![0.0; cols]; rows];
    #[allow(clippy::needless_range_loop)]
    for col in 0..cols {
        let slot = ParamSlot::from_flat(col, n);
        let params = unit_params(slot, 0, n);
        let g = apply_free_params(&params, zero_star)?;
        for k in 0..n {
            let j_mk = composite_jet_cc(&g, k, ChiCorner::Mk, &cc_chi);
            let j_vk = composite_jet_cc(&g, k, ChiCorner::Vk, &cc_chi);
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

/// Solve the CC-χ normal equations and inverse-apply to produce
/// the fully-populated `Gtri5`. Paper-faithful counterpart to
/// [`super::kp13_bi5_lsq::solve_g_tri_lsq`].
pub fn solve_g_tri_lsq_cc(target: &TensorBorder, cc_limit: Point3) -> Result<Gtri5, LsqError> {
    let sys = build_lsq_system_cc(target, cc_limit)?;
    let params = super::kp13_bi5_lsq::solve_normal_equations(&sys)?;
    let g = apply_free_params(&params, cc_limit)?;
    Ok(g)
}

/// Construct a `Bi5FreeParams` with one scalar slot set to `1.0`
/// in the given axis. Private duplicate of the helper in
/// [`super::kp13_bi5_lsq`]; kept here to avoid exposing it in
/// the bilinear module's public API.
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

#[cfg(test)]
#[path = "kp13_bi5_lsq_cc_tests.rs"]
mod tests;
