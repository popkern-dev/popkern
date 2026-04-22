//! Phase I.11.d.2 — Characteristic map `χ` for the KP13 bi-5 LSQ.
//!
//! KP13 §3 eq (4) minimises `Σ_k Σ_{t ∈ {m_k, v_k}} ‖H_t(g△^k ∘ χ^k)
//! − H_t(b^k)‖²` over free BB-coefs of `g△`. Evaluating the
//! composite `g△^k ∘ χ^k` at the corners of the unit square requires
//! the characteristic map `χ^k : [0, 1]² → Ω`.
//!
//! # Paper specification (KP13 §3)
//!
//! > `χ : [0..1]^{2×n} → R², normalized to edge-length 1 [...] so
//! > that `g△^k ∘ χ^k(1, 1) = v_k` (`χ^k(0, 0) = (0, 0)`, Fig. 5a)
//! > and `g△^k_{050} = m_k` (Fig. 5b).
//!
//! Standard reading: `χ` is the planar Catmull-Clark characteristic
//! map restricted to the n extraordinary-vertex sectors, normalised
//! so the first ring has edge length 1. Unit-square corners map to:
//!
//! | corner `(u, v)` | `χ^k(u, v)` in Ω |
//! |-----------------|-------------------|
//! | `(0, 0)`        | `o = (0, 0)` (star) |
//! | `(1, 0)`        | `m_k` (= `(c_k, s_k)`) |
//! | `(0, 1)`        | `m_{k+1}` (= `(c_{k+1}, s_{k+1})`) |
//! | `(1, 1)`        | `v_k` (outside the natural n-gon Ω) |
//!
//! # Current status — bilinear approximation
//!
//! Two variants ship side-by-side:
//!
//! **[`BilinearChi`]** — `χ^k(u, v) := u·m_k + v·m_{k+1}`. Satisfies
//! the first three corner conditions exactly; interprets `v_k :=
//! m_k + m_{k+1}` for the last corner. All partial derivatives are
//! constant (bilinear has no curvature). Fast, simple, and a useful
//! baseline for the LSQ machinery.
//!
//! **[`CatmullClarkChi`]** — a CC-consistent upgrade that matches
//! [`BilinearChi`] at the four corner positions but introduces the
//! **subdominant-eigenvalue scaling** of partials at the star. The
//! construction is a bi-cubic Bézier patch whose CPs are bilinear
//! everywhere except near the star (`P_{1,0}` and `P_{0,1}`), where
//! they are scaled by `λ_n / 3`. The resulting partials at the star
//! are:
//!
//! ```text
//!   ∂_u χ^k(0, 0)   = λ_n · m_k
//!   ∂_v χ^k(0, 0)   = λ_n · m_{k+1}
//!   ∂²_u χ^k(0, 0)  = 4(1 − λ_n) · m_k
//!   ∂²_v χ^k(0, 0)  = 4(1 − λ_n) · m_{k+1}
//!   ∂²_uv χ^k(0, 0) = 3(1 − λ_n) · (m_k + m_{k+1})
//! ```
//!
//! where `λ_n := (c + 5 + √((c + 1)(c + 9))) / 16`, `c := cos(2π/n)`
//! — the Catmull-Clark subdominant eigenvalue (see polar-kp-math.md
//! §3.6 and the standard CC literature).
//!
//! At the three non-star corners, partials match [`BilinearChi`] by
//! construction. This captures the dominant CC behaviour at the star
//! (where it matters most for the LSQ weighting) without requiring
//! the full Stam 1998 arbitrary-parameter evaluator.
//!
//! # Fidelity caveat (CatmullClarkChi)
//!
//! This is **not** a full reproduction of the CC characteristic map.
//! Specifically:
//!
//! * `v_k` is still `m_k + m_{k+1}`; the true CC `v_k` differs by a
//!   sub-subdominant-eigenvalue-dependent correction.
//! * Partials at `(1, 0)`, `(0, 1)`, `(1, 1)` are bilinear, not the
//!   smooth bi-cubic values a full CC surface would give.
//!
//! These remain additive upgrades tractable in a later session.
//!
//! See `docs/journal/polar-kp-math.md` §1.5 (CC limit point) and
//! §3.6 (λ formula).

#![allow(dead_code)]

use super::kp13_bi5_gtri5::sector_corner;
use std::f64::consts::PI;

/// Catmull-Clark subdominant eigenvalue at valence `n`:
///
/// ```text
///   λ_n := (c + 5 + √((c + 1)(c + 9))) / 16,   c := cos(2π/n)
/// ```
///
/// For `n = 4` (regular quad vertex) `λ = 1/2` — the standard
/// bicubic B-spline subdivision factor. For `n → ∞`, `λ → 1`.
///
/// Panics if `n < 3`.
pub fn cc_subdominant_eigenvalue(n: usize) -> f64 {
    assert!(n >= 3, "CC subdominant eigenvalue needs n ≥ 3, got {n}");
    let c = (2.0 * PI / n as f64).cos();
    (c + 5.0 + ((c + 1.0) * (c + 9.0)).sqrt()) / 16.0
}

/// Value of `χ^k` at a corner of the unit square, plus its first
/// and second partial derivatives. All quantities are in the
/// 2-D domain Ω (R²).
#[derive(Debug, Clone, Copy)]
pub struct ChiJet {
    /// `χ^k(u, v)` evaluated at the corner.
    pub pos: (f64, f64),
    /// `∂χ^k / ∂u` at the corner.
    pub d_u: (f64, f64),
    /// `∂χ^k / ∂v` at the corner.
    pub d_v: (f64, f64),
    /// `∂²χ^k / ∂u²` at the corner.
    pub d_uu: (f64, f64),
    /// `∂²χ^k / ∂v²` at the corner.
    pub d_vv: (f64, f64),
    /// `∂²χ^k / ∂u∂v` at the corner.
    pub d_uv: (f64, f64),
}

/// Bilinear approximation `χ^k(u, v) := u · m_k + v · m_{k+1}` to
/// the Catmull-Clark characteristic map.
///
/// Per `m_k = (cos(k · 2π/n), sin(k · 2π/n))`.
#[derive(Debug, Clone, Copy)]
pub struct BilinearChi {
    pub n: usize,
}

/// Which corner of the unit square to query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChiCorner {
    /// `(0, 0)` — star point.
    Star,
    /// `(1, 0)` — `m_k`.
    Mk,
    /// `(0, 1)` — `m_{k+1}`.
    Mk1,
    /// `(1, 1)` — `v_k`.
    Vk,
}

impl ChiCorner {
    #[inline]
    pub fn uv(self) -> (f64, f64) {
        match self {
            ChiCorner::Star => (0.0, 0.0),
            ChiCorner::Mk => (1.0, 0.0),
            ChiCorner::Mk1 => (0.0, 1.0),
            ChiCorner::Vk => (1.0, 1.0),
        }
    }
}

impl BilinearChi {
    /// New bilinear χ for valence `n ≥ 3`.
    pub fn new(n: usize) -> Self {
        assert!(n >= 3, "bilinear χ needs n ≥ 3, got {n}");
        Self { n }
    }

    /// Evaluate `χ^k(u, v) = u · m_k + v · m_{k+1}`.
    #[inline]
    pub fn eval(&self, k: usize, u: f64, v: f64) -> (f64, f64) {
        let (mkx, mky) = sector_corner(self.n, k);
        let (mk1x, mk1y) = sector_corner(self.n, k + 1);
        (u * mkx + v * mk1x, u * mky + v * mk1y)
    }

    /// Full `ChiJet` at the given corner of sector `k`. For the
    /// bilinear map, second partials are zero and first partials
    /// are constant.
    pub fn jet(&self, k: usize, corner: ChiCorner) -> ChiJet {
        let (mkx, mky) = sector_corner(self.n, k);
        let (mk1x, mk1y) = sector_corner(self.n, k + 1);
        let (u, v) = corner.uv();
        ChiJet {
            pos: (u * mkx + v * mk1x, u * mky + v * mk1y),
            d_u: (mkx, mky),
            d_v: (mk1x, mk1y),
            d_uu: (0.0, 0.0),
            d_vv: (0.0, 0.0),
            d_uv: (0.0, 0.0),
        }
    }
}

/// CC-consistent χ with subdominant-eigenvalue-scaled partials at
/// the star corner. See the module docs for the construction
/// details and fidelity caveats.
#[derive(Debug, Clone, Copy)]
pub struct CatmullClarkChi {
    pub n: usize,
    pub lambda: f64,
}

impl CatmullClarkChi {
    /// New CC-consistent χ for valence `n ≥ 3`. Computes `λ_n`.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            lambda: cc_subdominant_eigenvalue(n),
        }
    }

    /// Sector-`k` corner positions (same as bilinear).
    fn corner_pos(&self, k: usize, corner: ChiCorner) -> (f64, f64) {
        let (mkx, mky) = sector_corner(self.n, k);
        let (mk1x, mk1y) = sector_corner(self.n, k + 1);
        match corner {
            ChiCorner::Star => (0.0, 0.0),
            ChiCorner::Mk => (mkx, mky),
            ChiCorner::Mk1 => (mk1x, mk1y),
            ChiCorner::Vk => (mkx + mk1x, mky + mk1y),
        }
    }

    /// Evaluate `χ^k(u, v)` via the bi-cubic Bézier construction
    /// (subdominant-scaled at the star, bilinear elsewhere).
    pub fn eval(&self, k: usize, u: f64, v: f64) -> (f64, f64) {
        let (mkx, mky) = sector_corner(self.n, k);
        let (mk1x, mk1y) = sector_corner(self.n, k + 1);
        let l = self.lambda;
        // CP grid P_{i,j} for i, j in 0..=3. Bilinear baseline
        //   P_{i,j} = (i/3)·m_k + (j/3)·m_{k+1}
        // except (1, 0) and (0, 1) which use λ/3 instead of 1/3.
        let cp = |i: usize, j: usize| -> (f64, f64) {
            let ai = if i == 1 && j == 0 {
                l / 3.0
            } else {
                i as f64 / 3.0
            };
            let aj = if i == 0 && j == 1 {
                l / 3.0
            } else {
                j as f64 / 3.0
            };
            (ai * mkx + aj * mk1x, ai * mky + aj * mk1y)
        };
        let bu = cubic_bernstein(u);
        let bv = cubic_bernstein(v);
        let mut x = 0.0;
        let mut y = 0.0;
        for (i, bu_i) in bu.iter().enumerate() {
            for (j, bv_j) in bv.iter().enumerate() {
                let w = bu_i * bv_j;
                let (cx, cy) = cp(i, j);
                x += w * cx;
                y += w * cy;
            }
        }
        (x, y)
    }

    /// Full `ChiJet` at the given corner of sector `k`.
    pub fn jet(&self, k: usize, corner: ChiCorner) -> ChiJet {
        let (mkx, mky) = sector_corner(self.n, k);
        let (mk1x, mk1y) = sector_corner(self.n, k + 1);
        let l = self.lambda;
        let pos = self.corner_pos(k, corner);
        let (d_u, d_v, d_uu, d_vv, d_uv) = match corner {
            ChiCorner::Star => {
                // ∂u χ = λ·m_k, ∂v χ = λ·m_{k+1}
                // ∂²u χ = 4(1-λ)·m_k, ∂²v χ = 4(1-λ)·m_{k+1}
                // ∂²uv χ = 3(1-λ)·(m_k + m_{k+1})
                let s1 = 1.0 - l;
                (
                    (l * mkx, l * mky),
                    (l * mk1x, l * mk1y),
                    (4.0 * s1 * mkx, 4.0 * s1 * mky),
                    (4.0 * s1 * mk1x, 4.0 * s1 * mk1y),
                    (3.0 * s1 * (mkx + mk1x), 3.0 * s1 * (mky + mk1y)),
                )
            }
            // Non-star corners: match bilinear (partials = constant
            // m_k / m_{k+1}; second partials zero).
            _ => ((mkx, mky), (mk1x, mk1y), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)),
        };
        ChiJet {
            pos,
            d_u,
            d_v,
            d_uu,
            d_vv,
            d_uv,
        }
    }
}

/// Cubic Bernstein basis `[B_0^3(t), B_1^3(t), B_2^3(t), B_3^3(t)]`.
#[inline]
fn cubic_bernstein(t: f64) -> [f64; 4] {
    let s = 1.0 - t;
    [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t]
}

#[cfg(test)]
#[path = "kp13_bi5_chi_tests.rs"]
mod tests;
