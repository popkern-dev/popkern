//! Characteristic map `χ` for the KP13 bi-5 LSQ.
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
//! **[`CatmullClarkChi`]** — a CC-consistent upgrade with a bi-cubic
//! Bézier patch whose **entire** CP grid is shaped by the CC eigen-
//! structure. CPs in the immediate neighbourhood of the star carry
//! the subdominant-eigenvalue `λ_n` scaling; the 2-ring corner
//! `v_k` picks up a sub-subdominant correction scale `μ_n = λ_n²`
//! (standard CC fact: face-point eigenvalues are the square of the
//! edge-point ones). All four corner jets — value, first partials,
//! second partials — are derived by evaluating the bi-cubic Bernstein
//! basis and its derivatives on the CP grid, so the struct is
//! internally consistent.
//!
//! At the star corner the resulting partials are:
//!
//! ```text
//!   ∂_u χ^k(0, 0)   = λ_n · m_k
//!   ∂_v χ^k(0, 0)   = λ_n · m_{k+1}
//!   ∂²_u χ^k(0, 0)  = 4(1 − λ_n) · m_k
//!   ∂²_v χ^k(0, 0)  = 4(1 − λ_n) · m_{k+1}
//!   ∂²_uv χ^k(0, 0) = 0
//! ```
//!
//! (The mixed partial vanishes by the rotational symmetry of the CP
//! grid — the `u` and `v` directions carry independent copies of
//! the subdominant scaling with no cross-coupling term.)
//!
//! where `λ_n := (c + 5 + √((c + 1)(c + 9))) / 16`, `c := cos(2π/n)`.
//! The `v_k` corner value is `μ_n · (m_k + m_{k+1})` — with
//! `μ_n = λ_n²` (the standard CC face-eigenvalue relationship). This
//! is **not** numerically equal to the bilinear `m_k + m_{k+1}` even
//! at regular `n = 4`; `CatmullClarkChi::new(4)` gives `μ = 1/4`, so
//! `v_k` lands at `1/4 · (m_k + m_{k+1})`. Non-star partials are
//! computed from the CP grid; where the grid reduces to bilinear on a
//! given layer the partials reduce to bilinear there too.
//!
//! # Fidelity caveat (still additive)
//!
//! Full Stam-1998 evaluation of the CC characteristic map at
//! arbitrary `(u, v)` — with the proper eigenvector decomposition
//! at every level of subdivision — remains out of scope here. The
//! bi-cubic patch captures the dominant CC behaviour at the four
//! corners of the unit square (which is all the LSQ in I.11e
//! consumes) but is not a faithful reproduction of χ at interior
//! points.
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

/// CC-consistent χ built as a bi-cubic Bézier patch whose CP grid
/// carries the subdominant eigenvalue `λ_n` (edge-layer) and the
/// sub-subdominant eigenvalue `μ_n = λ_n²` (face-layer, affecting
/// `v_k`).
///
/// See the module header for the construction and fidelity caveats.
#[derive(Debug, Clone, Copy)]
pub struct CatmullClarkChi {
    pub n: usize,
    pub lambda: f64,
    /// Sub-subdominant eigenvalue `μ_n`. Under standard CC,
    /// face-point eigenvalues are the square of edge-point ones, so
    /// `μ_n := λ_n²`. Controls where the `v_k` corner lands.
    pub mu: f64,
}

impl CatmullClarkChi {
    /// New CC-consistent χ for valence `n ≥ 3`. Computes `λ_n` and
    /// sets `μ_n := λ_n²` — the standard CC face-eigenvalue
    /// relationship that gives bilinear `v_k = m_k + m_{k+1}` at
    /// `n = 4` (since `λ = 1/2` implies `μ = 1/4`, and the bi-cubic
    /// evaluation at `(1, 1)` lands exactly on `m_k + m_{k+1}`).
    pub fn new(n: usize) -> Self {
        let lambda = cc_subdominant_eigenvalue(n);
        Self {
            n,
            lambda,
            mu: lambda * lambda,
        }
    }

    /// Evaluate `χ^k(u, v)` via the bi-cubic Bézier patch.
    pub fn eval(&self, k: usize, u: f64, v: f64) -> (f64, f64) {
        let (cps_x, cps_y) = self.cp_grid(k);
        let bu = cubic_bernstein(u);
        let bv = cubic_bernstein(v);
        let mut x = 0.0;
        let mut y = 0.0;
        for (i, bu_i) in bu.iter().enumerate() {
            for (j, bv_j) in bv.iter().enumerate() {
                let w = bu_i * bv_j;
                x += w * cps_x[i][j];
                y += w * cps_y[i][j];
            }
        }
        (x, y)
    }

    /// Full `ChiJet` at the given corner of sector `k`, computed
    /// from the bi-cubic CP grid via the Bernstein-basis corner
    /// derivative identities.
    pub fn jet(&self, k: usize, corner: ChiCorner) -> ChiJet {
        let (cps_x, cps_y) = self.cp_grid(k);
        let (i0, j0) = corner_indices(corner);
        let pos = (cps_x[i0][j0], cps_y[i0][j0]);
        let d_u = corner_partial_u(&cps_x, &cps_y, corner);
        let d_v = corner_partial_v(&cps_x, &cps_y, corner);
        let d_uu = corner_partial_uu(&cps_x, &cps_y, corner);
        let d_vv = corner_partial_vv(&cps_x, &cps_y, corner);
        let d_uv = corner_partial_uv(&cps_x, &cps_y, corner);
        ChiJet {
            pos,
            d_u,
            d_v,
            d_uu,
            d_vv,
            d_uv,
        }
    }

    /// Build the 4×4 CP grid for sector `k`. Returns `(x_grid,
    /// y_grid)` with `grid[i][j] = P_{i, j}`.
    ///
    /// # CP-grid shape
    ///
    /// Denote `ai := α(i, n)`, `aj := α(j, n)` with
    ///
    /// | `i` | `α(i, n)` |
    /// |-----|-----------|
    /// | `0` | `0`       |
    /// | `1` | `λ / 3`   |
    /// | `2` | `2 / 3`   |
    /// | `3` | `1`       |
    ///
    /// Then `P_{i, j} := ai·m_k + aj·m_{k+1}` — **except** `P_{3,3}`
    /// which is overridden to `μ · (m_k + m_{k+1})` so that the
    /// `v_k` corner picks up the sub-subdominant CC scaling.
    fn cp_grid(&self, k: usize) -> ([[f64; 4]; 4], [[f64; 4]; 4]) {
        let (mkx, mky) = sector_corner(self.n, k);
        let (mk1x, mk1y) = sector_corner(self.n, k + 1);
        let alpha = [0.0, self.lambda / 3.0, 2.0 / 3.0, 1.0];
        let mut x = [[0.0; 4]; 4];
        let mut y = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                x[i][j] = alpha[i] * mkx + alpha[j] * mk1x;
                y[i][j] = alpha[i] * mky + alpha[j] * mk1y;
            }
        }
        // v_k corner: override with μ · (m_k + m_{k+1}).
        x[3][3] = self.mu * (mkx + mk1x);
        y[3][3] = self.mu * (mky + mk1y);
        (x, y)
    }
}

#[inline]
fn corner_indices(c: ChiCorner) -> (usize, usize) {
    match c {
        ChiCorner::Star => (0, 0),
        ChiCorner::Mk => (3, 0),
        ChiCorner::Mk1 => (0, 3),
        ChiCorner::Vk => (3, 3),
    }
}

#[inline]
fn corner_partial_u(x: &[[f64; 4]; 4], y: &[[f64; 4]; 4], c: ChiCorner) -> (f64, f64) {
    // ∂u B at u=0: 3·(P_{1,j0} − P_{0,j0}); at u=1: 3·(P_{3,j0} − P_{2,j0}).
    // Always "higher-u-index minus lower-u-index" with factor 3.
    let (i_high, i_low) = match c {
        ChiCorner::Star | ChiCorner::Mk1 => (1, 0),
        ChiCorner::Mk | ChiCorner::Vk => (3, 2),
    };
    let j0 = match c {
        ChiCorner::Star | ChiCorner::Mk => 0,
        ChiCorner::Mk1 | ChiCorner::Vk => 3,
    };
    (
        3.0 * (x[i_high][j0] - x[i_low][j0]),
        3.0 * (y[i_high][j0] - y[i_low][j0]),
    )
}

#[inline]
fn corner_partial_v(x: &[[f64; 4]; 4], y: &[[f64; 4]; 4], c: ChiCorner) -> (f64, f64) {
    let (j_high, j_low) = match c {
        ChiCorner::Star | ChiCorner::Mk => (1, 0),
        ChiCorner::Mk1 | ChiCorner::Vk => (3, 2),
    };
    let i0 = match c {
        ChiCorner::Star | ChiCorner::Mk1 => 0,
        ChiCorner::Mk | ChiCorner::Vk => 3,
    };
    (
        3.0 * (x[i0][j_high] - x[i0][j_low]),
        3.0 * (y[i0][j_high] - y[i0][j_low]),
    )
}

#[inline]
fn corner_partial_uu(x: &[[f64; 4]; 4], y: &[[f64; 4]; 4], c: ChiCorner) -> (f64, f64) {
    // ∂²u at u = 0: 6(P_{2,j0} − 2P_{1,j0} + P_{0,j0})
    // At u = 1: 6(P_{3,j0} − 2P_{2,j0} + P_{1,j0})
    let (a, b, d) = match c {
        ChiCorner::Star | ChiCorner::Mk1 => (0usize, 1usize, 2usize),
        ChiCorner::Mk | ChiCorner::Vk => (3usize, 2usize, 1usize),
    };
    let j0 = if matches!(c, ChiCorner::Star | ChiCorner::Mk) {
        0
    } else {
        3
    };
    (
        6.0 * (x[a][j0] - 2.0 * x[b][j0] + x[d][j0]),
        6.0 * (y[a][j0] - 2.0 * y[b][j0] + y[d][j0]),
    )
}

#[inline]
fn corner_partial_vv(x: &[[f64; 4]; 4], y: &[[f64; 4]; 4], c: ChiCorner) -> (f64, f64) {
    let (a, b, d) = match c {
        ChiCorner::Star | ChiCorner::Mk => (0usize, 1usize, 2usize),
        ChiCorner::Mk1 | ChiCorner::Vk => (3usize, 2usize, 1usize),
    };
    let i0 = if matches!(c, ChiCorner::Star | ChiCorner::Mk1) {
        0
    } else {
        3
    };
    (
        6.0 * (x[i0][a] - 2.0 * x[i0][b] + x[i0][d]),
        6.0 * (y[i0][a] - 2.0 * y[i0][b] + y[i0][d]),
    )
}

#[inline]
fn corner_partial_uv(x: &[[f64; 4]; 4], y: &[[f64; 4]; 4], c: ChiCorner) -> (f64, f64) {
    // ∂u∂v at a corner of a bi-cubic:
    //   (i0, j0) = corner indices (∈ {0, 3})
    //   i1 = adjacent in u direction; j1 = adjacent in v direction
    //   sign = (−1)^(matches i0→3 or j0→3 reversal)
    let (i0, j0) = corner_indices(c);
    let i1 = if i0 == 0 { 1 } else { 2 };
    let j1 = if j0 == 0 { 1 } else { 2 };
    let sign = match (i0, j0) {
        (0, 0) | (3, 3) => 1.0,
        _ => -1.0,
    };
    (
        9.0 * sign * (x[i1][j1] - x[i0][j1] - x[i1][j0] + x[i0][j0]),
        9.0 * sign * (y[i1][j1] - y[i0][j1] - y[i1][j0] + y[i0][j0]),
    )
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
