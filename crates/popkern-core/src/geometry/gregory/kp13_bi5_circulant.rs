//! KP09a §2.3(b) circulant system for `g^i_{21}`.
//!
//! # Summary
//!
//! After the star-point quadratic identity (§2.3(a), implemented in
//! `kp13_bi5_star_quadratic`) pins the 6 near-star coefs, one more
//! layer of BB-coefs needs to be set to enforce C² continuity across
//! sector boundaries. Paper [KP09a §2.3(b)]:
//!
//! > "The coefficients `g^i_{30}` can be chosen freely. If the
//! > coefficients `g^i_{21}` satisfy the circulant system (8), ...
//! > then setting the coefficients `g^i_{12}` according to (3)
//! > enforces the C² constraints."
//!
//! # Circulant system (KP09a eq 8)
//!
//! For each sector `i = 0..n−1`:
//!
//! ```text
//!   g^{i-1}_{21}  +  4c · g^i_{21}  +  g^{i+1}_{21}  =  2 · R^i
//! ```
//!
//! where
//!
//! ```text
//!   R^i := 2c² · g^i_{30}
//!        + c  · g^i_{03}
//!        + 4cc̄ · g^i_{20}
//!        − 2c̄ · g^i_{11}
//!        +  c̄ · g^i_{02}
//!        + 2c̄² · g^i_{10}
//!   c   := cos(2π / n)
//!   c̄   := 1 − c
//! ```
//!
//! All RHS CPs are in sector `i`'s paper 2-index convention:
//! `g^i_{a,b}` = our flat index `flat_index(5, a, b)`.
//!
//! | paper `(a, b)` | flat | where it sits               |
//! |----------------|-----:|-----------------------------|
//! | `(3, 0)`       |   15 | free param (interior)       |
//! | `(0, 3)`       |    3 | shared edge with sector i+1 |
//! | `(2, 0)`       |   11 | near-star (from quadratic)  |
//! | `(1, 1)`       |    7 | near-star                   |
//! | `(0, 2)`       |    2 | near-star                   |
//! | `(1, 0)`       |    6 | near-star                   |
//! | `(2, 1)` (LHS) |   12 | this module solves for      |
//!
//! # DFT closed form (Lemma 1)
//!
//! For `n ≥ 5`, `n ≠ 6`:
//!
//! ```text
//!   g^i_{21} = (1/n) · Σ_{l=0..n−1} Σ_{j=0..n−1}
//!                R^l · cos((i−l)·j·α) / (2c + cos(j·α))
//! ```
//!
//! The denominator `2c + cos(j·α)` vanishes at `j = n/2` iff
//! `n = 6` (so `c = 1/2`, `cos(jα) = cos(π) = −1`). Paper excludes
//! `n ∈ {3, 4, 6}` from this construction.
//!
//! # Back-fill `g^i_{12}` (KP09a eq 3)
//!
//! Once `g^{i+1}_{21}` is solved, KP09a eq (3) at `(j=2, k=1)` gives
//!
//! ```text
//!   g^i_{12} = ν₀ · g^i_{02}  +  ν₁ · g^i_{03}  −  g^{i+1}_{21}
//! ```
//!
//! (rewriting `p^{i+1}_{2, 1} = ν₀·p^i_{0, 2} + ν₁·p^i_{0, 3} + ν₂·p^i_{1, 2}`
//! and solving for `p^i_{1, 2}` with `ν₂ = −1`.)
//!
//! Paper labels `g^i_{12}` = flat 8.

#![allow(dead_code)]

use super::kp13_bi5_gtri5::Gtri5;
use crate::geometry::Point3;
use std::f64::consts::PI;

use super::kp13_bi5_connect::paper_flat;

/// Error from attempting to solve the circulant system.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CirculantError {
    /// Valence `n` is invalid (below 3) — the KP13 pipeline requires
    /// at least a valence-3 extraordinary vertex. Paper § 6 covers
    /// `n ∈ {3, 4}` via min-norm pseudo-inverse; this variant fires
    /// only for truly degenerate `n < 3`.
    UnsupportedValence(usize),
}

impl std::fmt::Display for CirculantError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CirculantError::UnsupportedValence(n) => {
                write!(f, "circulant: valence n = {n} < 3 is not a valid EV")
            }
        }
    }
}

impl std::error::Error for CirculantError {}

/// Compute `R^i` for a single sector per KP09a eq (8).
///
/// `sector` supplies the six input CPs at paper indices `(1, 0)`,
/// `(0, 2)`, `(1, 1)`, `(2, 0)`, `(0, 3)`, `(3, 0)`. The caller is
/// responsible for having populated these consistently (via the
/// star-quadratic for the near-star four, and via `g^i_{30}` /
/// `g^{i+1}_{30}` for the edge-layer two).
pub fn compute_r_i(sector: &super::kp13_bi5_gtri5::Gtri5Sector, n: usize) -> Point3 {
    let c = (2.0 * PI / n as f64).cos();
    let c_bar = 1.0 - c;
    let g_30 = sector.cps[paper_flat(3, 0)];
    let g_03 = sector.cps[paper_flat(0, 3)];
    let g_20 = sector.cps[paper_flat(2, 0)];
    let g_11 = sector.cps[paper_flat(1, 1)];
    let g_02 = sector.cps[paper_flat(0, 2)];
    let g_10 = sector.cps[paper_flat(1, 0)];
    let coefs = [
        (2.0 * c * c, g_30),
        (c, g_03),
        (4.0 * c * c_bar, g_20),
        (-2.0 * c_bar, g_11),
        (c_bar, g_02),
        (2.0 * c_bar * c_bar, g_10),
    ];
    weighted_sum(&coefs)
}

/// Solve the circulant system for all `g^i_{21}` and write results
/// into the supplied `Gtri5`.
///
/// Preconditions:
/// * Near-star coefs populated (e.g., via the star-quadratic).
/// * `g^i_{30}` set for each sector.
/// * `g^i_{03}` set to equal sector `i+1`'s `g^{i+1}_{30}` along the
///   shared edge (use [`sync_edge_layer_03_and_30`] if needed).
///
/// Writes `g^i_{21}` (paper `(2, 1)` = flat 12) for each sector.
///
/// # Valence support
///
/// * `n ∈ {5, 7, 8, 9, 10, ...}` — regular Lemma 1 closed form.
///   All DFT denominators `2c + cos(jα)` are non-zero so the
///   recurrence is uniquely solvable.
/// * `n ∈ {3, 4, 6}` — rank-deficient cases. Each has one or more
///   DFT modes where `2c + cos(jα) = 0`, making the circulant
///   matrix singular. We use the **minimum-norm pseudo-inverse**
///   solution, obtained by skipping the singular modes in
///   Lemma 1's sum. Specifically:
///   - `n = 3`: `c = −1/2`, singular mode `j = 0` (kernel =
///     constant sector).
///   - `n = 4`: `c = 0`, singular modes `j ∈ {1, 3}` (kernel =
///     the two alternating `(1, 0, −1, 0)` / `(0, 1, 0, −1)`
///     patterns).
///   - `n = 6`: `c = 1/2`, singular mode `j = 3` (kernel =
///     `(1, −1, 1, −1, 1, −1)`).
///
///   The pseudo-inverse gives the LSQ-equivalent minimum-norm
///   solution over the free DoFs. Paper KP09a § 6 prescribes
///   guided-surfacing / `F_5` minimisation for these free DoFs
///   (see `kp15_bi6_cap_f5_solve` for the bi-6 polar variant of
///   the same machinery); we default to min-norm here because it
///   is invariant, deterministic, and consistent with the LSQ
///   residual treatment at regular valences. If the RHS is
///   orthogonal to all null modes, the pseudo-inverse result is
///   the exact solution.
///
/// # Errors
///
/// Returns [`CirculantError::UnsupportedValence`] only for `n < 3`
/// — those are not valid extraordinary vertices.
pub fn solve_g21_circulant(gtri5: &mut Gtri5) -> Result<(), CirculantError> {
    let n = gtri5.n;
    if n < 3 {
        return Err(CirculantError::UnsupportedValence(n));
    }

    let c = (2.0 * PI / n as f64).cos();
    let alpha = 2.0 * PI / n as f64;
    // Threshold below which we treat `2c + cos(jα)` as a zero
    // (null-mode) eigenvalue. At n = 6 the zero is exact; other
    // valences give strictly nonzero denominators, so this only
    // kicks in for the true pseudo-inverse case.
    let zero_tol = 1e-12;

    // Compute R^i per sector up front.
    let r: Vec<Point3> = (0..n).map(|i| compute_r_i(gtri5.sector(i), n)).collect();

    // Apply Lemma 1 per-output to compute g^i_{21}. Because the
    // formula is linear and coordinate-wise, we compute each axis
    // independently. For singular modes (n = 6 at j = n/2), we
    // skip the term — the remaining partial sum is the min-norm
    // pseudo-inverse solution.
    let g21: Vec<Point3> = (0..n)
        .map(|i| {
            let mut x = 0.0;
            let mut y = 0.0;
            let mut z = 0.0;
            for (l, r_l) in r.iter().enumerate() {
                for j in 0..n {
                    let denom = 2.0 * c + (j as f64 * alpha).cos();
                    if denom.abs() < zero_tol {
                        continue;
                    }
                    let coef = ((i as f64 - l as f64) * j as f64 * alpha).cos() / denom;
                    x += r_l.x * coef;
                    y += r_l.y * coef;
                    z += r_l.z * coef;
                }
            }
            Point3::new(x / n as f64, y / n as f64, z / n as f64)
        })
        .collect();

    let idx_21 = paper_flat(2, 1);
    for (i, g) in g21.into_iter().enumerate() {
        gtri5.sector_mut(i).cps[idx_21] = g;
    }

    Ok(())
}

/// Residual of the circulant three-term recurrence at sector `i`:
///
/// ```text
///   residual^i := (g^{i-1}_{21} + 4c · g^i_{21} + g^{i+1}_{21}) − 2R^i
/// ```
///
/// For `n ∈ {5, 7, 8, 9, ...}` this is zero after [`solve_g21_circulant`].
///
/// For `n = 6` the min-norm pseudo-inverse leaves a residual of
/// `(−1)^i · λ` where `λ = (1/6) Σ (−1)^l R^l` is the null-space
/// component of R — that is, the residual is a pure j = 3 DFT mode.
/// Interpret a non-zero residual as "input data isn't compatible
/// with the n = 6 circulant" rather than a solver bug.
pub fn circulant_residual(gtri5: &Gtri5, i: usize) -> Point3 {
    let n = gtri5.n;
    let c = (2.0 * PI / n as f64).cos();
    let idx_21 = paper_flat(2, 1);
    let prev = (i + n - 1) % n;
    let next = (i + 1) % n;
    let gm = gtri5.sector(prev).cps[idx_21];
    let gi = gtri5.sector(i).cps[idx_21];
    let gp = gtri5.sector(next).cps[idx_21];
    let r = compute_r_i(gtri5.sector(i), n);
    Point3::new(
        gm.x + 4.0 * c * gi.x + gp.x - 2.0 * r.x,
        gm.y + 4.0 * c * gi.y + gp.y - 2.0 * r.y,
        gm.z + 4.0 * c * gi.z + gp.z - 2.0 * r.z,
    )
}

/// After `g^i_{21}` is solved, back-fill `g^i_{12}` (paper `(1, 2)`,
/// flat 8) via KP09a eq (3) rearranged:
///
/// ```text
///   g^i_{12} = ν₀ · g^i_{02}  +  ν₁ · g^i_{03}  −  g^{i+1}_{21}
/// ```
///
/// Assumes `g^i_{02}`, `g^i_{03}`, `g^{i+1}_{21}` already populated.
pub fn fill_g12_from_g21(gtri5: &mut Gtri5) {
    let n = gtri5.n;
    let c = (2.0 * PI / n as f64).cos();
    let c_bar = 1.0 - c;
    let nu0 = 2.0 * c_bar;
    let nu1 = 2.0 * c;

    let idx_02 = paper_flat(0, 2);
    let idx_03 = paper_flat(0, 3);
    let idx_12 = paper_flat(1, 2);
    let idx_21 = paper_flat(2, 1);

    // Read g^{i+1}_{21} before writing g^i_{12} to avoid aliasing.
    let g_next_21: Vec<Point3> = (0..n)
        .map(|i| gtri5.sector((i + 1) % n).cps[idx_21])
        .collect();

    for (i, g_next) in g_next_21.iter().enumerate() {
        let s = gtri5.sector(i);
        let g_02 = s.cps[idx_02];
        let g_03 = s.cps[idx_03];
        let g_12 = Point3::new(
            nu0 * g_02.x + nu1 * g_03.x - g_next.x,
            nu0 * g_02.y + nu1 * g_03.y - g_next.y,
            nu0 * g_02.z + nu1 * g_03.z - g_next.z,
        );
        gtri5.sector_mut(i).cps[idx_12] = g_12;
    }
}

/// Synchronise the shared-edge CPs at paper `(0, 3)` / `(3, 0)`
/// across adjacent sectors.
///
/// For each sector pair `(i, i+1 mod n)`, sector `i`'s `(0, 3)`
/// (flat 3) must equal sector `i+1`'s `(3, 0)` (flat 15) — they're
/// the same geometric CP on the shared edge.
///
/// Uses sector `i+1`'s `(3, 0)` as authoritative (the "free
/// parameter" of the paper lives on each sector as `g^i_{30}`).
pub fn sync_edge_layer_03_and_30(gtri5: &mut Gtri5) {
    sync_shared_edge_pair(gtri5, (0, 3), (3, 0));
}

/// Synchronise **every** depth-0 shared-edge CP pair across
/// adjacent sectors. The boundary between sector `i` and sector
/// `i+1` is the ray from the star to `m_{i+1}`. Along this edge,
/// sector `i`'s paper `(0, j)` coincides with sector `i+1`'s
/// paper `(j, 0)` for `j = 0..=5`.
///
/// Pairs:
///
/// | j | sector i flat | sector i+1 flat | in-module status |
/// |--:|---------------|-----------------|------------------|
/// | 0 | 0 (star)      | 0 (star)        | `set_star_point` |
/// | 1 | 1 (Q)         | 6 (Q)           | star-quadratic   |
/// | 2 | 2 (Q)         | 11 (Q)          | star-quadratic   |
/// | 3 | 3             | 15 (`g^i_{30}`) | `g^{i+1}_{30}` free |
/// | 4 | 4             | 18              | free outer-edge CP |
/// | 5 | 5 (m_{k+1})   | 20 (m_k)        | `set_mid_edge`   |
///
/// For `j ∈ {3, 4}`, this function copies sector `i+1`'s CP to
/// sector `i`'s matching slot (sector `i+1` is authoritative —
/// it treats that CP as a "free parameter on its own edge").
///
/// The star (`j = 0`) and the mid-edge (`j = 5`) are already
/// maintained by their dedicated invariants.
///
/// For `j ∈ {1, 2}` — near-star — both sides are pinned by the
/// globally-consistent star quadratic (assumed already applied).
/// This function does **not** overwrite them.
pub fn sync_shared_edge_all_depths(gtri5: &mut Gtri5) {
    // j = 3: sector i's (0, 3) ← sector i+1's (3, 0)
    sync_shared_edge_pair(gtri5, (0, 3), (3, 0));
    // j = 4: sector i's (0, 4) ← sector i+1's (4, 0)
    sync_shared_edge_pair(gtri5, (0, 4), (4, 0));
}

/// Copy sector `i+1`'s CP at paper position `authoritative` to
/// sector `i`'s slot at `derived`, for all `i`.
fn sync_shared_edge_pair(
    gtri5: &mut Gtri5,
    derived: (usize, usize),
    authoritative: (usize, usize),
) {
    let n = gtri5.n;
    let idx_derived = paper_flat(derived.0, derived.1);
    let idx_auth = paper_flat(authoritative.0, authoritative.1);
    for i in 0..n {
        let next = (i + 1) % n;
        let auth_cp = gtri5.sector(next).cps[idx_auth];
        gtri5.sector_mut(i).cps[idx_derived] = auth_cp;
    }
}

/// Full pipeline: given a `Gtri5` with near-star coefs and
/// `g^i_{30}` already populated, solve the circulant for `g^i_{21}`,
/// then back-fill `g^i_{12}`.
///
/// Calls [`sync_edge_layer_03_and_30`] before solving so the RHS
/// uses consistent shared-edge values.
pub fn complete_second_derivative_layer(gtri5: &mut Gtri5) -> Result<(), CirculantError> {
    sync_shared_edge_all_depths(gtri5);
    solve_g21_circulant(gtri5)?;
    fill_g12_from_g21(gtri5);
    Ok(())
}

/// Apply KP09a eq (3) at `j = 3` to derive sector `(i+1)`'s
/// `g^{i+1}_{31}` (paper `(3, 1)`, flat 16) from sector `i`'s
/// `g^i_{03}` (flat 3), `g^i_{04}` (flat 4), `g^i_{13}` (flat 9).
///
/// The paper formula `p^{i+1}_{j, 1} = ν₀·p^i_{0, j} +
/// ν₁·p^i_{0, j+1} + ν₂·p^i_{1, j}` at `j = 3` gives the depth-1
/// coefficient two layers inward from the star on the shared-edge
/// side of sector `(i+1)`. Flat 9 of sector `i` is a free parameter
/// under KP13 Fig 4's labelling (no marker), so this equation
/// **derives** flat 16 rather than constrains flat 9.
pub fn fill_g31_via_c1(gtri5: &mut Gtri5) {
    use super::kp13_bi5_connect::{ConnectWeights, c1_join};
    let n = gtri5.n;
    let w = ConnectWeights::for_valence(n);
    let idx_31 = super::kp13_bi5_connect::paper_flat(3, 1);
    let derived: Vec<Point3> = (0..n).map(|i| c1_join(&w, gtri5.sector(i), 3)).collect();
    for (i, p) in derived.into_iter().enumerate() {
        let next = (i + 1) % n;
        gtri5.sector_mut(next).cps[idx_31] = p;
    }
}

/// Apply KP09a eq (3) at `j = 4` to derive sector `(i+1)`'s
/// `g^{i+1}_{41}` (paper `(4, 1)`, flat 19) from sector `i`'s
/// `g^i_{04}` (flat 4), `g^i_{05}` (flat 5 = `m_{i+1}`),
/// `g^i_{14}` (flat 10).
///
/// Analogous to [`fill_g31_via_c1`], but at the m_k end of the
/// shared edge.
pub fn fill_g41_via_c1(gtri5: &mut Gtri5) {
    use super::kp13_bi5_connect::{ConnectWeights, c1_join};
    let n = gtri5.n;
    let w = ConnectWeights::for_valence(n);
    let idx_41 = super::kp13_bi5_connect::paper_flat(4, 1);
    let derived: Vec<Point3> = (0..n).map(|i| c1_join(&w, gtri5.sector(i), 4)).collect();
    for (i, p) in derived.into_iter().enumerate() {
        let next = (i + 1) % n;
        gtri5.sector_mut(next).cps[idx_41] = p;
    }
}

/// Fuller pipeline that also applies the higher-j C¹ derivations
/// (`g^i_{31}` and `g^i_{41}`) after [`complete_second_derivative_layer`].
///
/// This brings g△ to full C¹ along every shared edge (not just near
/// the star). C² along the whole shared edge would additionally
/// require eq (4) at `j ∈ {2, 3}`, which constrains flat 13 / flat
/// 14 coefs — these are marked as **free** in KP13 Fig 4, so we do
/// not enforce them here. Callers seeking the stricter "full-C²"
/// reading can apply eq (4) derivations explicitly via [`super::kp13_bi5_connect::c2_join`].
pub fn complete_derivative_layers_with_higher_j(gtri5: &mut Gtri5) -> Result<(), CirculantError> {
    complete_second_derivative_layer(gtri5)?;
    fill_g31_via_c1(gtri5);
    fill_g41_via_c1(gtri5);
    Ok(())
}

#[inline]
fn weighted_sum(pairs: &[(f64, Point3)]) -> Point3 {
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for &(w, p) in pairs {
        x += w * p.x;
        y += w * p.y;
        z += w * p.z;
    }
    Point3::new(x, y, z)
}

#[cfg(test)]
#[path = "kp13_bi5_circulant_tests.rs"]
mod tests;
