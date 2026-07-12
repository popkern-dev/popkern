//! KP15 bi-6 deep-interior population — the complete deterministic
//! resolution of every dependent BB-coefficient of the polar cap.
//!
//! This is the missing construction that the KP15 resume plan
//! identified as the fold root cause: "deep-interior CPs
//! unpopulated by any code path". It is a faithful port of the
//! KP15 companion Maple script (`docs/papers/n6.txt`), reproduced
//! and independently verified by the sympy oracle
//! `tools/reference/kp15-bi6-cap/n6_oracle.py`. That oracle proves
//! `G1 ≡ G2 ≡ 0` (symbolic) and closes to ~1e-12 numerically at
//! every `n ∈ {5..10}`; this port is validated bit-close against
//! its committed fixture `n6_cap_fixture.json`.
//!
//! # What it resolves
//!
//! Given the reparameterised tensor border (rows/cols of depth
//! `≤ 2`), the central quadratic data, `(c, τ)`, and the paper free
//! set, it assigns every one of the 16 deep-block CPs per sector
//! (script indices `i, j ∈ [3, 6]`) via, in order:
//! - (S1) the central quadratic plus the eq (6) C² sector glue;
//! - (S2) the eq (7) `2n×2n` ring (`p[r][5,4]`, `p[r][4,5]`) and
//!   the eq (8) `n×n` ring (`p[r][4,4]`; `lyc[0] = sal` deferred);
//! - eq (5) the reduced-degree sector curve (`p[r][6,3]`);
//! - the row-6 glue (`p[t][i,6] = p[r][6,i]`);
//! - (S3) the `µ̃`/`µ̄` rules (`p[r][5,3]`, `p[t][3,5]`, `p[t][3,4]`).
//!
//! # Indexing convention (critical)
//!
//! All grids and indices are in the **script** convention: script
//! `p[i, j]` equals paper `p[6 − i, 6 − j]`, and script `[6, 6]`
//! is the extraordinary point. See
//! `tools/reference/kp15-bi6-cap/README.md`.
//!
//! # Valence
//!
//! Valence-general in `(c, τ)` for `n ≥ 5`. The rank-deficient
//! eq (7) case (`n ∈ {3, 6}`, where `p[0][4,5] = P210` is a genuine
//! free parameter) is handled. `n < 5` (the paper §3.5 `n = 3`
//! special case and the rectangular `n = 4` path) declines with a
//! structured diagnostic — those valences are served by other
//! construction paths.

use super::kp_polar_shared::weighted_sum;
use super::kp15_bi6_reparam::Bi6Reparams;
use crate::geometry::Point3;
use border::{bspline_net_to_bi3_bb, reparam_tensor_border_row};
use rings::{Grid, apply_s3_rules, solve_eq7_ring, solve_eq8_ring};

// The tensor-border reparameterisation and the interior-ring
// solvers are implementation details of this construction, kept in
// sibling files only for the line-size budget.
#[path = "kp15_bi6_deep_interior_border.rs"]
mod border;
#[path = "kp15_bi6_deep_interior_rings.rs"]
mod rings;

/// The KP15 free-parameter set (paper §3.6 step 4), script
/// convention. Every other cap CP is determined from these plus
/// the reparameterised input net.
#[derive(Debug, Clone)]
pub struct Bi6DeepInteriorFree {
    /// `P00 = p[s][6,6]` — the extraordinary point (shared by all
    /// sectors).
    pub p00: Point3,
    /// `TAN0 = p[0][5,6]`, `TAN1 = p[0][6,5]` — central first-order
    /// tangent data on sector 0.
    pub tan: [Point3; 2],
    /// `TWO0 = p[0][4,6]`, `TWO1 = p[0][5,5]`, `TWO2 = p[0][6,4]` —
    /// central quadratic data on sector 0.
    pub two: [Point3; 3],
    /// `P220 = p[0][4,4]` — free at every valence.
    pub p220: Point3,
    /// `P210 = p[0][4,5]` — free only at `n ∈ {3, 6}` (eq (7) ring
    /// rank-deficient); ignored otherwise.
    pub p210: Point3,
    /// Per-sector `p[r][4,3]` (`r = 0..n`).
    pub p23: Vec<Point3>,
    /// Per-sector `p[r][3,3]` (`r = 0..n`).
    pub p33: Vec<Point3>,
}

/// Result of [`build_bi6_cap_from_net`]: the per-sector bi-6 CP
/// grids in script convention, plus the deferred eq (8) closure
/// residual `sal` (`lyc[0]`, which the paper proves vanishes — a
/// nonzero value here is a diagnostic that the input net was not
/// C²-consistent).
#[derive(Debug, Clone)]
pub struct Bi6CapScriptGrids {
    /// Per-sector `grid[s][i][j]`, `i, j ∈ 0..=6`, script convention.
    pub sectors: Vec<[[Point3; 7]; 7]>,
    /// Deferred eq (8) residual `sal = lyc[0]`; ~0 for a valid net.
    pub sal: Point3,
}

const ZERO: Point3 = Point3 {
    x: 0.0,
    y: 0.0,
    z: 0.0,
};

/// Assemble the full per-sector bi-6 cap grids (script convention)
/// from the input B-spline net and the paper free set, resolving
/// every determined CP per the KP15 construction.
///
/// `net[s]` is the `4 × 4` bicubic B-spline control net for sector
/// `s` (`net.len()` must equal `params.n`). `free` supplies the
/// paper free-parameter values (`p23`/`p33` length `params.n`).
///
/// Returns an error for `n < 5` (declined; see module docs), for a
/// net/free length mismatch, or if an interior ring system is
/// singular (should not happen for well-formed input).
pub fn build_bi6_cap_from_net(
    params: &Bi6Reparams,
    net: &[[[Point3; 4]; 4]],
    free: &Bi6DeepInteriorFree,
) -> Result<Bi6CapScriptGrids, String> {
    let n = params.n;
    if n < 5 {
        return Err(format!(
            "KP15 bi-6 deep-interior population is valence-general for n ≥ 5; \
             got n = {n} (n = 3 is the paper §3.5 special case, n = 4 is served \
             by the rectangular path — both decline here)"
        ));
    }
    params.validate()?;
    if net.len() != n {
        return Err(format!("net length {} ≠ n = {n}", net.len()));
    }
    if free.p23.len() != n || free.p33.len() != n {
        return Err(format!(
            "free.p23/p33 length ({}, {}) ≠ n = {n}",
            free.p23.len(),
            free.p33.len()
        ));
    }

    let c = params.c;
    let tau = params.tau;
    let d2 = params.d2;
    let sf = |k: usize| k % n;

    let mut g: Vec<Grid> = vec![[[ZERO; 7]; 7]; n];

    // ---- Border: reparameterised tensor border (depth ≤ 2) ----
    for (s, sector_net) in net.iter().enumerate() {
        let b = bspline_net_to_bi3_bb(sector_net);
        let tad = reparam_tensor_border_row(|i, j| b[i][j], c, d2, tau);
        let tak = reparam_tensor_border_row(|i, j| b[j][i], c, d2, tau);
        for i in 0..7 {
            for j in 0..3 {
                g[s][i][j] = tad[j][i];
                g[s][j][i] = tak[j][i];
            }
        }
    }

    // ---- (S1) central quadratic + eq (6) C² sector glue ----
    for grid in g.iter_mut() {
        grid[6][6] = free.p00;
    }
    g[0][5][6] = free.tan[0];
    g[0][6][5] = free.tan[1];
    g[0][4][6] = free.two[0];
    g[0][5][5] = free.two[1];
    g[0][6][4] = free.two[2];
    g[0][4][4] = free.p220;
    for (r, grid) in g.iter_mut().enumerate() {
        grid[4][3] = free.p23[r];
        grid[3][3] = free.p33[r];
    }
    for s in 0..n - 1 {
        let t = sf(s + 1);
        let (p56, p65, p66, p64, p55, p46) = (
            g[s][5][6], g[s][6][5], g[s][6][6], g[s][6][4], g[s][5][5], g[s][4][6],
        );
        g[t][5][6] = p65;
        g[t][4][6] = p64;
        g[t][6][5] = weighted_sum(&[(-1.0, p56), (2.0 * c, p65), (2.0 * (1.0 - c), p66)]);
        g[t][5][5] = weighted_sum(&[
            (5.0 * c / 3.0, p64),
            (c / 3.0, p66),
            (-2.0 * c, p65),
            (-1.0, p55),
            (2.0, p65),
        ]);
        g[t][6][4] = weighted_sum(&[
            (24.0 / 5.0 * tau * c - 4.0, p56),
            (
                -24.0 / 5.0 * tau * c * c + 44.0 / 5.0 * c - 24.0 / 5.0 * c * c,
                p65,
            ),
            (
                -24.0 / 5.0 * tau * c + 24.0 / 5.0 * tau * c * c + 4.0 - 4.0 * c
                    + 4.0 / 5.0 * c * c,
                p66,
            ),
            (1.0, p46),
            (4.0 * c * c, p64),
            (-24.0 / 5.0 * c, p55),
        ]);
    }

    // ---- (S2) eq (7): p[r][5,4], p[r][4,5] ring ----
    solve_eq7_ring(&mut g, params, free, n)?;

    // ---- (S2) eq (8): p[r][4,4] ring; lyc[0] = sal deferred ----
    let sal = solve_eq8_ring(&mut g, params, n)?;

    // ---- eq (5) reduced-degree sector curve p[r][6,3] ----
    for grid in g.iter_mut() {
        grid[6][3] = weighted_sum(&[
            (0.75, grid[6][4]),
            (0.75, grid[6][2]),
            (0.05, grid[6][6]),
            (-0.3, grid[6][5]),
            (-0.3, grid[6][1]),
            (0.05, grid[6][0]),
        ]);
    }

    // ---- row-6 glue: p[t][i,6] = p[r][6,i] ----
    for r in 0..n {
        let t = sf(r + 1);
        let row6: [Point3; 7] = g[r][6];
        for i in 0..7 {
            g[t][i][6] = row6[i];
        }
    }

    // ---- (S3) µ̃ / µ̄ rules ----
    apply_s3_rules(&mut g, params, n);

    Ok(Bi6CapScriptGrids { sectors: g, sal })
}

#[cfg(test)]
#[path = "kp15_bi6_deep_interior_tests.rs"]
mod tests;
