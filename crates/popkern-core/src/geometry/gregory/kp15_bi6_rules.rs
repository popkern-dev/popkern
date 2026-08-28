//! K-P 2015 bi-6 BB-coefficient assignment rules.
//!
//! Source: `docs/journal/polar-kp-math.md` §1.5 (Catmull-Clark
//! star point in `kp_polar_shared`), §2.2 (rules (5), (S1)–(S4)).
//! All rules are transcribed verbatim from [KP15] §3.1–§3.2.
//!
//! # Rule labeling matches the paper verbatim
//!
//! | Paper | This module | Produces |
//! |-------|-------------|----------|
//! | eq (5)  | [`rule_5_reduce_edge`]     | `p̀_{30}` (edge degree-reduction) |
//! | (S1)    | [`rule_s1_star_propagation`] | `ṕ_{00..02}` |
//! | (S2) eq (7)  | [`rule_s2_p21`] / [`rule_s2_p12`] | `ṕ_{21}`, `ṕ_{12}` |
//! | (S2) eq (8)  | [`rule_s2_p22`]            | `ṕ_{22}` |
//! | (S3) eq (9)  | [`rule_s3`]                | `p̀_{31}`, `ṕ_{31}`, `ṕ_{32}` |
//! | (S4) eq (10), (11) | [`rule_s4`]          | `p̀_{i1}`, `ṕ_{i1}` (i=5,6); `p̀_{41}`, `ṕ_{41}` |
//!
//! The earlier revision grouped paper's (S2) eq (8) as a
//! separate "(S3)" and rolled paper's (S3) and (S4) into a
//! single "(S4)"; the current labels match [KP15] exactly so
//! cross-referencing with the paper reads 1:1.
//!
//! # Notation
//!
//! Pairs of adjacent bi-6 patches `p̀ := p_{k−1}` (previous
//! sector) and `ṕ := p_k` (current sector) share an edge. The
//! BB indices `ij` range over `0..=6` per direction. The shared
//! edge of `p̀` is at `j = 0` (`p̀_{i0}` for `i = 0..6`); the
//! star point `p_{00}` is identical in all sectors (the cap's
//! hub).
//!
//! # Coefficient conventions
//!
//! Rules mix affine combinations (coefficients summing to 1) and
//! polynomial-derived linear combinations (coefficients summing
//! to some `c`- / `τ`- dependent expression). In all cases the
//! assignment operates coordinate-wise on the `Point3` `x, y, z`
//! fields — this is the paper's convention. The resulting points
//! are valid BB-coefficients because the formulas are derived
//! from Bernstein polynomial identities that constrain the
//! patches' G² joining, not from affine-geometric constraints on
//! individual CPs.
//!
//! # API
//!
//! Each rule is a pure function taking the specific input
//! BB-CPs (by value, or bundled in a struct when the arg count
//! exceeds ten) and returning the output BB-CP(s). The caller
//! (I.10 cap assembly) orchestrates the dependency order.
//! No hidden state; no allocation on the hot path; deterministic.

use super::kp_polar_shared::weighted_sum;
use super::kp15_bi6_reparam::Bi6Reparams;
use crate::geometry::Point3;

/// Rule (5) — degree-reduction of the shared edge from 6 to 5.
/// Determines `p̀_{30}` from the other six edge CPs
/// `p̀_{i0}, i ∈ {0, 1, 2, 4, 5, 6}` [KP15 eq (5)]:
///
/// ```text
///   p̀_{30} = ¾·(p̀_{20} + p̀_{40})
///          − 3/10·(p̀_{10} + p̀_{50})
///          + 1/20·(p̀_{00} + p̀_{60})
/// ```
pub fn rule_5_reduce_edge(
    p00: Point3,
    p10: Point3,
    p20: Point3,
    p40: Point3,
    p50: Point3,
    p60: Point3,
) -> Point3 {
    weighted_sum(&[
        (0.05, p00),
        (-0.3, p10),
        (0.75, p20),
        (0.75, p40),
        (-0.3, p50),
        (0.05, p60),
    ])
}

/// Outputs of [`rule_s1_star_propagation`] — the six `ṕ_{ij}`
/// CPs with `i + j ≤ 2`.
#[derive(Debug, Clone, Copy)]
pub struct S1Outputs {
    /// `ṕ_{00}` = `p̀_{00}` (identity).
    pub p00: Point3,
    /// `ṕ_{10}` = `p̀_{10}` (identity).
    pub p10: Point3,
    /// `ṕ_{20}` = `p̀_{20}` (identity).
    pub p20: Point3,
    /// `ṕ_{01}` via `2c·p̀_{10} + 2(1−c)·p̀_{00} − p̀_{01}`.
    pub p01: Point3,
    /// `ṕ_{11}` via [KP15 eq (6)] line 3.
    pub p11: Point3,
    /// `ṕ_{02}` via [KP15 eq (6)] line 4 (the long formula).
    pub p02: Point3,
}

/// Rule (S1) — star-point BB propagation from the previous
/// sector to the current one [KP15 eq (6)] (math summary §2.2):
///
/// ```text
///   ṕ_{00} := p̀_{00},    ṕ_{10} := p̀_{10},    ṕ_{20} := p̀_{20}
///   ṕ_{01} := −p̀_{01} + 2c·p̀_{10} + 2(1 − c)·p̀_{00}
///   ṕ_{11} := −p̀_{11} + (5c/3)·p̀_{20} + 2(1 − c)·p̀_{10} + (c/3)·p̀_{00}
///   ṕ_{02} :=   p̀_{02}
///            − (24c/5)·p̀_{11}
///            + 4c²·p̀_{20}
///            + 4·(6τc/5 − 1)·p̀_{01}
///            + (4c/5)·(11 − 6(1 + τ)c)·p̀_{10}
///            + 4·(1 − c − 6τc/5 + ((6τ + 1)/5)c²)·p̀_{00}
/// ```
pub fn rule_s1_star_propagation(
    params: &Bi6Reparams,
    p_prev_00: Point3,
    p_prev_01: Point3,
    p_prev_02: Point3,
    p_prev_10: Point3,
    p_prev_11: Point3,
    p_prev_20: Point3,
) -> S1Outputs {
    let c = params.c;
    let tau = params.tau;

    let p01 = weighted_sum(&[
        (2.0 * (1.0 - c), p_prev_00),
        (2.0 * c, p_prev_10),
        (-1.0, p_prev_01),
    ]);

    let p11 = weighted_sum(&[
        (c / 3.0, p_prev_00),
        (2.0 * (1.0 - c), p_prev_10),
        (5.0 * c / 3.0, p_prev_20),
        (-1.0, p_prev_11),
    ]);

    // ṕ_{02} weights per [KP15 eq (6)]:
    //   p̀_{00}: 4·(1 − c − 6τc/5 + ((6τ + 1)/5)c²)
    //   p̀_{10}: (4c/5)·(11 − 6(1 + τ)c)
    //   p̀_{20}: 4c²
    //   p̀_{01}: 4·(6τc/5 − 1)
    //   p̀_{11}: −24c/5
    //   p̀_{02}: 1
    let w00 = 4.0 * (1.0 - c - 6.0 * tau * c / 5.0 + (6.0 * tau + 1.0) * c * c / 5.0);
    let w10 = (4.0 * c / 5.0) * (11.0 - 6.0 * (1.0 + tau) * c);
    let w20 = 4.0 * c * c;
    let w01 = 4.0 * (6.0 * tau * c / 5.0 - 1.0);
    let w11 = -24.0 * c / 5.0;
    let p02 = weighted_sum(&[
        (w00, p_prev_00),
        (w10, p_prev_10),
        (w20, p_prev_20),
        (w01, p_prev_01),
        (w11, p_prev_11),
        (1.0, p_prev_02),
    ]);

    S1Outputs {
        p00: p_prev_00,
        p10: p_prev_10,
        p20: p_prev_20,
        p01,
        p11,
        p02,
    }
}

/// Rule (S2) eq (7) — determines `ṕ_{21}`.
///
/// ```text
///   ṕ_{21} := −p̀_{21} + (2 − Σ µ̇_i)·p̀_{00}
///           + µ̇_1·p̀_{10} + µ̇_2·p̀_{20}
///           + µ̇_3·p̀_{40} + µ̇_4·p̀_{50} + µ̇_5·p̀_{60}
/// ```
#[allow(clippy::too_many_arguments)]
pub fn rule_s2_p21(
    mu_dot: [f64; 5],
    p_prev_00: Point3,
    p_prev_10: Point3,
    p_prev_20: Point3,
    p_prev_40: Point3,
    p_prev_50: Point3,
    p_prev_60: Point3,
    p_prev_21: Point3,
) -> Point3 {
    let sum_mu: f64 = mu_dot.iter().sum();
    weighted_sum(&[
        (2.0 - sum_mu, p_prev_00),
        (mu_dot[0], p_prev_10),
        (mu_dot[1], p_prev_20),
        (mu_dot[2], p_prev_40),
        (mu_dot[3], p_prev_50),
        (mu_dot[4], p_prev_60),
        (-1.0, p_prev_21),
    ])
}

/// Rule (S2) eq (7) — determines `ṕ_{12}`.
///
/// ```text
///   ṕ_{12} :=   p̀_{12} + (− Σ µ̈_i)·p̀_{00}
///           + µ̈_1·p̀_{01} + µ̈_2·p̀_{10} + µ̈_3·p̀_{11}
///           + µ̈_4·p̀_{20} + µ̈_5·p̀_{21}
///           + µ̈_6·p̀_{40} + µ̈_7·p̀_{50} + µ̈_8·p̀_{60}
/// ```
#[allow(clippy::too_many_arguments)]
pub fn rule_s2_p12(
    mu_ddot: [f64; 8],
    p_prev_00: Point3,
    p_prev_01: Point3,
    p_prev_10: Point3,
    p_prev_11: Point3,
    p_prev_20: Point3,
    p_prev_21: Point3,
    p_prev_40: Point3,
    p_prev_50: Point3,
    p_prev_60: Point3,
    p_prev_12: Point3,
) -> Point3 {
    let sum_mu: f64 = mu_ddot.iter().sum();
    weighted_sum(&[
        (-sum_mu, p_prev_00),
        (mu_ddot[0], p_prev_01),
        (mu_ddot[1], p_prev_10),
        (mu_ddot[2], p_prev_11),
        (mu_ddot[3], p_prev_20),
        (mu_ddot[4], p_prev_21),
        (mu_ddot[5], p_prev_40),
        (mu_ddot[6], p_prev_50),
        (mu_ddot[7], p_prev_60),
        (1.0, p_prev_12),
    ])
}

/// Rule (S2) eq (8) — determines `ṕ_{22}` via the `n × n`
/// circulant system.
///
/// ```text
///   ṕ_{22} :=   p̀_{22} + (− Σ_{i=1..8} ⃛µ_i)·p̀_{00}
///           + ⃛µ_1·p̀_{01} + ⃛µ_2·p̀_{10} + ⃛µ_3·p̀_{11}
///           + ⃛µ_4·p̀_{20} + ⃛µ_5·p̀_{21}
///           + ⃛µ_6·p̀_{40} + ⃛µ_7·p̀_{50} + ⃛µ_8·p̀_{60}
///           + ⃛µ_9·(p̀_{42} − ṕ_{42})
///           + ⃛µ_{10}·(p̀_{52} − ṕ_{52})
///           + ⃛µ_{11}·(p̀_{62} − ṕ_{62})
/// ```
///
/// The `p̀_{00}` weight sums **only indices 1..8** (⃛µ_{9,10,11}
/// multiply the diff_*2 terms, not `p̀_{00}`).
#[allow(clippy::too_many_arguments)]
pub fn rule_s2_p22(
    mu_tdot: [f64; 11],
    p_prev_00: Point3,
    p_prev_01: Point3,
    p_prev_10: Point3,
    p_prev_11: Point3,
    p_prev_20: Point3,
    p_prev_21: Point3,
    p_prev_40: Point3,
    p_prev_50: Point3,
    p_prev_60: Point3,
    p_prev_22: Point3,
    diff_42: Point3,
    diff_52: Point3,
    diff_62: Point3,
) -> Point3 {
    let sum_mu_1_to_8: f64 = mu_tdot.iter().take(8).sum();
    weighted_sum(&[
        (-sum_mu_1_to_8, p_prev_00),
        (mu_tdot[0], p_prev_01),
        (mu_tdot[1], p_prev_10),
        (mu_tdot[2], p_prev_11),
        (mu_tdot[3], p_prev_20),
        (mu_tdot[4], p_prev_21),
        (mu_tdot[5], p_prev_40),
        (mu_tdot[6], p_prev_50),
        (mu_tdot[7], p_prev_60),
        (mu_tdot[8], diff_42),
        (mu_tdot[9], diff_52),
        (mu_tdot[10], diff_62),
        (1.0, p_prev_22),
    ])
}

/// Inputs to [`rule_s4`]. Paper eqs (10) and (11) consume only
/// row-0 CPs at i ∈ {4, 5, 6} and row-2 differences at the same.
#[derive(Debug, Clone, Copy)]
pub struct S4Inputs {
    /// `p̀_{40}` (shared edge).
    pub p_prev_40: Point3,
    /// `p̀_{50}` (shared edge).
    pub p_prev_50: Point3,
    /// `p̀_{60}` (shared edge).
    pub p_prev_60: Point3,
    /// `p̀_{42}`, `ṕ_{42}` — row 2 at i=4.
    pub p_prev_42: Point3,
    /// See [`p_prev_42`](Self::p_prev_42).
    pub p_cur_42: Point3,
    /// `p̀_{52}`, `ṕ_{52}` — row 2 at i=5.
    pub p_prev_52: Point3,
    /// See [`p_prev_52`](Self::p_prev_52).
    pub p_cur_52: Point3,
    /// `p̀_{62}`, `ṕ_{62}` — row 2 at i=6.
    pub p_prev_62: Point3,
    /// See [`p_prev_62`](Self::p_prev_62).
    pub p_cur_62: Point3,
}

/// Outputs of [`rule_s4`] — the six remaining inner-CP
/// assignments at i ∈ {4, 5, 6}.
#[derive(Debug, Clone, Copy)]
pub struct S4Outputs {
    /// Local `p̀_{51}`; previous sector rotational `p_{15}`.
    pub p_prev_51: Point3,
    /// Local `ṕ_{51}`; current sector rotational `p_{51}`.
    pub p_cur_51: Point3,
    /// Local `p̀_{61}`; previous sector rotational `p_{16}`.
    pub p_prev_61: Point3,
    /// Local `ṕ_{61}`; current sector rotational `p_{61}`.
    pub p_cur_61: Point3,
    /// Local `p̀_{41}`; previous sector rotational `p_{14}`.
    pub p_prev_41: Point3,
    /// Local `ṕ_{41}`; current sector rotational `p_{41}`.
    pub p_cur_41: Point3,
}

/// Rule (S4) — simple differential assignments [KP15 eq (10),
/// (11)].
///
/// ```text
///   p̀_{i1} := p̀_{i0} + ¼·(p̀_{i2} − ṕ_{i2})         i ∈ {5, 6}
///   ṕ_{i1} := p̀_{i0} − ¼·(p̀_{i2} − ṕ_{i2})         i ∈ {5, 6}
///
///   w₁ := −c/(50τ),  w₂ := c/50,  w₃ := c/15
///   p̀_{41} := p̀_{40} + ¼·(p̀_{42} − ṕ_{42})
///                    + w₁·(p̀_{62} − ṕ_{62})
///                    + w₂·(p̀_{52} − ṕ_{52})
///                    + w₃·(p̀_{60} − p̀_{50})
///   ṕ_{41} := p̀_{40} − ¼·(p̀_{42} − ṕ_{42})
///                    − w₁·(p̀_{62} − ṕ_{62})
///                    − w₂·(p̀_{52} − ṕ_{52})
///                    + w₃·(p̀_{60} − p̀_{50})
/// ```
pub fn rule_s4(params: &Bi6Reparams, inp: &S4Inputs) -> S4Outputs {
    let c = params.c;
    let tau = params.tau;
    let w1 = -c / (50.0 * tau);
    let w2 = c / 50.0;
    let w3 = c / 15.0;

    let diff_42 = inp.p_prev_42 - inp.p_cur_42;
    let diff_52 = inp.p_prev_52 - inp.p_cur_52;
    let diff_62 = inp.p_prev_62 - inp.p_cur_62;
    let diff_6050 = inp.p_prev_60 - inp.p_prev_50;

    let p_prev_51 = inp.p_prev_50 + diff_52 * 0.25;
    let p_cur_51 = inp.p_prev_50 - diff_52 * 0.25;
    let p_prev_61 = inp.p_prev_60 + diff_62 * 0.25;
    let p_cur_61 = inp.p_prev_60 - diff_62 * 0.25;

    let p_prev_41 = inp.p_prev_40 + diff_42 * 0.25 + diff_62 * w1 + diff_52 * w2 + diff_6050 * w3;
    let p_cur_41 = inp.p_prev_40 - diff_42 * 0.25 - diff_62 * w1 - diff_52 * w2 + diff_6050 * w3;

    S4Outputs {
        p_prev_51,
        p_cur_51,
        p_prev_61,
        p_cur_61,
        p_prev_41,
        p_cur_41,
    }
}

#[cfg(test)]
#[path = "kp15_bi6_rules_tests.rs"]
mod tests;
