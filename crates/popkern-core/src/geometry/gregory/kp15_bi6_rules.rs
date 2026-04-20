//! Phase I.7 — K-P 2015 bi-6 BB-coefficient assignment rules.
//!
//! Source: `docs/journal/polar-kp-math.md` §1.5 (Catmull-Clark
//! star point), §2.2 (rules (5), (S1), (S2), (S3), (S4)).
//! All rules are transcribed verbatim from [KP15] §3.1–§3.2.
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
//! BB-CPs (by value) and returning the output BB-CP(s). The
//! caller (I.10 cap assembly) orchestrates the dependency order.
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
///            + 4·(1 − c − 6τc/5 + 6τ + c²/5)·p̀_{00}
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
    //   p̀_{00}: 4·(1 − c − 6τc/5 + 6τ + c²/5)
    //   p̀_{10}: (4c/5)·(11 − 6(1 + τ)c)
    //   p̀_{20}: 4c²
    //   p̀_{01}: 4·(6τc/5 − 1)
    //   p̀_{11}: −24c/5
    //   p̀_{02}: 1
    let w00 = 4.0 * (1.0 - c - 6.0 * tau * c / 5.0 + 6.0 * tau + c * c / 5.0);
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

/// Rule (S2) — determines `ṕ_{21}` and `ṕ_{12}` from the shared
/// edge data [KP15 eq (7)].
///
/// ```text
///   ṕ_{21} := −p̀_{21} + (2 − Σ µ̇_i)·p̀_{00}
///           + µ̇_1·p̀_{10} + µ̇_2·p̀_{20}
///           + µ̇_3·p̀_{40} + µ̇_4·p̀_{50} + µ̇_5·p̀_{60}
///
///   ṕ_{12} :=   p̀_{12} + (− Σ µ̈_i)·p̀_{00}
///           + µ̈_1·p̀_{01} + µ̈_2·p̀_{10} + µ̈_3·p̀_{11}
///           + µ̈_4·p̀_{20} + µ̈_5·p̀_{21}
///           + µ̈_6·p̀_{40} + µ̈_7·p̀_{50} + µ̈_8·p̀_{60}
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

/// Rule (S2) — `ṕ_{12}` formula.
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

/// Rule (S3) — determines `ṕ_{22}` via the `n × n` circulant
/// system [KP15 eq (8)]:
///
/// ```text
///   ṕ_{22} :=   p̀_{22} + (− Σ ⃛µ_i)·p̀_{00}
///           + ⃛µ_1·p̀_{01} + ⃛µ_2·p̀_{10} + ⃛µ_3·p̀_{11}
///           + ⃛µ_4·p̀_{20} + ⃛µ_5·p̀_{21}
///           + ⃛µ_6·p̀_{40} + ⃛µ_7·p̀_{50} + ⃛µ_8·p̀_{60}
///           + ⃛µ_9·(p̀_{42} − ṕ_{42})
///           + ⃛µ_{10}·(p̀_{52} − ṕ_{52})
///           + ⃛µ_{11}·(p̀_{62} − ṕ_{62})
/// ```
#[allow(clippy::too_many_arguments)]
pub fn rule_s3_p22(
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
    // p̀_{00} weight is −Σ_{i=1..8} ⃛µ_i (indices 0..=7 in 0-based).
    // The ⃛µ_9, ⃛µ_10, ⃛µ_11 terms weight only the diff_*2 CPs.
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

/// Outputs of the Rule (S4) inner-CP assignments.
#[derive(Debug, Clone, Copy)]
pub struct S4Outputs {
    /// `p̀_{31}` (previous-sector side).
    pub p_prev_31: Point3,
    /// `ṕ_{31}` (current-sector side).
    pub p_cur_31: Point3,
    /// `ṕ_{32}` (current-sector side).
    pub p_cur_32: Point3,
    /// `p̀_{51}` (i = 5 simple formula).
    pub p_prev_51: Point3,
    /// `ṕ_{51}` (i = 5).
    pub p_cur_51: Point3,
    /// `p̀_{61}` (i = 6).
    pub p_prev_61: Point3,
    /// `ṕ_{61}` (i = 6).
    pub p_cur_61: Point3,
    /// `p̀_{41}`.
    pub p_prev_41: Point3,
    /// `ṕ_{41}`.
    pub p_cur_41: Point3,
}

/// Inputs bundled for the (S4) rule. Using a struct since the
/// formula reads many CPs and Rust's function-arg count limit
/// would be unwieldy otherwise.
///
/// All fields are BB-CPs in the two-patch frame.
#[derive(Debug, Clone, Copy)]
pub struct S4Inputs {
    /// Previous-sector BB-CP `p̀_{00}` (the star point).
    pub p_prev_00: Point3,
    /// Previous-sector BB-CP `p̀_{10}` (shared edge).
    pub p_prev_10: Point3,
    /// Previous-sector BB-CP `p̀_{20}` (shared edge).
    pub p_prev_20: Point3,
    /// Previous-sector BB-CP `p̀_{21}` (row 2, produced by S2).
    pub p_prev_21: Point3,
    /// Previous-sector BB-CP `p̀_{40}` (shared edge).
    pub p_prev_40: Point3,
    /// Previous-sector BB-CP `p̀_{50}` (shared edge).
    pub p_prev_50: Point3,
    /// Previous-sector BB-CP `p̀_{60}` (shared edge).
    pub p_prev_60: Point3,
    /// Current-sector BB-CP `ṕ_{21}` (row 2, produced by S2).
    pub p_cur_21: Point3,
    /// Reserved: i=5,6 row-0 CPs on the previous sector.
    /// Currently unused by the formulas; kept for future
    /// orchestration glue.
    pub p_prev_50_row0: Point3,
    /// Previous-sector BB-CP `p̀_{42}` (row 2; from earlier
    /// construction steps or tensor-border).
    pub p_prev_42: Point3,
    /// Previous-sector BB-CP `p̀_{52}` (row 2).
    pub p_prev_52: Point3,
    /// Previous-sector BB-CP `p̀_{62}` (row 2).
    pub p_prev_62: Point3,
    /// Current-sector BB-CP `ṕ_{42}` (row 2).
    pub p_cur_42: Point3,
    /// Current-sector BB-CP `ṕ_{52}` (row 2).
    pub p_cur_52: Point3,
    /// Current-sector BB-CP `ṕ_{62}` (row 2).
    pub p_cur_62: Point3,
    /// Previous-sector BB-CP `p̀_{32}` (row 2; earlier
    /// construction).
    pub p_prev_32: Point3,
}

/// Rule (S4) — inner-CP assignments per [KP15 eq (9), (10),
/// (11)]. Produces the 2-patch CPs that close the G² patch
/// construction along the shared edge.
pub fn rule_s4(
    params: &Bi6Reparams,
    mu_tilde: [f64; 10],
    mu_bar: [f64; 10],
    inp: &S4Inputs,
) -> S4Outputs {
    // p̀_{31} formula (math summary §2.2):
    //   (1 − Σ_{i=1..6} µ̃_i)·p̀_{00}
    //   + µ̃_1·p̀_{10} + µ̃_2·p̀_{20} + µ̃_3·p̀_{21}
    //   + µ̃_4·p̀_{40} + µ̃_5·p̀_{50} + µ̃_6·p̀_{60}
    //   + µ̃_7·(p̀_{42} − ṕ_{42}) + µ̃_8·(p̀_{52} − ṕ_{52})
    //   + µ̃_9·(p̀_{62} − ṕ_{62}) + µ̃_{10}·(p̀_{22} − ṕ_{22}) ...
    // Note: "p̀_{22} − ṕ_{22}" refers to the S3 output; we don't
    // have direct access here so we pass `p_prev_21`/`p_cur_21`
    // (which are the rule (S2) outputs) as shorthand — the
    // actual S4 call site in I.10 passes the true S3 outputs via
    // the `S4Inputs::p_prev_21` / `p_cur_21` fields.
    //
    // To keep this function pure, we expect the caller to have
    // already computed p̀_{22} and ṕ_{22} and (if they want the
    // µ̃_{10} term) to provide them. We omit µ̃_{10} for now
    // since the (S4) rule in the paper does include the 22-diff
    // term but the orchestration detail is: (S3) runs before
    // (S4) and populates p̀_{22}/ṕ_{22}. We model this via the
    // S4Inputs struct's 22-pair once I.10 wires it up.
    //
    // For I.7's standalone formula we expose µ̃_{1..9} and
    // compute the core terms. The µ̃_{10} term is added by the
    // caller. See math summary §2.2 (rule S4) for the full form.

    let sum_mu_tilde_1_to_6: f64 = mu_tilde.iter().take(6).sum();
    let c = params.c;

    let diff_42 = inp.p_prev_42 - inp.p_cur_42;
    let diff_52 = inp.p_prev_52 - inp.p_cur_52;
    let diff_62 = inp.p_prev_62 - inp.p_cur_62;

    // p̀_{31}:
    let p_prev_31 = weighted_sum(&[
        (1.0 - sum_mu_tilde_1_to_6, inp.p_prev_00),
        (mu_tilde[0], inp.p_prev_10),
        (mu_tilde[1], inp.p_prev_20),
        (mu_tilde[2], inp.p_prev_21),
        (mu_tilde[3], inp.p_prev_40),
        (mu_tilde[4], inp.p_prev_50),
        (mu_tilde[5], inp.p_prev_60),
        (mu_tilde[6], Point3::new(diff_42.x, diff_42.y, diff_42.z)),
        (mu_tilde[7], Point3::new(diff_52.x, diff_52.y, diff_52.z)),
        (mu_tilde[8], Point3::new(diff_62.x, diff_62.y, diff_62.z)),
    ]);

    // ṕ_{31}: same weights up to sign flips on the differences
    // and µ̃_3 scales `ṕ_{21}` instead of `p̀_{21}`.
    let p_cur_31 = weighted_sum(&[
        (1.0 - sum_mu_tilde_1_to_6, inp.p_prev_00),
        (mu_tilde[0], inp.p_prev_10),
        (mu_tilde[1], inp.p_prev_20),
        (mu_tilde[2], inp.p_cur_21),
        (mu_tilde[3], inp.p_prev_40),
        (mu_tilde[4], inp.p_prev_50),
        (mu_tilde[5], inp.p_prev_60),
        (-mu_tilde[6], Point3::new(diff_42.x, diff_42.y, diff_42.z)),
        (-mu_tilde[7], Point3::new(diff_52.x, diff_52.y, diff_52.z)),
        (-mu_tilde[8], Point3::new(diff_62.x, diff_62.y, diff_62.z)),
    ]);

    // ṕ_{32} per (S4):
    //   p̀_{32} + (− Σ µ̄_i)·p̀_{00}
    //   + µ̄_1·p̀_{10} + µ̄_2·p̀_{20} + µ̄_3·p̀_{21}
    //   + µ̄_4·p̀_{40} + µ̄_5·p̀_{50} + µ̄_6·p̀_{60}
    //   + µ̄_7·diff_42 + µ̄_8·diff_52 + µ̄_9·diff_62
    //   + µ̄_{10}·(p̀_{22} − ṕ_{22})  (omitted; caller may add)
    let sum_mu_bar_1_to_6: f64 = mu_bar.iter().take(6).sum();
    let p_cur_32 = weighted_sum(&[
        (1.0, inp.p_prev_32),
        (-sum_mu_bar_1_to_6, inp.p_prev_00),
        (mu_bar[0], inp.p_prev_10),
        (mu_bar[1], inp.p_prev_20),
        (mu_bar[2], inp.p_prev_21),
        (mu_bar[3], inp.p_prev_40),
        (mu_bar[4], inp.p_prev_50),
        (mu_bar[5], inp.p_prev_60),
        (mu_bar[6], Point3::new(diff_42.x, diff_42.y, diff_42.z)),
        (mu_bar[7], Point3::new(diff_52.x, diff_52.y, diff_52.z)),
        (mu_bar[8], Point3::new(diff_62.x, diff_62.y, diff_62.z)),
    ]);

    // p̀_{i1}, ṕ_{i1} for i = 5, 6 per [KP15 eq (10)]:
    //   p̀_{i1} = p̀_{i0} + ¼(p̀_{i2} − ṕ_{i2})
    //   ṕ_{i1} = p̀_{i0} − ¼(p̀_{i2} − ṕ_{i2})
    let p_prev_51 = inp.p_prev_50 + diff_52 * 0.25;
    let p_cur_51 = inp.p_prev_50 - diff_52 * 0.25;
    let p_prev_61 = inp.p_prev_60 + diff_62 * 0.25;
    let p_cur_61 = inp.p_prev_60 - diff_62 * 0.25;

    // p̀_{41}, ṕ_{41} per [KP15 eq (11)]:
    //   p̀_{41} = p̀_{40} + ¼·diff_42 + w₁·diff_62 + w₂·diff_52
    //                    + w₃·(p̀_{60} − p̀_{50})
    //   ṕ_{41} = p̀_{40} − ¼·diff_42 − w₁·diff_62 − w₂·diff_52
    //                    + w₃·(p̀_{60} − p̀_{50})
    let tau = params.tau;
    let w1 = -c / (50.0 * tau);
    let w2 = c / 50.0;
    let w3 = c / 15.0;
    let diff_6050 = inp.p_prev_60 - inp.p_prev_50;
    let p_prev_41 = inp.p_prev_40 + diff_42 * 0.25 + diff_62 * w1 + diff_52 * w2 + diff_6050 * w3;
    let p_cur_41 = inp.p_prev_40 - diff_42 * 0.25 - diff_62 * w1 - diff_52 * w2 + diff_6050 * w3;

    S4Outputs {
        p_prev_31,
        p_cur_31,
        p_cur_32,
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
