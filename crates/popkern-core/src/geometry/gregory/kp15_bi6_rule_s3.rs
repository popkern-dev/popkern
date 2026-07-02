//! Phase I.7 — K-P 2015 bi-6 Rule (S3) per [KP15 eq (9)].
//!
//! Separated from [`kp15_bi6_rules`](super::kp15_bi6_rules)
//! because the Rule (S3) input struct is large (17 BB-CPs) and
//! the output formula has three 11-term linear combinations —
//! bundling the whole thing with the simpler rules would push
//! `kp15_bi6_rules.rs` over the 500-LoC cap.
//!
//! Produces the three per-sector-pair inner CPs
//! `p̀_{31}`, `ṕ_{31}`, `ṕ_{32}` using µ̃ and µ̄ (see
//! [`super::kp15_bi6_mu`]). Includes the µ̃_{10} / µ̄_{10}
//! terms that weight `(p̀_{22} − ṕ_{22})` — the paper lists
//! them as eq (9) lines 3 / 6.

use super::kp_polar_shared::weighted_sum;
use crate::geometry::Point3;

/// Inputs to [`rule_s3`]. Paper eq (9) consumes 17 BB-CPs per
/// sector-pair; bundling them avoids an unwieldy arg list.
#[derive(Debug, Clone, Copy)]
pub struct S3Inputs {
    /// `p̀_{00}` (the star point, shared by both sectors).
    pub p_prev_00: Point3,
    /// `p̀_{10}` (shared edge).
    pub p_prev_10: Point3,
    /// `p̀_{20}` (shared edge).
    pub p_prev_20: Point3,
    /// `p̀_{40}` (shared edge).
    pub p_prev_40: Point3,
    /// `p̀_{50}` (shared edge).
    pub p_prev_50: Point3,
    /// `p̀_{60}` (shared edge).
    pub p_prev_60: Point3,
    /// `p̀_{21}` (row 2; output of rule_s2_p21 at previous k).
    pub p_prev_21: Point3,
    /// `ṕ_{21}` (row 2; output of rule_s2_p21 at current k).
    pub p_cur_21: Point3,
    /// `p̀_{22}` (row 2; output of rule_s2_p22 at previous k).
    pub p_prev_22: Point3,
    /// `ṕ_{22}` (row 2; output of rule_s2_p22 at current k).
    pub p_cur_22: Point3,
    /// `p̀_{42}` — row 2 at i=4 (tensor-border).
    pub p_prev_42: Point3,
    /// `ṕ_{42}` — row 2 at i=4.
    pub p_cur_42: Point3,
    /// `p̀_{52}` — row 2 at i=5.
    pub p_prev_52: Point3,
    /// `ṕ_{52}` — row 2 at i=5.
    pub p_cur_52: Point3,
    /// `p̀_{62}` — row 2 at i=6.
    pub p_prev_62: Point3,
    /// `ṕ_{62}` — row 2 at i=6.
    pub p_cur_62: Point3,
    /// Local `p̀_{32}`. Under Fig. 2's swap between the previous
    /// sector's rotational indexing and the sector-symmetric
    /// pair-frame, this is the previous sector's `p_{23}` free
    /// parameter from KP15 §3.6.
    pub p_prev_32: Point3,
}

/// Outputs of [`rule_s3`].
#[derive(Debug, Clone, Copy)]
pub struct S3Outputs {
    /// `p̀_{31}` (previous-sector side).
    pub p_prev_31: Point3,
    /// `ṕ_{31}` (current-sector side).
    pub p_cur_31: Point3,
    /// `ṕ_{32}` (current-sector side).
    pub p_cur_32: Point3,
}

/// Rule (S3) eq (9) — inner-CP assignments `p̀_{31}`, `ṕ_{31}`,
/// `ṕ_{32}` using µ̃ and µ̄ [KP15 eq (9)]. Includes the full
/// µ̃_{10}/µ̄_{10} terms that weight `(p̀_{22} − ṕ_{22})`.
///
/// ```text
///   p̀_{31} := (1 − Σ_{i=1..6} µ̃_i)·p̀_{00}
///           + µ̃_1·p̀_{10} + µ̃_2·p̀_{20} + µ̃_3·p̀_{21}
///           + µ̃_4·p̀_{40} + µ̃_5·p̀_{50} + µ̃_6·p̀_{60}
///           + µ̃_7·(p̀_{42} − ṕ_{42}) + µ̃_8·(p̀_{52} − ṕ_{52})
///           + µ̃_9·(p̀_{62} − ṕ_{62}) + µ̃_{10}·(p̀_{22} − ṕ_{22})
///
///   ṕ_{31} := (1 − Σ_{i=1..6} µ̃_i)·p̀_{00}
///           + µ̃_1·p̀_{10} + µ̃_2·p̀_{20} + µ̃_3·ṕ_{21}
///           + µ̃_4·p̀_{40} + µ̃_5·p̀_{50} + µ̃_6·p̀_{60}
///           − µ̃_7·(p̀_{42} − ṕ_{42}) − µ̃_8·(p̀_{52} − ṕ_{52})
///           − µ̃_9·(p̀_{62} − ṕ_{62}) − µ̃_{10}·(p̀_{22} − ṕ_{22})
///
///   ṕ_{32} := p̀_{32} + (− Σ_{i=1..6} µ̄_i)·p̀_{00}
///           + µ̄_1·p̀_{10} + µ̄_2·p̀_{20} + µ̄_3·p̀_{21}
///           + µ̄_4·p̀_{40} + µ̄_5·p̀_{50} + µ̄_6·p̀_{60}
///           + µ̄_7·(p̀_{42} − ṕ_{42}) + µ̄_8·(p̀_{52} − ṕ_{52})
///           + µ̄_9·(p̀_{62} − ṕ_{62}) + µ̄_{10}·(p̀_{22} − ṕ_{22})
/// ```
///
/// µ̃ and µ̄ come from I.6 (`kp15_bi6_mu::mu_tilde`, `mu_bar`).
pub fn rule_s3(mu_tilde: [f64; 10], mu_bar: [f64; 10], inp: &S3Inputs) -> S3Outputs {
    let sum_mu_tilde_1_to_6: f64 = mu_tilde.iter().take(6).sum();
    let sum_mu_bar_1_to_6: f64 = mu_bar.iter().take(6).sum();

    // Differences flow into the sum as Point3 for the
    // coordinate-wise `weighted_sum` convention.
    let diff_42_p = {
        let v = inp.p_prev_42 - inp.p_cur_42;
        Point3::new(v.x, v.y, v.z)
    };
    let diff_52_p = {
        let v = inp.p_prev_52 - inp.p_cur_52;
        Point3::new(v.x, v.y, v.z)
    };
    let diff_62_p = {
        let v = inp.p_prev_62 - inp.p_cur_62;
        Point3::new(v.x, v.y, v.z)
    };
    let diff_22_p = {
        let v = inp.p_prev_22 - inp.p_cur_22;
        Point3::new(v.x, v.y, v.z)
    };

    let p_prev_31 = weighted_sum(&[
        (1.0 - sum_mu_tilde_1_to_6, inp.p_prev_00),
        (mu_tilde[0], inp.p_prev_10),
        (mu_tilde[1], inp.p_prev_20),
        (mu_tilde[2], inp.p_prev_21),
        (mu_tilde[3], inp.p_prev_40),
        (mu_tilde[4], inp.p_prev_50),
        (mu_tilde[5], inp.p_prev_60),
        (mu_tilde[6], diff_42_p),
        (mu_tilde[7], diff_52_p),
        (mu_tilde[8], diff_62_p),
        (mu_tilde[9], diff_22_p),
    ]);

    let p_cur_31 = weighted_sum(&[
        (1.0 - sum_mu_tilde_1_to_6, inp.p_prev_00),
        (mu_tilde[0], inp.p_prev_10),
        (mu_tilde[1], inp.p_prev_20),
        (mu_tilde[2], inp.p_cur_21),
        (mu_tilde[3], inp.p_prev_40),
        (mu_tilde[4], inp.p_prev_50),
        (mu_tilde[5], inp.p_prev_60),
        (-mu_tilde[6], diff_42_p),
        (-mu_tilde[7], diff_52_p),
        (-mu_tilde[8], diff_62_p),
        (-mu_tilde[9], diff_22_p),
    ]);

    let p_cur_32 = weighted_sum(&[
        (1.0, inp.p_prev_32),
        (-sum_mu_bar_1_to_6, inp.p_prev_00),
        (mu_bar[0], inp.p_prev_10),
        (mu_bar[1], inp.p_prev_20),
        (mu_bar[2], inp.p_prev_21),
        (mu_bar[3], inp.p_prev_40),
        (mu_bar[4], inp.p_prev_50),
        (mu_bar[5], inp.p_prev_60),
        (mu_bar[6], diff_42_p),
        (mu_bar[7], diff_52_p),
        (mu_bar[8], diff_62_p),
        (mu_bar[9], diff_22_p),
    ]);

    S3Outputs {
        p_prev_31,
        p_cur_31,
        p_cur_32,
    }
}

#[cfg(test)]
#[path = "kp15_bi6_rule_s3_tests.rs"]
mod tests;
