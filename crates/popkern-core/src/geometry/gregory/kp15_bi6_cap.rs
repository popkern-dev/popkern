//! Phase I.10a — bi-6 polar cap scaffolding.
//!
//! Assembles N bi-6 patches per the K-P 2015 polar construction
//! ([KP15 §3.6] algorithm) given:
//! 1. The Catmull-Clark cap-net inputs (star point + two outer
//!    rings).
//! 2. Per-sector boundary data at rows 0 / 1 / 2 (the tensor-
//!    border already reparameterised to bi-6 by the caller).
//! 3. Free-parameter values for `p⁰_{22}` (pinned free) and
//!    optionally `p⁰_{21}` at `n ∈ {3, 6}`, plus per-sector
//!    interior CPs at j = 3.
//!
//! # Scope
//!
//! This is the **scaffolding** commit for I.10 — runs all the
//! I.1–I.9 primitives end-to-end and produces bi-9 grids via
//! degree elevation. Does **not** yet perform F_5 minimisation
//! over the free parameters (that's I.10b) or wire into the
//! fillet crate (I.10c).
//!
//! # Pipeline
//!
//! 1. Catmull-Clark star point = `p^k_{00}` for all sectors
//!    (shared).
//! 2. Circulant (S2) eq (7): solve 2n×2n for
//!    `(p^k_{21}, p^k_{12})`.
//! 3. Circulant (S2) eq (8): solve n×n for `p^k_{22}` (skipped
//!    at `n = 3`; paper §3.5 marks all `p^k_{22}` free).
//! 4. Per sector, apply (S3) to get `p^k_{31}, p^k_{32}` and
//!    apply `p̀_{31}` to the previous sector.
//! 5. Per sector, apply (S4) to get `p^k_{51}, p^k_{61},
//!    p^k_{41}` and analogues on the previous sector.
//! 6. Populate remaining "interior" BB-CPs (rows 3..6 at
//!    j > 2) from the caller-supplied `interior_filler`; zero is
//!    a safe default for smoke testing.
//! 7. Degree-elevate each bi-6 grid to bi-9 via
//!    [`elevate_bi6_grid_to_bi9`](super::elevate_bi6_grid_to_bi9).

use super::kp_degree_elevate::elevate_bi6_grid_to_bi9;
use super::kp_polar_shared::catmull_clark_star_point;
use super::kp15_bi6_circulant::{
    CirculantSectorRhs, solve_bi6_circulant_s2_eq7, solve_bi6_circulant_s2_eq8,
};
use super::kp15_bi6_mu::{mu_bar, mu_ddot, mu_dot, mu_tdot, mu_tilde};
use super::kp15_bi6_reparam::Bi6Reparams;
use super::kp15_bi6_rule_s3::{S3Inputs, rule_s3};
use super::kp15_bi6_rules::{S4Inputs, rule_s4};
use crate::geometry::Point3;

/// Cap input. All caller-supplied data; no optimisation yet.
#[derive(Debug, Clone)]
pub struct Bi6CapInput {
    /// Valence; must be ≥ 3.
    pub n: usize,
    /// Per-valence K-P 2015 parameters (τ, c, d₂).
    pub params: Bi6Reparams,

    /// Catmull-Clark inputs: star point of the cap-net.
    pub c7: Point3,
    /// Outer-ring nodes (inner — "c_5" in the paper's Fig. 1a),
    /// length `n`.
    pub c5s: Vec<Point3>,
    /// Outer-ring nodes (outer — "c_6"), length `n`.
    pub c6s: Vec<Point3>,

    /// Per-sector circulant RHS data (the row-0 and row-1
    /// BB-CPs the circulant solvers need). Length `n`.
    pub sectors: Vec<CirculantSectorRhs>,

    /// Free parameter `p⁰_{22}` — pinned for the eq (8) solve.
    /// Paper §2.5 step 4 places this in the free-parameter set
    /// `P`; this scaffold takes it as input. F_5 min in I.10b
    /// will choose the optimal value.
    pub p0_22: Point3,

    /// Required at `n ∈ {3, 6}` (eq (7) system rank-deficient);
    /// optional otherwise.
    pub p0_21_pinned: Option<Point3>,

    /// Optional per-sector fillers for "deep interior" BB-CPs
    /// (rows 3..6 at j > 2). `interior_filler[k][i][j]` is the
    /// value to use for `p^k_{i, j}` at `i, j ∈ 3..=6`. `None`
    /// leaves the cell at the cap star.
    ///
    /// Paper §2.5 step 4 places `p^k_{23}` and `p^k_{33}` in
    /// the free parameter set; others are also effectively
    /// free (the F_5 min chooses them). This scaffold delegates
    /// to the caller; I.10b will supply the optimised values.
    pub interior_filler: Option<Vec<[[Option<Point3>; 7]; 7]>>,
}

/// Cap output — N bi-6 patches (49 CPs each) and their bi-9
/// degree-elevations (100 CPs each).
#[derive(Debug, Clone)]
pub struct Bi6CapOutput {
    /// Star point (common `p^k_{00}` across all sectors).
    pub star_point: Point3,
    /// Per-sector bi-6 CP grids `[i, j] ∈ 0..7 × 0..7`.
    pub patches_bi6: Vec<[[Point3; 7]; 7]>,
    /// Per-sector bi-9 CP grids (degree-elevated for output
    /// uniformity with Phase H rectangular path).
    pub patches_bi9: Vec<[[Point3; 10]; 10]>,
}

/// Assemble the bi-6 cap end-to-end. Does **not** optimise free
/// parameters; the caller is responsible for supplying them.
///
/// Returns an error if the circulant systems are singular for
/// the given inputs (shouldn't happen for well-formed inputs
/// with valid `p0_21_pinned` at `n ∈ {3, 6}`).
pub fn assemble_bi6_polar_cap(input: &Bi6CapInput) -> Result<Bi6CapOutput, String> {
    let n = input.n;
    if n < 3 {
        return Err(format!("bi-6 cap needs n ≥ 3, got {n}"));
    }
    input.params.validate()?;
    if input.sectors.len() != n {
        return Err(format!("sectors length {} ≠ n = {n}", input.sectors.len()));
    }
    if input.c5s.len() != n || input.c6s.len() != n {
        return Err(format!(
            "c5s / c6s length mismatch (got {}, {}); expected {n}",
            input.c5s.len(),
            input.c6s.len()
        ));
    }

    // 1. Star point.
    let star = catmull_clark_star_point(n, input.c7, &input.c5s, &input.c6s);

    // 2. Circulant (S2) eq (7): (p_21, p_12) ring.
    let mu_d = mu_dot(&input.params);
    let mu_dd = mu_ddot(&input.params);
    let s2_eq7 = solve_bi6_circulant_s2_eq7(
        &input.params,
        mu_d,
        mu_dd,
        &input.sectors,
        input.p0_21_pinned,
    )?;
    let p_21 = s2_eq7.p_21;
    let p_12 = s2_eq7.p_12;

    // 3. Circulant (S2) eq (8): p_22 ring. Skipped at n = 3.
    let mu_td = mu_tdot(&input.params);
    let p_22: Vec<Point3> = if n == 3 {
        // Paper §3.5: all three p^k_{22} are free. Scaffold
        // uses the caller's `p0_22` for k = 0 and zero for
        // k = 1, 2. I.10b will hook this up to the F_5 min.
        vec![
            input.p0_22,
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
        ]
    } else {
        solve_bi6_circulant_s2_eq8(&input.params, mu_td, &input.sectors, &p_12, input.p0_22)?
    };

    // 4. Rule (S3): p^k_{31}, p^{k+1}_{31}, p^{k+1}_{32} for each
    //    sector-pair (k, k+1). Collect the per-sector p_31, p_32
    //    assignments.
    let mt = mu_tilde(&input.params);
    let mb = mu_bar(&input.params);
    let mut p_31: Vec<Point3> = vec![Point3::new(0.0, 0.0, 0.0); n];
    let mut p_32: Vec<Point3> = vec![Point3::new(0.0, 0.0, 0.0); n];
    // For (S3), we need `p̀_{32}` = p_32 of the previous sector.
    // Start with a default value (will be overwritten by F_5 min
    // in I.10b, or by interior_filler here if supplied).
    let p_prev_32_default = Point3::new(0.0, 0.0, 0.0);
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        let s_prev = &input.sectors[km1];
        let s_cur = &input.sectors[k];
        let s3_in = S3Inputs {
            p_prev_00: s_prev.p_00,
            p_prev_10: s_prev.p_10,
            p_prev_20: s_prev.p_20,
            p_prev_40: s_prev.p_40,
            p_prev_50: s_prev.p_50,
            p_prev_60: s_prev.p_60,
            p_prev_21: p_21[km1],
            p_cur_21: p_21[k],
            p_prev_22: p_22[km1],
            p_cur_22: p_22[k],
            p_prev_42: s_prev.p_42,
            p_cur_42: s_cur.p_42,
            p_prev_52: s_prev.p_52,
            p_cur_52: s_cur.p_52,
            p_prev_62: s_prev.p_62,
            p_cur_62: s_cur.p_62,
            p_prev_32: p_prev_32_default,
        };
        let s3_out = rule_s3(mt, mb, &s3_in);
        p_31[km1] = s3_out.p_prev_31;
        p_31[k] = s3_out.p_cur_31;
        p_32[k] = s3_out.p_cur_32;
    }

    // 5. Rule (S4): p^k_{51}, p^k_{61}, p^k_{41} for each pair.
    let mut p_51: Vec<Point3> = vec![Point3::new(0.0, 0.0, 0.0); n];
    let mut p_61: Vec<Point3> = vec![Point3::new(0.0, 0.0, 0.0); n];
    let mut p_41: Vec<Point3> = vec![Point3::new(0.0, 0.0, 0.0); n];
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        let s_prev = &input.sectors[km1];
        let s_cur = &input.sectors[k];
        let s4_in = S4Inputs {
            p_prev_40: s_prev.p_40,
            p_prev_50: s_prev.p_50,
            p_prev_60: s_prev.p_60,
            p_prev_42: s_prev.p_42,
            p_cur_42: s_cur.p_42,
            p_prev_52: s_prev.p_52,
            p_cur_52: s_cur.p_52,
            p_prev_62: s_prev.p_62,
            p_cur_62: s_cur.p_62,
        };
        let s4_out = rule_s4(&input.params, &s4_in);
        p_51[km1] = s4_out.p_prev_51;
        p_51[k] = s4_out.p_cur_51;
        p_61[km1] = s4_out.p_prev_61;
        p_61[k] = s4_out.p_cur_61;
        p_41[km1] = s4_out.p_prev_41;
        p_41[k] = s4_out.p_cur_41;
    }

    // 6. Populate each per-sector bi-6 grid.
    //
    // Unpopulated CPs (in particular, deep interior at rows 3..7
    // × cols 3..7 not touched by interior_filler or the S3/S4
    // rules) default to the **star point** rather than origin.
    // Origin is a geometrically-unrelated 3D location; putting
    // unset CPs there causes the bi-6 surface to bend through
    // origin and produces folds at (u=1, v=1) on every sub-patch.
    // Star-point default at least keeps the surface near the cap's
    // natural centre. Fixes the fold diagnosed via
    // `POPKERN_KP_POLAR_DEEP_TRACE=1` (see
    // `f5-tuning-plan.md` finding 2026-04-21). Option B' of the
    // plan (explicit deep-interior derivation / Coons fill) can
    // later replace this pragmatic default with a paper-correct
    // construction.
    let mut patches_bi6: Vec<[[Point3; 7]; 7]> = Vec::with_capacity(n);
    for k in 0..n {
        let s = &input.sectors[k];
        let mut grid = [[star; 7]; 7];
        // Row 0 (shared edge, from caller).
        grid[1][0] = s.p_10;
        grid[2][0] = s.p_20;
        grid[4][0] = s.p_40;
        grid[5][0] = s.p_50;
        grid[6][0] = s.p_60;
        // grid[3][0] via rule (5) — degree-reduction of the edge.
        grid[3][0] = super::kp15_bi6_rules::rule_5_reduce_edge(
            s.p_00, s.p_10, s.p_20, s.p_40, s.p_50, s.p_60,
        );
        // Row 1 CPs from (S1) propagation + circulant + rules.
        grid[0][1] = s.p_01;
        grid[1][1] = s.p_11;
        grid[2][1] = p_21[k];
        grid[3][1] = p_31[k];
        grid[4][1] = p_41[k];
        grid[5][1] = p_51[k];
        grid[6][1] = p_61[k];
        // Row 2.
        grid[0][2] = s.p_02;
        grid[1][2] = p_12[k];
        grid[2][2] = p_22[k];
        grid[3][2] = p_32[k];
        grid[4][2] = s.p_42;
        grid[5][2] = s.p_52;
        grid[6][2] = s.p_62;
        // Rows 3..6: deep interior. Use caller-supplied filler
        // or leave at the star default.
        if let Some(filler) = input.interior_filler.as_ref()
            && let Some(per_sector) = filler.get(k)
        {
            for i in 0..7 {
                for j in 3..7 {
                    if let Some(cp) = per_sector[i][j] {
                        grid[i][j] = cp;
                    }
                }
            }
        }
        patches_bi6.push(grid);
    }

    // 7. Degree-elevate bi-6 → bi-9.
    let patches_bi9: Vec<[[Point3; 10]; 10]> = patches_bi6
        .iter()
        .map(|g| elevate_bi6_grid_to_bi9(*g))
        .collect();

    Ok(Bi6CapOutput {
        star_point: star,
        patches_bi6,
        patches_bi9,
    })
}

/// Phase E.3 production-side entry point — assemble a bi-6 polar
/// cap via the [KP15] constraint-matrix table with per-coordinate
/// F_κ minimisation over the deep-interior BB-coefs.
///
/// Equivalent input/output to [`assemble_bi6_polar_cap`], but
/// instead of populating rows 3..6 from `input.interior_filler`
/// (defaulting to the star), this path:
/// 1. Builds the cross-sector identification-augmented constraint
///    matrix `M(n, c, τ)` via
///    [`super::kp15_bi6_local_constraint::Kp15Bi6GeneratingTable::build_with_identifications`].
/// 2. Seeds the free-column subset from `input.sectors` and the
///    Catmull-Clark star at `(0, 0)` per sector.
/// 3. Minimises F_κ (κ from `kappa`) over the deep-interior free
///    columns via three independent SPD solves, one per
///    coordinate (paper §3.7: "works for each coordinate
///    separately").
/// 4. Applies the table to derive every other BB-coef from the
///    solved free values.
/// 5. Degree-elevates each per-sector bi-6 grid to bi-9.
///
/// Compared with `assemble_bi6_polar_cap`'s star-default chain,
/// this path delivers:
/// - Machine-zero residuals on [KP15 eq (2)] (G¹) and
///   [KP15 eq (3)] (G²) at every cross-sector edge.
/// - F_κ values 3.5–5.5 orders of magnitude smaller than the
///   zero-default cap (see
///   `phase_e3_table_driven_f5_satisfies_eq2_eq3_and_beats_existing_f5`
///   in the local constraint module's tests).
///
/// The `interior_filler` field of `Bi6CapInput` is **ignored** on
/// this path. The `p0_22` and `p0_21_pinned` free parameters are
/// also currently ignored — the F_κ solve picks them along with
/// the deep interior. Caller-side overrides for those will land
/// when this path becomes the default.
pub fn assemble_bi6_polar_cap_via_table(
    input: &Bi6CapInput,
    kappa: usize,
) -> Result<Bi6CapOutput, String> {
    let n = input.n;
    if n < 3 {
        return Err(format!("bi-6 cap needs n ≥ 3, got {n}"));
    }
    input.params.validate()?;
    if input.sectors.len() != n {
        return Err(format!("sectors length {} ≠ n = {n}", input.sectors.len()));
    }
    if input.c5s.len() != n || input.c6s.len() != n {
        return Err(format!(
            "c5s / c6s length mismatch (got {}, {}); expected {n}",
            input.c5s.len(),
            input.c6s.len()
        ));
    }
    let star = catmull_clark_star_point(n, input.c7, &input.c5s, &input.c6s);
    let table =
        super::kp15_bi6_local_constraint::Kp15Bi6GeneratingTable::build_with_identifications(
            n,
            input.params,
        );
    let q_local = super::kp15_bi6_f5::f_kappa_quadratic_form(kappa);
    let (patches_bi6, _f_values) =
        table.populate_cap_grids_from_input_with_f_kappa(&input.sectors, star, &q_local)?;
    let patches_bi9: Vec<[[Point3; 10]; 10]> = patches_bi6
        .iter()
        .map(|g| elevate_bi6_grid_to_bi9(*g))
        .collect();
    Ok(Bi6CapOutput {
        star_point: star,
        patches_bi6,
        patches_bi9,
    })
}

#[cfg(test)]
#[path = "kp15_bi6_cap_tests.rs"]
mod tests;
