#![allow(dead_code)]
//! the polar-cap campaignI F_5 tuning option (F) — Levenberg-Marquardt refinement
//! on the polar bi-6 free parameter set.
//!
//! When the F_κ fallback ladder (option C) exhausts without
//! finding a non-folding cap, L-M is the next escape hatch. It
//! uses the same infrastructure that got the bi-9 cap campaign's binonic path
//! from G¹ ≈ 1.06e-1 (naive gradient descent floor) down to
//! 1.12e-2 on the `cube_extreme` fixture.
//!
//! # Why L-M helps where F_κ doesn't
//!
//! F_κ is a quadratic smoothness functional. Its SPD minimum is
//! unique given fixed boundary constraints. If every κ ∈ {1, 3,
//! 4, 5} produces a folded cap (as observed on sharp-apex
//! pentagons), the boundary constraints are forcing folds
//! regardless of the solver's freedom — **within the convex F_κ
//! framework**.
//!
//! L-M uses a **different objective**: per-sample G¹ coplanarity
//! at shared spokes, with self-intersection as a hard-constraint
//! gate. Starting from a (possibly folded) F_κ warm start, L-M
//! steps away from smoothness-optimal configurations toward
//! G¹-optimal ones, rejecting steps that fold. If a nearby non-
//! folded configuration exists with acceptable G¹, L-M finds it.
//!
//! If L-M also cannot find one, the fold is genuinely structural
//! (bi-6 can't hold that curvature) and K-P 2013 bi-5 escalation
//! is the paper's prescribed answer.

use super::kp_lm_solver::{LmConfig, LmEval, run_lm};
use super::kp_polar_deep_fill::boundary_owned_deep_fill_input;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{Bi6CapFreeParams, Bi6CapInput, assemble_bi6_polar_cap};

const FD_EPS: f64 = 1e-5;
const LAMBDA_INIT: f64 = 1e-3;
const LAMBDA_MAX: f64 = 1e10;
const MAX_ITER: usize = 40;
const MAX_FINITE_DIFFERENCE_EVALUATIONS: usize = 2_500;
const TARGET_RESIDUAL: f64 = 1e-8;
const N_SAMPLES_PER_SPOKE: usize = 10;
const SELF_INT_TOL: f64 = 1e-6;

/// Option D — fold-penalty threshold. Each sub-patch samples
/// `|∂u × ∂v|` at a grid of interior parameters; if any sample
/// falls below this threshold, an additional quadratic penalty
/// residual is emitted (shape `(threshold − sample)²`) to push
/// L-M away from near-fold regions. A non-zero penalty doesn't
/// forbid the configuration (the hard fold-safety check still
/// governs accept/reject), but it shapes the L-M descent
/// direction.
const FOLD_PENALTY_THRESHOLD: f64 = 0.01;

/// Weight multiplier for the fold-penalty residual relative to
/// the G¹ coplanarity residual. Chosen large enough to dominate
/// when fold-proximity is severe but small enough that L-M
/// doesn't trade G¹ quality for "fold-distance" once well clear
/// of the fold boundary.
const FOLD_PENALTY_WEIGHT: f64 = 10.0;

/// Sample count per patch for fold-proximity evaluation.
const FOLD_PENALTY_SAMPLES: usize = 5;
const EXTRA_DEEP_CELLS: [(usize, usize); 20] = [
    (1, 3),
    (2, 3),
    (3, 3),
    (4, 3),
    (5, 3),
    (1, 4),
    (2, 4),
    (3, 4),
    (4, 4),
    (5, 4),
    (1, 5),
    (2, 5),
    (3, 5),
    (4, 5),
    (5, 5),
    (1, 6),
    (2, 6),
    (3, 6),
    (4, 6),
    (5, 6),
];

pub(super) struct PolarLmReport {
    pub patches: Option<Vec<[[Point3; 10]; 10]>>,
    pub rejection_reason: Option<&'static str>,
}

enum SingleLmOutcome {
    Accepted(Vec<[[Point3; 10]; 10]>),
    InfeasibleStart,
    InfeasibleStartRejected,
    RejectedAfterFeasibleStart,
}

/// L-M evaluator wrapping a `Bi6CapInput`. Each residual is the
/// coplanarity of the surface normals at one sample on one shared
/// spoke between adjacent bi-9 patches in the polar ring.
struct PolarLmEval<'a> {
    base: &'a Bi6CapInput,
}

struct PolarLmVars {
    free: Bi6CapFreeParams,
    extra_deep: Vec<[Point3; EXTRA_DEEP_CELLS.len()]>,
}

impl PolarLmEval<'_> {
    fn n(&self) -> usize {
        self.base.n
    }

    /// Unpack a flat parameter vector into polar L-M variables.
    /// Layout matches `Bi6CapFreeParams` probe convention: x/y/z
    /// for `p0_01, p0_02, p0_11, p0_22`, then per-sector
    /// `p_23[k], p_33[k]` triplets, then selected deep-interior
    /// cells per sector. The selected cells deliberately exclude
    /// the exact interface edge, which remains owned by the
    /// boundary construction.
    fn unpack(&self, x: &[f64]) -> Option<PolarLmVars> {
        let n = self.n();
        let free_len = 2 * n + 4;
        let expected = 3 * (free_len + n * EXTRA_DEEP_CELLS.len());
        if x.len() != expected {
            return None;
        }
        let pt = |i: usize| Point3::new(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
        let p0_01 = pt(0);
        let p0_02 = pt(1);
        let p0_11 = pt(2);
        let p0_22 = pt(3);
        let mut p_23 = Vec::with_capacity(n);
        let mut p_33 = Vec::with_capacity(n);
        for k in 0..n {
            p_23.push(pt(4 + 2 * k));
            p_33.push(pt(4 + 2 * k + 1));
        }
        let free = Bi6CapFreeParams {
            p0_01,
            p0_02,
            p0_11,
            p0_22,
            p_23,
            p_33,
        };
        let mut extra_deep = Vec::with_capacity(n);
        let mut offset = free_len;
        for _ in 0..n {
            let mut cells = [Point3::new(0.0, 0.0, 0.0); EXTRA_DEEP_CELLS.len()];
            for cell in &mut cells {
                *cell = pt(offset);
                offset += 1;
            }
            extra_deep.push(cells);
        }
        Some(PolarLmVars { free, extra_deep })
    }

    /// Apply free params + (S1) propagation (mirroring
    /// `kp15_bi6_cap_f5_solve::apply_free_params` for consistency)
    /// and assemble the cap.
    fn build(&self, vars: &PolarLmVars) -> Option<Vec<[[Point3; 10]; 10]>> {
        let mut input = boundary_owned_deep_fill_input(self.base, &vars.free, 1.0)?;
        if let Some(filler) = input.interior_filler.as_mut() {
            for (k, cells) in vars.extra_deep.iter().enumerate() {
                for (idx, &(i, j)) in EXTRA_DEEP_CELLS.iter().enumerate() {
                    filler[k][i][j] = Some(cells[idx]);
                }
            }
        }
        assemble_bi6_polar_cap(&input)
            .ok()
            .map(|out| out.patches_bi9)
    }

    /// Per-sample G¹ coplanarity residuals for the ring, plus
    /// option D fold-penalty residuals. Layout:
    /// `[coplan_0, ..., coplan_{10n-1}, fold_0, ..., fold_{5n-1}]`.
    ///
    /// Each coplanarity term is `|(xtan_a × xtan_b) · t_spoke| /
    /// (|xtan_a||xtan_b||t_spoke|)` at one spoke sample. Zero = G¹.
    ///
    /// Each fold-penalty term is `max(0, threshold - |∂u × ∂v|)`
    /// at one patch interior sample, scaled by
    /// `FOLD_PENALTY_WEIGHT`. Non-zero when a sample approaches
    /// a fold boundary; L-M pushes away from such configurations.
    fn compute_residuals(&self, patches: &[[[Point3; 10]; 10]]) -> Vec<f64> {
        use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
        let n = patches.len();
        let ts: Vec<f64> = (0..N_SAMPLES_PER_SPOKE)
            .map(|k| 0.05 + 0.90 * (k as f64) / ((N_SAMPLES_PER_SPOKE - 1) as f64))
            .collect();
        let mut out = Vec::with_capacity(
            n * N_SAMPLES_PER_SPOKE + n * FOLD_PENALTY_SAMPLES * FOLD_PENALTY_SAMPLES,
        );
        let wrapped: Vec<GregoryPatch> = patches
            .iter()
            .map(|g| {
                GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid: *g }))
            })
            .collect();
        // G¹ coplanarity residuals.
        for k in 0..n {
            let a = &wrapped[k];
            let b = &wrapped[(k + 1) % n];
            for &t in &ts {
                let t_b = 1.0 - t;
                let xtan_a = a.partial_u(1.0, t);
                let xtan_b = b.partial_v(t_b, 0.0);
                let t_spoke = a.partial_v(1.0, t);
                let cross = xtan_a.cross(&xtan_b);
                let numer = (cross.x * t_spoke.x + cross.y * t_spoke.y + cross.z * t_spoke.z).abs();
                let na = (xtan_a.x * xtan_a.x + xtan_a.y * xtan_a.y + xtan_a.z * xtan_a.z).sqrt();
                let nb = (xtan_b.x * xtan_b.x + xtan_b.y * xtan_b.y + xtan_b.z * xtan_b.z).sqrt();
                let ns =
                    (t_spoke.x * t_spoke.x + t_spoke.y * t_spoke.y + t_spoke.z * t_spoke.z).sqrt();
                let denom = na * nb * ns;
                let coplan = if denom > 1e-300 { numer / denom } else { 0.0 };
                out.push(coplan);
            }
        }
        // Option D: fold-penalty residuals. Sample |∂u × ∂v| at
        // an interior grid per patch; emit `max(0, threshold - v)`
        // scaled by FOLD_PENALTY_WEIGHT per sample.
        let grid_ts: Vec<f64> = (0..FOLD_PENALTY_SAMPLES)
            .map(|k| 0.1 + 0.8 * (k as f64) / ((FOLD_PENALTY_SAMPLES - 1) as f64))
            .collect();
        for patch in &wrapped {
            for &u in &grid_ts {
                for &v in &grid_ts {
                    let du = patch.partial_u(u, v);
                    let dv = patch.partial_v(u, v);
                    let cross = du.cross(&dv);
                    let mag = (cross.x * cross.x + cross.y * cross.y + cross.z * cross.z).sqrt();
                    let penalty = (FOLD_PENALTY_THRESHOLD - mag).max(0.0);
                    out.push(FOLD_PENALTY_WEIGHT * penalty);
                }
            }
        }
        out
    }
}

impl LmEval for PolarLmEval<'_> {
    fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>> {
        let vars = self.unpack(x)?;
        let patches = self.build(&vars)?;
        Some(self.compute_residuals(&patches))
    }

    fn fold_safe(&mut self, x: &[f64]) -> bool {
        use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
        let Some(vars) = self.unpack(x) else {
            return false;
        };
        let Some(patches) = self.build(&vars) else {
            return false;
        };
        for g in &patches {
            let gp =
                GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid: *g }));
            if gp.self_intersects_with_tolerance(SELF_INT_TOL) {
                return false;
            }
        }
        true
    }
}

/// Pack a `Bi6CapFreeParams` into the flat parameter vector layout
/// expected by [`PolarLmEval::unpack`].
fn pack(base: &Bi6CapInput, params: &Bi6CapFreeParams) -> Vec<f64> {
    let n = params.p_23.len();
    let expected = 3 * (2 * n + 4 + n * EXTRA_DEEP_CELLS.len());
    let mut out = Vec::with_capacity(expected);
    for p in [&params.p0_01, &params.p0_02, &params.p0_11, &params.p0_22] {
        out.push(p.x);
        out.push(p.y);
        out.push(p.z);
    }
    for k in 0..n {
        out.push(params.p_23[k].x);
        out.push(params.p_23[k].y);
        out.push(params.p_23[k].z);
        out.push(params.p_33[k].x);
        out.push(params.p_33[k].y);
        out.push(params.p_33[k].z);
    }
    let extra = initial_extra_deep_cells(base, params);
    for cells in &extra {
        for p in cells {
            out.push(p.x);
            out.push(p.y);
            out.push(p.z);
        }
    }
    out
}

fn initial_extra_deep_cells(
    base: &Bi6CapInput,
    params: &Bi6CapFreeParams,
) -> Vec<[Point3; EXTRA_DEEP_CELLS.len()]> {
    let mut out = vec![[base.c7; EXTRA_DEEP_CELLS.len()]; base.n];
    if let Some(input) = boundary_owned_deep_fill_input(base, params, 1.0)
        && let Some(filler) = input.interior_filler
    {
        for (k, cells) in out.iter_mut().enumerate() {
            for (idx, &(i, j)) in EXTRA_DEEP_CELLS.iter().enumerate() {
                if let Some(p) = filler[k][i][j] {
                    cells[idx] = p;
                }
            }
        }
    }
    out
}

/// Attempt L-M refinement. Tries multiple warm-starts:
///
/// 1. The caller-supplied F_κ output (best guess for where L-M
///    should converge — but may be already folded).
/// 2. Pulled-back variants at α = 0.5, 0.25, 0.05 — interior CPs
///    contracted toward the cap-local star to reduce curvature and
///    likely find a fold-safe starting point.
/// 3. All star — maximally-degenerate cap-local start; usually
///    fold-safe because interior CPs collapse to the cap centre.
///
/// Accepts the first warm-start that is itself fold-safe; then
/// lets L-M refine. If none of the warm-starts are fold-safe,
/// returns `None`.
///
/// Returns `Some(patches_bi9)` if L-M converges to a fold-safe
/// configuration; `None` if every warm-start produces a folded
/// initial config OR L-M fails to converge.
///
/// This is the entry point used by
/// [`promote_decomp_to_kp_polar_bi6`] as a post-ladder escape
/// when F_κ fails.
///
/// [`promote_decomp_to_kp_polar_bi6`]: super::kp_polar_bi6_promote::promote_decomp_to_kp_polar_bi6
pub(super) fn try_polar_lm(
    base: &Bi6CapInput,
    warm_start: &Bi6CapFreeParams,
) -> Option<Vec<[[Point3; 10]; 10]>> {
    try_polar_lm_with_report(base, warm_start).patches
}

pub(super) fn try_polar_lm_with_report(
    base: &Bi6CapInput,
    warm_start: &Bi6CapFreeParams,
) -> PolarLmReport {
    let anchor = base.c7;
    let contract = |s: f64| contract_params(warm_start, anchor, s);
    let collapsed = collapsed_params(warm_start.p_23.len(), anchor);
    let candidates = [
        contract(1.0),
        contract(0.5),
        contract(0.25),
        contract(0.05),
        collapsed,
    ];
    let n_dof = pack(base, warm_start).len();
    let estimated_evaluations = candidates.len() * MAX_ITER * (n_dof + 1);
    if estimated_evaluations > MAX_FINITE_DIFFERENCE_EVALUATIONS {
        return PolarLmReport {
            patches: None,
            rejection_reason: Some(
                "full_deep_lm_no_fold_safe_iterate;lm_fd_evaluation_budget_exceeded",
            ),
        };
    }

    let mut saw_feasible_start = false;
    for (idx, candidate) in candidates.iter().enumerate() {
        match try_single_lm(base, candidate) {
            SingleLmOutcome::Accepted(patches) => {
                if popkern_core::trace::trace_enabled("kp-polar") {
                    eprintln!("[kp_polar_bi6 lm] warm-start #{idx} produced fold-safe iterate");
                }
                return PolarLmReport {
                    patches: Some(patches),
                    rejection_reason: None,
                };
            }
            SingleLmOutcome::InfeasibleStart => {}
            SingleLmOutcome::InfeasibleStartRejected => {}
            SingleLmOutcome::RejectedAfterFeasibleStart => {
                saw_feasible_start = true;
            }
        }
    }
    PolarLmReport {
        patches: None,
        rejection_reason: Some(if saw_feasible_start {
            "feasible_warm_start_optimized_to_fold"
        } else {
            "full_deep_lm_no_fold_safe_iterate"
        }),
    }
}

/// Run L-M once with a single warm-start. Returns `None` if the
/// warm-start itself isn't fold-safe and the infeasible-start
/// solve cannot produce a fold-safe final iterate, OR L-M doesn't
/// produce a fold-safe final iterate.
fn try_single_lm(base: &Bi6CapInput, warm_start: &Bi6CapFreeParams) -> SingleLmOutcome {
    let eval_check = PolarLmEval { base };
    let x0 = pack(base, warm_start);

    let mut eval_mut = PolarLmEval { base };
    let feasible_start = eval_mut.fold_safe(&x0);

    let config = LmConfig {
        fd_eps: FD_EPS,
        lambda_init: LAMBDA_INIT,
        lambda_max: LAMBDA_MAX,
        max_iter: MAX_ITER,
        target_residual: TARGET_RESIDUAL,
        require_feasible_start: false,
    };
    let n = eval_check.n();
    let result = run_lm(PolarLmEval { base }, x0, &config, &format!("polar N={n}"));
    let Some(vars) = eval_check.unpack(&result.x) else {
        return SingleLmOutcome::RejectedAfterFeasibleStart;
    };
    let Some(patches) = eval_check.build(&vars) else {
        return SingleLmOutcome::RejectedAfterFeasibleStart;
    };
    use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
    for g in &patches {
        let gp = GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid: *g }));
        if gp.self_intersects_with_tolerance(SELF_INT_TOL) {
            return if feasible_start {
                SingleLmOutcome::RejectedAfterFeasibleStart
            } else {
                SingleLmOutcome::InfeasibleStartRejected
            };
        }
    }
    SingleLmOutcome::Accepted(patches)
}

fn contract_params(p: &Bi6CapFreeParams, anchor: Point3, s: f64) -> Bi6CapFreeParams {
    let sp = |q: Point3| {
        let v = q - anchor;
        anchor + v * s
    };
    Bi6CapFreeParams {
        p0_01: sp(p.p0_01),
        p0_02: sp(p.p0_02),
        p0_11: sp(p.p0_11),
        p0_22: sp(p.p0_22),
        p_23: p.p_23.iter().map(|q| sp(*q)).collect(),
        p_33: p.p_33.iter().map(|q| sp(*q)).collect(),
    }
}

fn collapsed_params(n: usize, anchor: Point3) -> Bi6CapFreeParams {
    Bi6CapFreeParams {
        p0_01: anchor,
        p0_02: anchor,
        p0_11: anchor,
        p0_22: anchor,
        p_23: vec![anchor; n],
        p_33: vec![anchor; n],
    }
}
