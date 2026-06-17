//! Active-set normal-curvature solve for constructive row bands.

use popkern_core::topology::VertexHandle;

use super::super::chain_discovery::ChainTrim;
use super::super::chain_fit::ChainRowKey;
use super::DerivedRowBand;
use super::curvature_common::{
    apply_deltas, limited_row_deltas, max_normal_aligned_g2_residual, max_trim_disagreement,
    solve_curvature_system, top_curvature_samples, write_repaired_rows,
};

const CURVATURE_SAMPLE_LIMIT: usize = 24;
const CURVATURE_SAMPLE_SUBDIVISIONS: usize = 8;
const MIN_NORMAL_ALIGNED_G2_RESIDUAL: f64 = 5.0e-2;
const MARGIN_NORMAL_ALIGNED_G2_RESIDUAL: f64 = 4.0e-2;
const MARGIN_CURVATURE_SAMPLE_LIMIT: usize = 2;
const COMMERCIAL_G1_TARGET: f64 = 1.25e-3;
const MAX_CURVATURE_CONTROL_FRACTION: f64 = 0.03;
const CURVATURE_MAX_PASSES: usize = 8;
const EXPANDED_ROW_BAND_MAX_CURVATURE_CONTROL_FRACTION: f64 = 0.06;
const EXPANDED_ROW_BAND_CURVATURE_MAX_PASSES: usize = 12;
const CURVATURE_DAMPING_FACTORS: [f64; 6] = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125];
const CURVATURE_IMPROVEMENT_EPS: f64 = 1.0e-12;

#[derive(Debug, Clone, Copy)]
struct CurvatureSolveBudget {
    max_control_fraction: f64,
    max_passes: usize,
}

pub(super) fn repair_curvature_solve(
    bands: &mut [DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> String {
    repair_curvature_solve_with_budget(
        bands,
        chain_trims,
        CurvatureSolveBudget {
            max_control_fraction: MAX_CURVATURE_CONTROL_FRACTION,
            max_passes: CURVATURE_MAX_PASSES,
        },
    )
}

pub(super) fn repair_curvature_solve_with_expanded_row_band_budget(
    bands: &mut [DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> String {
    repair_curvature_solve_with_budget(
        bands,
        chain_trims,
        CurvatureSolveBudget {
            max_control_fraction: EXPANDED_ROW_BAND_MAX_CURVATURE_CONTROL_FRACTION,
            max_passes: EXPANDED_ROW_BAND_CURVATURE_MAX_PASSES,
        },
    )
}

fn repair_curvature_solve_with_budget(
    bands: &mut [DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
    budget: CurvatureSolveBudget,
) -> String {
    if bands.len() < 2 {
        return declined_detail("insufficient_bands", 0.0, 0.0);
    }
    let mut before_g2 =
        max_normal_aligned_g2_residual(bands, chain_trims, CURVATURE_SAMPLE_SUBDIVISIONS);
    let mut before_g1 = max_trim_disagreement(bands, chain_trims, CURVATURE_SAMPLE_SUBDIVISIONS);
    let initial_g2 = before_g2;
    let initial_g1 = before_g1;
    let mut pass_count = 0usize;
    let mut total_samples = 0usize;
    let mut total_constraints = 0usize;
    let mut total_repaired_rows = 0usize;
    let mut total_limited_rows = 0usize;
    let mut max_control_delta = 0.0_f64;
    let mut best_detail = top_curvature_samples(
        bands,
        chain_trims,
        MIN_NORMAL_ALIGNED_G2_RESIDUAL,
        CURVATURE_SAMPLE_LIMIT,
        CURVATURE_SAMPLE_SUBDIVISIONS,
    )
    .first()
    .map(|sample| sample.detail.clone())
    .unwrap_or_else(|| "none".to_string());
    let mut decline_reason = "no_active_samples";

    for _ in 0..budget.max_passes {
        let (active_residual, sample_limit) = if before_g2 <= MIN_NORMAL_ALIGNED_G2_RESIDUAL {
            (
                MARGIN_NORMAL_ALIGNED_G2_RESIDUAL,
                MARGIN_CURVATURE_SAMPLE_LIMIT,
            )
        } else {
            (MIN_NORMAL_ALIGNED_G2_RESIDUAL, CURVATURE_SAMPLE_LIMIT)
        };
        let active_samples = top_curvature_samples(
            bands,
            chain_trims,
            active_residual,
            sample_limit,
            CURVATURE_SAMPLE_SUBDIVISIONS,
        );
        total_samples += active_samples.len();
        total_constraints += active_samples.len();
        if active_samples.is_empty() {
            decline_reason = if pass_count == 0 {
                "no_active_samples"
            } else {
                "converged_no_active_samples"
            };
            break;
        }
        if let Some(sample) = active_samples.first() {
            best_detail = sample.detail.clone();
        }
        let Some(solution) = solve_curvature_system(bands, &active_samples) else {
            decline_reason = "curvature_system_singular";
            break;
        };
        let (deltas, limited_rows, max_requested_delta) =
            limited_row_deltas(bands, &solution, budget.max_control_fraction);
        let mut best: Option<(f64, f64, f64, Vec<DerivedRowBand>)> = None;
        for damping in CURVATURE_DAMPING_FACTORS {
            let mut candidate = bands.to_vec();
            apply_deltas(&mut candidate, &deltas, damping);
            let after_g2 = max_normal_aligned_g2_residual(
                &candidate,
                chain_trims,
                CURVATURE_SAMPLE_SUBDIVISIONS,
            );
            if after_g2 + CURVATURE_IMPROVEMENT_EPS >= before_g2 {
                continue;
            }
            let after_g1 =
                max_trim_disagreement(&candidate, chain_trims, CURVATURE_SAMPLE_SUBDIVISIONS);
            if after_g1 > before_g1.max(COMMERCIAL_G1_TARGET) + 1.0e-10 {
                continue;
            }
            if best
                .as_ref()
                .is_none_or(|(_, best_g2, _, _)| after_g2 < *best_g2)
            {
                best = Some((damping, after_g2, after_g1, candidate));
            }
        }
        let Some((damping, after_g2, after_g1, candidate)) = best else {
            decline_reason = "no_quality_safe_improvement";
            total_limited_rows += limited_rows;
            max_control_delta = max_control_delta.max(max_requested_delta);
            break;
        };
        total_repaired_rows += write_repaired_rows(bands, candidate);
        total_limited_rows += limited_rows;
        max_control_delta = max_control_delta.max(max_requested_delta * damping);
        before_g2 = after_g2;
        before_g1 = after_g1;
        pass_count += 1;
        decline_reason = "kept";
    }

    format!(
        "curvature_solve_status={};curvature_sample_count={};curvature_pass_count={};curvature_constraint_count={};curvature_repaired_row_count={};curvature_limited_row_count={};max_curvature_control_delta={:.6e};max_curvature_residual_before={:.15e};max_curvature_residual_after={:.15e};max_curvature_g1_before={:.15e};max_curvature_g1_after={:.15e};worst_curvature_detail=[{}]",
        decline_reason,
        total_samples,
        pass_count,
        total_constraints,
        total_repaired_rows,
        total_limited_rows,
        max_control_delta,
        initial_g2,
        before_g2,
        initial_g1,
        before_g1,
        best_detail
    )
}

fn declined_detail(reason: &'static str, initial_g2: f64, initial_g1: f64) -> String {
    format!(
        "curvature_solve_status={reason};curvature_sample_count=0;curvature_pass_count=0;curvature_constraint_count=0;curvature_repaired_row_count=0;curvature_limited_row_count=0;max_curvature_control_delta=0.000000e0;max_curvature_residual_before={initial_g2:.15e};max_curvature_residual_after={initial_g2:.15e};max_curvature_g1_before={initial_g1:.15e};max_curvature_g1_after={initial_g1:.15e};worst_curvature_detail=[none]"
    )
}
