//! Final-gate micro curvature solve for near-target row-band witnesses.

use popkern_core::topology::VertexHandle;

use super::super::chain_discovery::ChainTrim;
use super::super::chain_fit::ChainRowKey;
use super::DerivedRowBand;
use super::curvature_common::{
    apply_deltas, limited_row_deltas, max_normal_aligned_g2_residual, max_trim_disagreement,
    solve_curvature_system, top_curvature_samples, write_repaired_rows,
};

const GATE_MICRO_SAMPLE_SUBDIVISIONS: usize = 8;
const COMMERCIAL_G1_TARGET: f64 = 1.25e-3;
const COMMERCIAL_G2_TARGET: f64 = 5.0e-2;
const GATE_MICRO_ACTIVATION_CEILING: f64 = 5.1e-2;
const GATE_MICRO_MAX_CONTROL_FRACTION: f64 = 0.015;
const GATE_MICRO_MAX_PASSES: usize = 4;
const GATE_MICRO_IMPROVEMENT_EPS: f64 = 1.0e-12;
const GATE_MICRO_DAMPING_FACTORS: [f64; 9] = [
    0.5,
    0.25,
    0.125,
    0.0625,
    0.03125,
    0.015625,
    0.0078125,
    0.00390625,
    0.001953125,
];

pub(super) fn repair_curvature_gate_micro_solve(
    bands: &mut [DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> String {
    if bands.len() < 2 {
        return declined_detail("insufficient_bands", 0.0, 0.0);
    }
    let mut before_g2 =
        max_normal_aligned_g2_residual(bands, chain_trims, GATE_MICRO_SAMPLE_SUBDIVISIONS);
    let mut before_g1 = max_trim_disagreement(bands, chain_trims, GATE_MICRO_SAMPLE_SUBDIVISIONS);
    let initial_g2 = before_g2;
    let initial_g1 = before_g1;
    let mut pass_count = 0usize;
    let mut total_samples = 0usize;
    let mut total_repaired_rows = 0usize;
    let mut total_limited_rows = 0usize;
    let mut max_control_delta = 0.0_f64;
    let mut worst_detail = "none".to_string();
    let mut status = "outside_activation_band";

    for _ in 0..GATE_MICRO_MAX_PASSES {
        if before_g2 <= COMMERCIAL_G2_TARGET {
            status = if pass_count == 0 {
                "commercial_g2_already_met"
            } else {
                "commercial_g2_met"
            };
            break;
        }
        if before_g1 > COMMERCIAL_G1_TARGET || before_g2 > GATE_MICRO_ACTIVATION_CEILING {
            status = "outside_activation_band";
            break;
        }
        let active_samples = top_curvature_samples(
            bands,
            chain_trims,
            COMMERCIAL_G2_TARGET,
            1,
            GATE_MICRO_SAMPLE_SUBDIVISIONS,
        );
        total_samples += active_samples.len();
        let Some(sample) = active_samples.first() else {
            status = "no_active_gate_sample";
            break;
        };
        worst_detail = sample.detail.clone();
        let Some(solution) = solve_curvature_system(bands, &active_samples) else {
            status = "gate_micro_system_singular";
            break;
        };
        let (deltas, limited_rows, max_requested_delta) =
            limited_row_deltas(bands, &solution, GATE_MICRO_MAX_CONTROL_FRACTION);
        let mut best: Option<(f64, f64, f64, Vec<DerivedRowBand>)> = None;
        for damping in GATE_MICRO_DAMPING_FACTORS {
            let mut candidate = bands.to_vec();
            apply_deltas(&mut candidate, &deltas, damping);
            let after_g2 = max_normal_aligned_g2_residual(
                &candidate,
                chain_trims,
                GATE_MICRO_SAMPLE_SUBDIVISIONS,
            );
            if after_g2 + GATE_MICRO_IMPROVEMENT_EPS >= before_g2 {
                continue;
            }
            let after_g1 =
                max_trim_disagreement(&candidate, chain_trims, GATE_MICRO_SAMPLE_SUBDIVISIONS);
            if after_g1 > COMMERCIAL_G1_TARGET + 1.0e-10 {
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
            status = "no_quality_safe_gate_micro_improvement";
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
        status = "kept";
    }

    format!(
        "curvature_gate_micro_solve_status={};curvature_gate_micro_sample_count={};curvature_gate_micro_pass_count={};curvature_gate_micro_repaired_row_count={};curvature_gate_micro_limited_row_count={};max_curvature_gate_micro_control_delta={:.6e};max_curvature_gate_micro_residual_before={:.15e};max_curvature_gate_micro_residual_after={:.15e};max_curvature_gate_micro_g1_before={:.15e};max_curvature_gate_micro_g1_after={:.15e};worst_curvature_gate_micro_detail=[{}]",
        status,
        total_samples,
        pass_count,
        total_repaired_rows,
        total_limited_rows,
        max_control_delta,
        initial_g2,
        before_g2,
        initial_g1,
        before_g1,
        worst_detail
    )
}

fn declined_detail(reason: &'static str, initial_g2: f64, initial_g1: f64) -> String {
    format!(
        "curvature_gate_micro_solve_status={reason};curvature_gate_micro_sample_count=0;curvature_gate_micro_pass_count=0;curvature_gate_micro_repaired_row_count=0;curvature_gate_micro_limited_row_count=0;max_curvature_gate_micro_control_delta=0.000000e0;max_curvature_gate_micro_residual_before={initial_g2:.15e};max_curvature_gate_micro_residual_after={initial_g2:.15e};max_curvature_gate_micro_g1_before={initial_g1:.15e};max_curvature_gate_micro_g1_after={initial_g1:.15e};worst_curvature_gate_micro_detail=[none]"
    )
}
