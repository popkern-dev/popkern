//! Active final-gate normal solve for the top shared-trim G1 samples.

use std::collections::{HashMap, HashSet};

use popkern_core::error::CoreError;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;

use super::chain_solver::solve_dense;
use super::top_g1_normal_math::{
    add_normal_equation, apply_control_deltas, component, limited_control_deltas, normal_sample,
    normal_sensitivity,
};
use super::worst_g1_normal_solve::WorstG1NormalSolveReport;
use super::worst_g1_normal_target::{
    ControlAuthority, ControlVariable, MatchedWorstTrim, control_variables, local_control_limit,
    local_sample_fractions, match_worst_trim, uv_on_side,
};
use super::{FullBoundaryTrim, rectangular_binonic_on_face};

const TOP_G1_NORMAL_REGULARIZATION: f64 = 1.0e-8;
const TOP_G1_NORMAL_DAMPING_FACTORS: [f64; 4] = [1.0, 0.5, 0.25, 0.125];
const QUALITY_COMPARISON_EPS: f64 = 1.0e-12;

#[derive(Debug, Clone, Copy)]
struct ActiveGateSample {
    matched: MatchedWorstTrim,
    fraction: f64,
}

#[derive(Debug, Default, Clone)]
struct TrialSummary {
    trial_count: usize,
    validation_reject_count: usize,
    quality_reject_count: usize,
    best_trial_damping: Option<f64>,
    best_trial_worst_g1: Option<f64>,
    best_trial_local_disagreement: Option<f64>,
}

impl TrialSummary {
    fn record_quality_trial(&mut self, damping: f64, worst_g1: f64, local_disagreement: f64) {
        if self.best_trial_worst_g1.is_none_or(|best| worst_g1 < best) {
            self.best_trial_damping = Some(damping);
            self.best_trial_worst_g1 = Some(worst_g1);
            self.best_trial_local_disagreement = Some(local_disagreement);
        }
    }

    fn decline_reason(&self, before_worst_g1: f64) -> &'static str {
        if self.trial_count == 0 {
            return "no_trials";
        }
        if self.validation_reject_count == self.trial_count {
            return "all_trials_invalid";
        }
        if self
            .best_trial_worst_g1
            .is_some_and(|best| best + QUALITY_COMPARISON_EPS >= before_worst_g1)
        {
            return "no_final_g1_improvement";
        }
        "quality_gate_declined"
    }
}

pub(super) fn reconcile_top_g1_normal_samples_if_improves(
    brep: &mut BRep,
    faces: &[FaceHandle],
    trims: &[FullBoundaryTrim],
) -> Result<WorstG1NormalSolveReport, CoreError> {
    let authority = ControlAuthority::EndpointAdjacentRows;
    let before_quality =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    let samples = active_gate_samples(&before_quality, trims);
    if samples.is_empty() {
        return Ok(declined(
            "missing_top_g1_samples",
            before_quality.worst_normal_disagreement,
        ));
    }
    let variables = active_control_variables(&samples, authority);
    if variables.is_empty() {
        return Ok(declined(
            "no_owned_controls",
            before_quality.worst_normal_disagreement,
        ));
    }
    let patches = load_sample_patches(brep, &samples)?;
    let local_before = max_active_disagreement(&patches, &samples);
    let Some(solution) = solve_normal_system(&patches, &samples, &variables) else {
        let mut report = declined(
            "normal_system_singular",
            before_quality.worst_normal_disagreement,
        );
        report.sample_count = expanded_sample_count(&samples);
        report.variable_count = variables.len();
        report.local_disagreement_before = local_before;
        report.local_disagreement_after = local_before;
        return Ok(report);
    };
    let control_limit = samples
        .iter()
        .map(|sample| {
            local_control_limit(
                &patches[&sample.matched.face_a],
                &patches[&sample.matched.face_b],
                sample.matched,
            )
        })
        .fold(0.0_f64, f64::max);
    let (control_deltas, limited_control_count, max_control_delta) =
        limited_control_deltas(&variables, &solution, control_limit);
    let mut best: Option<(
        f64,
        shared_trim_quality::Kp20SupportCellSharedTrimQuality,
        f64,
        BRep,
    )> = None;
    let mut trial_summary = TrialSummary::default();
    for damping in TOP_G1_NORMAL_DAMPING_FACTORS {
        trial_summary.trial_count += 1;
        let mut candidate = brep.clone();
        apply_control_deltas(&mut candidate, &control_deltas, damping)?;
        if candidate.validate().is_err() {
            trial_summary.validation_reject_count += 1;
            continue;
        }
        let after_quality =
            shared_trim_quality::measure_support_cell_shared_trim_quality(&candidate, faces, 9);
        let candidate_patches = load_sample_patches(&candidate, &samples)?;
        let local_after = max_active_disagreement(&candidate_patches, &samples);
        trial_summary.record_quality_trial(
            damping,
            after_quality.worst_normal_disagreement,
            local_after,
        );
        if after_quality.worst_g0_gap > before_quality.worst_g0_gap + QUALITY_COMPARISON_EPS
            || after_quality.pcurve_out_of_unit_count > before_quality.pcurve_out_of_unit_count
            || after_quality.degenerate_pcurve_count > before_quality.degenerate_pcurve_count
            || after_quality.worst_normal_disagreement + QUALITY_COMPARISON_EPS
                >= before_quality.worst_normal_disagreement
        {
            trial_summary.quality_reject_count += 1;
            continue;
        }
        if best.as_ref().is_none_or(|(_, quality, _, _)| {
            after_quality.worst_normal_disagreement < quality.worst_normal_disagreement
        }) {
            best = Some((damping, after_quality, local_after, candidate));
        }
    }
    let Some((damping, after_quality, local_after, candidate)) = best else {
        return Ok(WorstG1NormalSolveReport {
            status: "declined",
            decline_reason: trial_summary.decline_reason(before_quality.worst_normal_disagreement),
            control_authority: "top_g1_endpoint_adjacent_rows",
            sample_count: expanded_sample_count(&samples),
            variable_count: variables.len(),
            trial_count: trial_summary.trial_count,
            validation_reject_count: trial_summary.validation_reject_count,
            quality_reject_count: trial_summary.quality_reject_count,
            limited_control_count,
            max_control_delta,
            damping: None,
            best_trial_damping: trial_summary.best_trial_damping,
            best_trial_worst_g1: trial_summary.best_trial_worst_g1,
            best_trial_local_disagreement: trial_summary.best_trial_local_disagreement,
            before_worst_g1: before_quality.worst_normal_disagreement,
            after_worst_g1: before_quality.worst_normal_disagreement,
            local_disagreement_before: local_before,
            local_disagreement_after: local_before,
        });
    };
    *brep = candidate;
    Ok(WorstG1NormalSolveReport {
        status: "kept",
        decline_reason: "kept",
        control_authority: "top_g1_endpoint_adjacent_rows",
        sample_count: expanded_sample_count(&samples),
        variable_count: variables.len(),
        trial_count: trial_summary.trial_count,
        validation_reject_count: trial_summary.validation_reject_count,
        quality_reject_count: trial_summary.quality_reject_count,
        limited_control_count,
        max_control_delta: max_control_delta * damping,
        damping: Some(damping),
        best_trial_damping: Some(damping),
        best_trial_worst_g1: Some(after_quality.worst_normal_disagreement),
        best_trial_local_disagreement: Some(local_after),
        before_worst_g1: before_quality.worst_normal_disagreement,
        after_worst_g1: after_quality.worst_normal_disagreement,
        local_disagreement_before: local_before,
        local_disagreement_after: local_after,
    })
}

fn declined(status: &'static str, before_worst_g1: f64) -> WorstG1NormalSolveReport {
    WorstG1NormalSolveReport {
        status,
        decline_reason: status,
        control_authority: "top_g1_endpoint_adjacent_rows",
        sample_count: 0,
        variable_count: 0,
        trial_count: 0,
        validation_reject_count: 0,
        quality_reject_count: 0,
        limited_control_count: 0,
        max_control_delta: 0.0,
        damping: None,
        best_trial_damping: None,
        best_trial_worst_g1: None,
        best_trial_local_disagreement: None,
        before_worst_g1,
        after_worst_g1: before_worst_g1,
        local_disagreement_before: 0.0,
        local_disagreement_after: 0.0,
    }
}

fn active_gate_samples(
    quality: &shared_trim_quality::Kp20SupportCellSharedTrimQuality,
    trims: &[FullBoundaryTrim],
) -> Vec<ActiveGateSample> {
    let mut samples = Vec::new();
    for sample in quality.top_g1_samples.iter().flatten() {
        if let Some(matched) = trims
            .iter()
            .find_map(|trim| match_worst_trim(*trim, sample.trim))
        {
            samples.push(ActiveGateSample {
                matched,
                fraction: sample.trim.sample_fraction,
            });
        }
    }
    samples
}

fn active_control_variables(
    samples: &[ActiveGateSample],
    authority: ControlAuthority,
) -> Vec<ControlVariable> {
    let mut seen = HashSet::new();
    let mut variables = Vec::new();
    for sample in samples {
        for variable in control_variables(sample.matched, sample.fraction, authority) {
            if seen.insert((variable.key, variable.axis)) {
                variables.push(variable);
            }
        }
    }
    variables
}

fn load_sample_patches(
    brep: &BRep,
    samples: &[ActiveGateSample],
) -> Result<HashMap<FaceHandle, RectangularGregoryBiNonic>, CoreError> {
    let mut patches = HashMap::new();
    for sample in samples {
        for face in [sample.matched.face_a, sample.matched.face_b] {
            if let std::collections::hash_map::Entry::Vacant(entry) = patches.entry(face) {
                entry.insert(rectangular_binonic_on_face(brep, face)?);
            }
        }
    }
    Ok(patches)
}

fn expanded_sample_count(samples: &[ActiveGateSample]) -> usize {
    samples
        .iter()
        .map(|sample| local_sample_fractions(sample.fraction).len())
        .sum()
}

fn solve_normal_system(
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
    samples: &[ActiveGateSample],
    variables: &[ControlVariable],
) -> Option<Vec<f64>> {
    let mut matrix = vec![vec![0.0_f64; variables.len()]; variables.len()];
    let mut rhs = vec![0.0_f64; variables.len()];
    for active in samples {
        let patch_a = &patches[&active.matched.face_a];
        let patch_b = &patches[&active.matched.face_b];
        for (fraction, weight) in local_sample_fractions(active.fraction) {
            add_sample_equations(
                &mut matrix,
                &mut rhs,
                variables,
                active.matched,
                patch_a,
                patch_b,
                fraction,
                weight,
            )?;
        }
    }
    for (index, row) in matrix.iter_mut().enumerate() {
        row[index] += TOP_G1_NORMAL_REGULARIZATION;
    }
    solve_dense(matrix, rhs)
}

#[expect(clippy::too_many_arguments, reason = "one normal-equation sample")]
fn add_sample_equations(
    matrix: &mut [Vec<f64>],
    rhs: &mut [f64],
    variables: &[ControlVariable],
    matched: MatchedWorstTrim,
    patch_a: &RectangularGregoryBiNonic,
    patch_b: &RectangularGregoryBiNonic,
    fraction: f64,
    weight: f64,
) -> Option<()> {
    let sample_a = normal_sample(
        patch_a,
        uv_on_side(matched.side_a, matched.start_a, matched.end_a, fraction),
    )?;
    let sample_b = normal_sample(
        patch_b,
        uv_on_side(
            matched.side_b,
            matched.start_b,
            matched.end_b,
            1.0 - fraction,
        ),
    )?;
    let orientation = if sample_a.normal.dot(&sample_b.normal) < 0.0 {
        -1.0
    } else {
        1.0
    };
    let diff = sample_a.normal - sample_b.normal * orientation;
    let weight_sqrt = weight.sqrt();
    for axis in 0..3 {
        let mut row = vec![0.0_f64; variables.len()];
        for (column, variable) in variables.iter().enumerate() {
            let sensitivity = if variable.key.face == matched.face_a {
                normal_sensitivity(sample_a, *variable)
            } else if variable.key.face == matched.face_b {
                normal_sensitivity(sample_b, *variable) * -orientation
            } else {
                Vec3::new(0.0, 0.0, 0.0)
            };
            row[column] = component(sensitivity, axis) * weight_sqrt;
        }
        add_normal_equation(matrix, rhs, &row, -component(diff, axis) * weight_sqrt);
    }
    Some(())
}

fn max_active_disagreement(
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
    samples: &[ActiveGateSample],
) -> f64 {
    let mut worst = 0.0_f64;
    for active in samples {
        let patch_a = &patches[&active.matched.face_a];
        let patch_b = &patches[&active.matched.face_b];
        for (fraction, _) in local_sample_fractions(active.fraction) {
            let Some(sample_a) = normal_sample(
                patch_a,
                uv_on_side(
                    active.matched.side_a,
                    active.matched.start_a,
                    active.matched.end_a,
                    fraction,
                ),
            ) else {
                continue;
            };
            let Some(sample_b) = normal_sample(
                patch_b,
                uv_on_side(
                    active.matched.side_b,
                    active.matched.start_b,
                    active.matched.end_b,
                    1.0 - fraction,
                ),
            ) else {
                continue;
            };
            worst = worst.max(1.0 - sample_a.normal.dot(&sample_b.normal).abs());
        }
    }
    worst
}
