//! Local normal solve for the measured worst support-cell G1 trim.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    RectangularGregoryBiNonic, bernstein_binonic_all, bernstein_binonic_deriv_all,
};
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;

use super::chain_solver::solve_dense;
use super::worst_g1_normal_target::{
    ControlAuthority, ControlKey, ControlVariable, MatchedWorstTrim, control_variables,
    local_control_limit, local_sample_fractions, match_worst_trim, uv_on_side,
};
use super::{FullBoundaryTrim, rectangular_binonic_on_face};

const NORMAL_SOLVE_REGULARIZATION: f64 = 1.0e-8;
const LOCAL_NORMAL_DAMPING_FACTORS: [f64; 4] = [1.0, 0.5, 0.25, 0.125];
const QUALITY_COMPARISON_EPS: f64 = 1.0e-12;

#[derive(Debug, Clone)]
pub(super) struct WorstG1NormalSolveReport {
    pub(super) status: &'static str,
    pub(super) decline_reason: &'static str,
    pub(super) control_authority: &'static str,
    pub(super) sample_count: usize,
    pub(super) variable_count: usize,
    pub(super) trial_count: usize,
    pub(super) validation_reject_count: usize,
    pub(super) quality_reject_count: usize,
    pub(super) limited_control_count: usize,
    pub(super) max_control_delta: f64,
    pub(super) damping: Option<f64>,
    pub(super) best_trial_damping: Option<f64>,
    pub(super) best_trial_worst_g1: Option<f64>,
    pub(super) best_trial_local_disagreement: Option<f64>,
    pub(super) before_worst_g1: f64,
    pub(super) after_worst_g1: f64,
    pub(super) local_disagreement_before: f64,
    pub(super) local_disagreement_after: f64,
}

impl WorstG1NormalSolveReport {
    fn declined(status: &'static str, authority: ControlAuthority, before_worst_g1: f64) -> Self {
        Self {
            status,
            decline_reason: status,
            control_authority: authority.label(),
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
}

#[derive(Debug, Clone, Copy)]
struct NormalSample {
    normal: Vec3,
    raw_normal_magnitude: f64,
    tangent_u: Vec3,
    tangent_v: Vec3,
    u_basis: [f64; 10],
    v_basis: [f64; 10],
    u_deriv_basis: [f64; 10],
    v_deriv_basis: [f64; 10],
}

#[derive(Debug, Default, Clone)]
struct NormalTrialSummary {
    trial_count: usize,
    validation_reject_count: usize,
    quality_reject_count: usize,
    best_trial_damping: Option<f64>,
    best_trial_worst_g1: Option<f64>,
    best_trial_local_disagreement: Option<f64>,
}

impl NormalTrialSummary {
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

pub(super) fn reconcile_endpoint_adjacent_worst_g1_normal_sample_if_improves(
    brep: &mut BRep,
    faces: &[FaceHandle],
    trims: &[FullBoundaryTrim],
) -> Result<WorstG1NormalSolveReport, CoreError> {
    reconcile_worst_g1_normal_sample_with_authority_if_improves(
        brep,
        faces,
        trims,
        ControlAuthority::EndpointAdjacentRows,
    )
}

fn reconcile_worst_g1_normal_sample_with_authority_if_improves(
    brep: &mut BRep,
    faces: &[FaceHandle],
    trims: &[FullBoundaryTrim],
    authority: ControlAuthority,
) -> Result<WorstG1NormalSolveReport, CoreError> {
    let before_quality =
        shared_trim_quality::measure_support_cell_shared_trim_quality(brep, faces, 9);
    let Some(worst) = before_quality.worst_g1_trim else {
        return Ok(WorstG1NormalSolveReport::declined(
            "missing_worst_g1_trim",
            authority,
            before_quality.worst_normal_disagreement,
        ));
    };
    let Some(matched) = trims.iter().find_map(|trim| match_worst_trim(*trim, worst)) else {
        return Ok(WorstG1NormalSolveReport::declined(
            "matching_trim_not_found",
            authority,
            before_quality.worst_normal_disagreement,
        ));
    };
    let variables = control_variables(matched, worst.sample_fraction, authority);
    if variables.is_empty() {
        return Ok(WorstG1NormalSolveReport::declined(
            "no_owned_controls",
            authority,
            before_quality.worst_normal_disagreement,
        ));
    }
    let patch_a = rectangular_binonic_on_face(brep, matched.face_a)?;
    let patch_b = rectangular_binonic_on_face(brep, matched.face_b)?;
    let samples = local_sample_fractions(worst.sample_fraction);
    let local_before = max_local_disagreement(&patch_a, &patch_b, matched, &samples);
    let Some(solution) = solve_normal_system(&patch_a, &patch_b, matched, &variables, &samples)
    else {
        let mut report = WorstG1NormalSolveReport::declined(
            "normal_system_singular",
            authority,
            before_quality.worst_normal_disagreement,
        );
        report.sample_count = samples.len();
        report.variable_count = variables.len();
        report.local_disagreement_before = local_before;
        report.local_disagreement_after = local_before;
        return Ok(report);
    };
    let (control_deltas, limited_control_count, max_control_delta) = limited_control_deltas(
        &variables,
        &solution,
        local_control_limit(&patch_a, &patch_b, matched),
    );
    let mut best: Option<(
        f64,
        shared_trim_quality::Kp20SupportCellSharedTrimQuality,
        f64,
        BRep,
    )> = None;
    let mut trial_summary = NormalTrialSummary::default();
    for damping in LOCAL_NORMAL_DAMPING_FACTORS {
        trial_summary.trial_count += 1;
        let mut candidate = brep.clone();
        apply_control_deltas(&mut candidate, &control_deltas, damping)?;
        if candidate.validate().is_err() {
            trial_summary.validation_reject_count += 1;
            continue;
        }
        let after_quality =
            shared_trim_quality::measure_support_cell_shared_trim_quality(&candidate, faces, 9);
        let candidate_patch_a = rectangular_binonic_on_face(&candidate, matched.face_a)?;
        let candidate_patch_b = rectangular_binonic_on_face(&candidate, matched.face_b)?;
        let local_after =
            max_local_disagreement(&candidate_patch_a, &candidate_patch_b, matched, &samples);
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
            control_authority: authority.label(),
            sample_count: samples.len(),
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
        control_authority: authority.label(),
        sample_count: samples.len(),
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

fn solve_normal_system(
    patch_a: &RectangularGregoryBiNonic,
    patch_b: &RectangularGregoryBiNonic,
    matched: MatchedWorstTrim,
    variables: &[ControlVariable],
    samples: &[(f64, f64)],
) -> Option<Vec<f64>> {
    let mut matrix = vec![vec![0.0_f64; variables.len()]; variables.len()];
    let mut rhs = vec![0.0_f64; variables.len()];
    for &(fraction, weight) in samples {
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
                } else {
                    normal_sensitivity(sample_b, *variable) * -orientation
                };
                row[column] = component(sensitivity, axis) * weight_sqrt;
            }
            add_normal_equation(
                &mut matrix,
                &mut rhs,
                &row,
                -component(diff, axis) * weight_sqrt,
            );
        }
    }
    for (index, row) in matrix.iter_mut().enumerate() {
        row[index] += NORMAL_SOLVE_REGULARIZATION;
    }
    solve_dense(matrix, rhs)
}

fn normal_sample(patch: &RectangularGregoryBiNonic, uv: (f64, f64)) -> Option<NormalSample> {
    let u_basis = bernstein_binonic_all(uv.0);
    let v_basis = bernstein_binonic_all(uv.1);
    let u_deriv_basis = bernstein_binonic_deriv_all(uv.0);
    let v_deriv_basis = bernstein_binonic_deriv_all(uv.1);
    let tangent_u = patch.partial_u(uv.0, uv.1);
    let tangent_v = patch.partial_v(uv.0, uv.1);
    let raw_normal = tangent_u.cross(&tangent_v);
    let raw_normal_magnitude = raw_normal.dot(&raw_normal).sqrt();
    Some(NormalSample {
        normal: raw_normal.normalized()?,
        raw_normal_magnitude,
        tangent_u,
        tangent_v,
        u_basis,
        v_basis,
        u_deriv_basis,
        v_deriv_basis,
    })
}

fn normal_sensitivity(sample: NormalSample, variable: ControlVariable) -> Vec3 {
    let direction = axis_direction(variable.axis);
    let delta_u =
        direction * (sample.u_deriv_basis[variable.key.u] * sample.v_basis[variable.key.v]);
    let delta_v =
        direction * (sample.u_basis[variable.key.u] * sample.v_deriv_basis[variable.key.v]);
    let raw_delta = delta_u.cross(&sample.tangent_v) + sample.tangent_u.cross(&delta_v);
    (raw_delta - sample.normal * raw_delta.dot(&sample.normal))
        * (1.0 / sample.raw_normal_magnitude)
}

fn limited_control_deltas(
    variables: &[ControlVariable],
    solution: &[f64],
    limit: f64,
) -> (HashMap<ControlKey, Vec3>, usize, f64) {
    let mut deltas: HashMap<ControlKey, Vec3> = HashMap::new();
    for (variable, value) in variables.iter().zip(solution) {
        deltas
            .entry(variable.key)
            .and_modify(|delta| *delta = *delta + axis_direction(variable.axis) * *value)
            .or_insert(axis_direction(variable.axis) * *value);
    }
    let mut limited_count = 0usize;
    let mut max_delta = 0.0_f64;
    for delta in deltas.values_mut() {
        let magnitude = delta.dot(delta).sqrt();
        if magnitude > limit {
            *delta = *delta * (limit / magnitude);
            limited_count += 1;
            max_delta = max_delta.max(limit);
        } else {
            max_delta = max_delta.max(magnitude);
        }
    }
    (deltas, limited_count, max_delta)
}

fn apply_control_deltas(
    brep: &mut BRep,
    deltas: &HashMap<ControlKey, Vec3>,
    damping: f64,
) -> Result<(), CoreError> {
    let mut patches: HashMap<FaceHandle, RectangularGregoryBiNonic> = HashMap::new();
    for (&key, delta) in deltas {
        if let std::collections::hash_map::Entry::Vacant(entry) = patches.entry(key.face) {
            entry.insert(rectangular_binonic_on_face(brep, key.face)?);
        }
        let Some(patch) = patches.get_mut(&key.face) else {
            return Err(CoreError::PreconditionFailed(
                "kp20 worst-G1 normal solve control patch missing".into(),
            ));
        };
        patch.grid[key.u][key.v] = patch.grid[key.u][key.v] + *delta * damping;
    }
    for (face, patch) in patches {
        super::set_rectangular_binonic_on_face(brep, face, patch)?;
    }
    Ok(())
}

fn max_local_disagreement(
    patch_a: &RectangularGregoryBiNonic,
    patch_b: &RectangularGregoryBiNonic,
    matched: MatchedWorstTrim,
    samples: &[(f64, f64)],
) -> f64 {
    let mut worst = 0.0_f64;
    for &(fraction, _) in samples {
        let Some(sample_a) = normal_sample(
            patch_a,
            uv_on_side(matched.side_a, matched.start_a, matched.end_a, fraction),
        ) else {
            continue;
        };
        let Some(sample_b) = normal_sample(
            patch_b,
            uv_on_side(
                matched.side_b,
                matched.start_b,
                matched.end_b,
                1.0 - fraction,
            ),
        ) else {
            continue;
        };
        worst = worst.max(1.0 - sample_a.normal.dot(&sample_b.normal).abs());
    }
    worst
}

fn add_normal_equation(
    matrix: &mut [Vec<f64>],
    rhs: &mut [f64],
    coefficients: &[f64],
    target: f64,
) {
    for i in 0..coefficients.len() {
        if coefficients[i].abs() <= 1.0e-18 {
            continue;
        }
        rhs[i] += coefficients[i] * target;
        for j in 0..coefficients.len() {
            if coefficients[j].abs() > 1.0e-18 {
                matrix[i][j] += coefficients[i] * coefficients[j];
            }
        }
    }
}

fn axis_direction(axis: usize) -> Vec3 {
    match axis {
        0 => Vec3::new(1.0, 0.0, 0.0),
        1 => Vec3::new(0.0, 1.0, 0.0),
        _ => Vec3::new(0.0, 0.0, 1.0),
    }
}

fn component(vector: Vec3, axis: usize) -> f64 {
    match axis {
        0 => vector.x,
        1 => vector.y,
        _ => vector.z,
    }
}
