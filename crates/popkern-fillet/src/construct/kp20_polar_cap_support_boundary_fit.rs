//! Displacement-constrained row-fit selection for external KP20 cap/support seams.

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;

use super::{BoundarySide, row_fit};

pub(super) const MAX_BOUNDARY_ROW_DISPLACEMENT: f64 = 0.1;
const ROW_FIT_REGULARIZATION_CANDIDATES: [f64; 10] = [
    1.0e2, 1.0e1, 1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-6, 1.0e-8, 0.0,
];

pub(super) struct BoundaryFitCandidate {
    pub(super) patch: RectangularGregoryBiNonic,
    pub(super) sample_gap: f64,
    pub(super) sample_parameter: f64,
    pub(super) displacement: f64,
    pub(super) regularization: f64,
    pub(super) candidate_count: usize,
    pub(super) constrained_candidate_count: usize,
    pub(super) constrained_singular_candidates: usize,
    pub(super) hard_constraint_count: usize,
    pub(super) hard_constraint_residual: f64,
    pub(super) displacement_scaled_candidates: usize,
    pub(super) displacement_rejected_candidates: usize,
    pub(super) within_displacement: bool,
}

pub(super) fn best_boundary_fit_candidate(
    before: &RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[row_fit::BoundaryRowFitSample],
    hard_constraints: &[row_fit::BoundaryRowFitConstraint],
) -> BoundaryFitCandidate {
    let mut best_bounded = None;
    let mut candidate_count = 0;
    let mut constrained_candidate_count = 0;
    let mut constrained_singular_candidates = 0;
    let mut scaled = 0;
    let mut rejected = 0;
    for regularization in ROW_FIT_REGULARIZATION_CANDIDATES {
        let mut patch = before.clone();
        row_fit::fit_boundary_samples_with_regularization(
            &mut patch,
            side,
            samples,
            regularization,
        );
        row_fit::snap_boundary_endpoints_to_samples(&mut patch, side, samples, 1.0e-9);
        candidate_count += 1;
        let candidate = measured_candidate(
            before,
            side,
            samples,
            hard_constraints,
            patch,
            regularization,
        );
        record_bounded_candidate(
            before,
            side,
            samples,
            hard_constraints,
            candidate,
            &mut candidate_count,
            &mut scaled,
            &mut rejected,
            &mut best_bounded,
        );
        if !hard_constraints.is_empty() {
            let mut constrained_patch = before.clone();
            if row_fit::fit_boundary_samples_with_constraints(
                &mut constrained_patch,
                side,
                samples,
                hard_constraints,
                regularization,
            ) {
                row_fit::snap_boundary_endpoints_to_samples(
                    &mut constrained_patch,
                    side,
                    samples,
                    1.0e-9,
                );
                candidate_count += 1;
                constrained_candidate_count += 1;
                let constrained_candidate = measured_candidate(
                    before,
                    side,
                    samples,
                    hard_constraints,
                    constrained_patch,
                    regularization,
                );
                record_bounded_candidate(
                    before,
                    side,
                    samples,
                    hard_constraints,
                    constrained_candidate,
                    &mut candidate_count,
                    &mut scaled,
                    &mut rejected,
                    &mut best_bounded,
                );
            } else {
                constrained_singular_candidates += 1;
            }
        }
    }
    let mut selected = best_bounded.unwrap_or_else(|| BoundaryFitCandidate {
        patch: before.clone(),
        sample_gap: row_sample_gap(before, side, samples),
        sample_parameter: row_sample_gap_parameter(before, side, samples),
        displacement: 0.0,
        regularization: f64::INFINITY,
        candidate_count,
        constrained_candidate_count: 0,
        constrained_singular_candidates: 0,
        hard_constraint_count: hard_constraints.len(),
        hard_constraint_residual: row_fit::max_boundary_constraint_residual(
            before,
            side,
            hard_constraints,
        ),
        displacement_scaled_candidates: 0,
        displacement_rejected_candidates: 0,
        within_displacement: false,
    });
    selected.candidate_count = candidate_count;
    selected.constrained_candidate_count = constrained_candidate_count;
    selected.constrained_singular_candidates = constrained_singular_candidates;
    selected.displacement_scaled_candidates = scaled;
    selected.displacement_rejected_candidates = rejected;
    selected
}

#[allow(clippy::too_many_arguments)]
fn record_bounded_candidate(
    before: &RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[row_fit::BoundaryRowFitSample],
    hard_constraints: &[row_fit::BoundaryRowFitConstraint],
    candidate: BoundaryFitCandidate,
    candidate_count: &mut usize,
    scaled: &mut usize,
    rejected: &mut usize,
    best_bounded: &mut Option<BoundaryFitCandidate>,
) {
    if candidate.displacement > MAX_BOUNDARY_ROW_DISPLACEMENT {
        *rejected += 1;
        let scaled_patch = row_fit::scale_boundary_delta_to_displacement(
            before,
            &candidate.patch,
            side,
            MAX_BOUNDARY_ROW_DISPLACEMENT,
        );
        *candidate_count += 1;
        *scaled += 1;
        let scaled_candidate = measured_candidate(
            before,
            side,
            samples,
            hard_constraints,
            scaled_patch,
            candidate.regularization,
        );
        if scaled_candidate.displacement > MAX_BOUNDARY_ROW_DISPLACEMENT {
            *rejected += 1;
            return;
        }
        if best_bounded
            .as_ref()
            .is_none_or(|best| boundary_fit_candidate_is_better(&scaled_candidate, best))
        {
            *best_bounded = Some(scaled_candidate);
        }
        return;
    }
    if best_bounded
        .as_ref()
        .is_none_or(|best| boundary_fit_candidate_is_better(&candidate, best))
    {
        *best_bounded = Some(candidate);
    }
}

fn measured_candidate(
    before: &RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[row_fit::BoundaryRowFitSample],
    hard_constraints: &[row_fit::BoundaryRowFitConstraint],
    patch: RectangularGregoryBiNonic,
    regularization: f64,
) -> BoundaryFitCandidate {
    let (sample_gap, sample_parameter) = worst_row_sample_gap(&patch, side, samples);
    let hard_constraint_residual =
        row_fit::max_boundary_constraint_residual(&patch, side, hard_constraints);
    BoundaryFitCandidate {
        sample_gap,
        sample_parameter,
        displacement: row_fit::max_boundary_displacement(before, &patch, side),
        regularization,
        patch,
        candidate_count: 0,
        constrained_candidate_count: 0,
        constrained_singular_candidates: 0,
        hard_constraint_count: hard_constraints.len(),
        hard_constraint_residual,
        displacement_scaled_candidates: 0,
        displacement_rejected_candidates: 0,
        within_displacement: true,
    }
}

pub(super) fn row_sample_gap(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[row_fit::BoundaryRowFitSample],
) -> f64 {
    worst_row_sample_gap(patch, side, samples).0
}

fn row_sample_gap_parameter(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[row_fit::BoundaryRowFitSample],
) -> f64 {
    worst_row_sample_gap(patch, side, samples).1
}

fn worst_row_sample_gap(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[row_fit::BoundaryRowFitSample],
) -> (f64, f64) {
    samples.iter().fold((0.0_f64, 0.0_f64), |worst, sample| {
        let gap = row_fit::eval_boundary(patch, side, sample.parameter).distance_to(&sample.target);
        if gap > worst.0 {
            (gap, sample.parameter)
        } else {
            worst
        }
    })
}

fn boundary_fit_candidate_is_better(
    candidate: &BoundaryFitCandidate,
    current: &BoundaryFitCandidate,
) -> bool {
    (
        candidate.hard_constraint_residual,
        candidate.sample_gap,
        candidate.displacement,
        candidate.regularization,
    ) < (
        current.hard_constraint_residual,
        current.sample_gap,
        current.displacement,
        current.regularization,
    )
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::RectangularGregoryBiNonic;

    use super::*;

    #[test]
    fn bounded_constrained_candidate_honors_interior_vertex_target() {
        let before = RectangularGregoryBiNonic::zeros();
        let samples = (0..10)
            .map(|index| row_fit::BoundaryRowFitSample {
                parameter: index as f64 / 9.0,
                target: Point3::new(0.0, 0.0, 0.0),
            })
            .collect::<Vec<_>>();
        let constraints = [row_fit::BoundaryRowFitConstraint {
            parameter: 0.5,
            target: Point3::new(0.01, 0.0, 0.0),
        }];

        let candidate =
            best_boundary_fit_candidate(&before, BoundarySide::UMin, &samples, &constraints);

        assert!(candidate.constrained_candidate_count > 0);
        assert!(candidate.within_displacement);
        assert!(
            candidate.hard_constraint_residual < 1.0e-7,
            "hard constraint residual stayed {:.3e}",
            candidate.hard_constraint_residual
        );
        assert!(candidate.displacement <= MAX_BOUNDARY_ROW_DISPLACEMENT);
    }
}
