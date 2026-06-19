//! Displacement-constrained row-fit selection for external KP20 cap/support seams.

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;

use super::{BoundarySide, row_fit};

pub(super) const MAX_BOUNDARY_ROW_DISPLACEMENT: f64 = 0.1;
const ROW_FIT_REGULARIZATION_CANDIDATES: [f64; 4] = [1.0e-4, 1.0e-6, 1.0e-8, 0.0];

pub(super) struct BoundaryFitCandidate {
    pub(super) patch: RectangularGregoryBiNonic,
    pub(super) sample_gap: f64,
    pub(super) displacement: f64,
    pub(super) regularization: f64,
    pub(super) candidate_count: usize,
    pub(super) displacement_rejected_candidates: usize,
    pub(super) within_displacement: bool,
}

pub(super) fn best_boundary_fit_candidate(
    before: &RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[row_fit::BoundaryRowFitSample],
) -> BoundaryFitCandidate {
    let mut best_bounded = None;
    let mut rejected = 0;
    for regularization in ROW_FIT_REGULARIZATION_CANDIDATES {
        let mut patch = before.clone();
        row_fit::fit_boundary_samples_with_regularization(
            &mut patch,
            side,
            samples,
            regularization,
        );
        let candidate = BoundaryFitCandidate {
            sample_gap: row_sample_gap(&patch, side, samples),
            displacement: row_fit::max_boundary_displacement(before, &patch, side),
            regularization,
            patch,
            candidate_count: ROW_FIT_REGULARIZATION_CANDIDATES.len(),
            displacement_rejected_candidates: 0,
            within_displacement: true,
        };
        if candidate.displacement > MAX_BOUNDARY_ROW_DISPLACEMENT {
            rejected += 1;
            continue;
        }
        if best_bounded
            .as_ref()
            .is_none_or(|best| boundary_fit_candidate_is_better(&candidate, best))
        {
            best_bounded = Some(candidate);
        }
    }
    let mut selected = best_bounded.unwrap_or_else(|| BoundaryFitCandidate {
        patch: before.clone(),
        sample_gap: row_sample_gap(before, side, samples),
        displacement: 0.0,
        regularization: f64::INFINITY,
        candidate_count: ROW_FIT_REGULARIZATION_CANDIDATES.len(),
        displacement_rejected_candidates: 0,
        within_displacement: false,
    });
    selected.displacement_rejected_candidates = rejected;
    selected
}

pub(super) fn row_sample_gap(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[row_fit::BoundaryRowFitSample],
) -> f64 {
    samples.iter().fold(0.0_f64, |worst, sample| {
        worst.max(row_fit::eval_boundary(patch, side, sample.parameter).distance_to(&sample.target))
    })
}

fn boundary_fit_candidate_is_better(
    candidate: &BoundaryFitCandidate,
    current: &BoundaryFitCandidate,
) -> bool {
    (
        candidate.sample_gap,
        candidate.displacement,
        candidate.regularization,
    ) < (
        current.sample_gap,
        current.displacement,
        current.regularization,
    )
}
