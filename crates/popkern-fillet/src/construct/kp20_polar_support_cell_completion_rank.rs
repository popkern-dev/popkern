use super::shared_trim_quality::Kp20SupportCellSharedTrimQuality;

const BOUNDARY_SAFE_G0_PREFERENCE_RELATIVE_BAND: f64 = 5.0e-3;
const SHARED_TRIM_G0_TIE_BREAK_EPS: f64 = 1.0e-12;
const SHARED_TRIM_G1_TIE_BREAK_EPS: f64 = 1.0e-12;
const SHARED_TRIM_G2_TIE_BREAK_EPS: f64 = 1.0e-12;

#[derive(Debug, Clone, Copy)]
pub(super) struct SupportCellCompletionRank {
    pub(super) boundary_safe: bool,
    pub(super) quality: Kp20SupportCellSharedTrimQuality,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct SupportCellRankDecision {
    pub(super) candidate_is_better: bool,
    pub(super) reason: &'static str,
}

pub(super) fn selection_decision(
    candidate: SupportCellCompletionRank,
    current: Option<SupportCellCompletionRank>,
) -> SupportCellRankDecision {
    let Some(current) = current else {
        return SupportCellRankDecision {
            candidate_is_better: true,
            reason: "first_accepted_candidate",
        };
    };
    let a = candidate.quality;
    let b = current.quality;
    match (a.pcurve_out_of_unit_count, a.degenerate_pcurve_count)
        .cmp(&(b.pcurve_out_of_unit_count, b.degenerate_pcurve_count))
    {
        std::cmp::Ordering::Less => SupportCellRankDecision {
            candidate_is_better: true,
            reason: "candidate_pcurve_domain_counts_preferred",
        },
        std::cmp::Ordering::Greater => SupportCellRankDecision {
            candidate_is_better: false,
            reason: "current_pcurve_domain_counts_preferred",
        },
        std::cmp::Ordering::Equal => g0_tie_break_decision(candidate, current),
    }
}

fn g0_tie_break_decision(
    candidate: SupportCellCompletionRank,
    current: SupportCellCompletionRank,
) -> SupportCellRankDecision {
    if boundary_safe_quality_is_equivalent(candidate, current) {
        return SupportCellRankDecision {
            candidate_is_better: true,
            reason: "candidate_boundary_safe_quality_preferred",
        };
    }
    if boundary_safe_quality_is_equivalent(current, candidate) {
        return SupportCellRankDecision {
            candidate_is_better: false,
            reason: "current_boundary_safe_quality_preferred",
        };
    }
    if shared_trim_g0_is_equivalent(candidate, current) {
        if shared_trim_g1_is_equivalent(candidate, current) {
            if let Some(candidate_is_better) = commercial_g2_tie_break(candidate, current) {
                if candidate_is_better {
                    return SupportCellRankDecision {
                        candidate_is_better: true,
                        reason: "candidate_commercial_g0_g1_g2_tiebreak_preferred",
                    };
                }
                return SupportCellRankDecision {
                    candidate_is_better: false,
                    reason: "current_commercial_g0_g1_g2_tiebreak_preferred",
                };
            }
        }
        if candidate.quality.worst_normal_disagreement + SHARED_TRIM_G1_TIE_BREAK_EPS
            < current.quality.worst_normal_disagreement
        {
            return SupportCellRankDecision {
                candidate_is_better: true,
                reason: "candidate_commercial_g0_g1_tiebreak_preferred",
            };
        }
        return SupportCellRankDecision {
            candidate_is_better: false,
            reason: "current_commercial_g0_g1_tiebreak_preferred",
        };
    }
    if candidate.quality.worst_g0_gap + SHARED_TRIM_G0_TIE_BREAK_EPS < current.quality.worst_g0_gap
    {
        return SupportCellRankDecision {
            candidate_is_better: true,
            reason: "candidate_g0_preferred",
        };
    }
    SupportCellRankDecision {
        candidate_is_better: false,
        reason: "current_g0_preferred",
    }
}

/// G0 readings are equivalent when numerically tied OR when both
/// are below the commercial G0 target: the gate is met on both
/// sides, so selection falls through to the G1 comparison
/// instead of preferring sub-picometre G0 noise over
/// order-of-magnitude G1 differences.
fn shared_trim_g0_is_equivalent(
    candidate: SupportCellCompletionRank,
    current: SupportCellCompletionRank,
) -> bool {
    if candidate.quality.worst_g0_gap
        <= super::shared_trim_quality::COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET
        && current.quality.worst_g0_gap
            <= super::shared_trim_quality::COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G0_TARGET
    {
        return true;
    }
    (candidate.quality.worst_g0_gap - current.quality.worst_g0_gap).abs()
        <= SHARED_TRIM_G0_TIE_BREAK_EPS
}

/// G1 is equivalent when both readings are already below the
/// commercial G1 target or are numerically tied. At that point
/// selecting on the next continuity order is more meaningful
/// than preserving small below-target G1 differences.
fn shared_trim_g1_is_equivalent(
    candidate: SupportCellCompletionRank,
    current: SupportCellCompletionRank,
) -> bool {
    if candidate.quality.worst_normal_disagreement
        <= super::shared_trim_quality::COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G1_TARGET
        && current.quality.worst_normal_disagreement
            <= super::shared_trim_quality::COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G1_TARGET
    {
        return true;
    }
    (candidate.quality.worst_normal_disagreement - current.quality.worst_normal_disagreement).abs()
        <= SHARED_TRIM_G1_TIE_BREAK_EPS
}

fn commercial_g2_tie_break(
    candidate: SupportCellCompletionRank,
    current: SupportCellCompletionRank,
) -> Option<bool> {
    if candidate.quality.g2_unsupported_sample_count != 0
        || current.quality.g2_unsupported_sample_count != 0
        || candidate.quality.g2_degenerate_sample_count
            != current.quality.g2_degenerate_sample_count
    {
        return None;
    }
    let candidate_g2 = candidate.quality.worst_normal_aligned_g2_residual;
    let current_g2 = current.quality.worst_normal_aligned_g2_residual;
    let candidate_commercial =
        candidate_g2 <= super::shared_trim_quality::COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G2_TARGET;
    let current_commercial =
        current_g2 <= super::shared_trim_quality::COMMERCIAL_SUPPORT_CELL_SHARED_TRIM_G2_TARGET;
    match (candidate_commercial, current_commercial) {
        (true, false) => Some(true),
        (false, true) => Some(false),
        (true, true) => {
            if candidate_g2 + SHARED_TRIM_G2_TIE_BREAK_EPS < current_g2 {
                Some(true)
            } else if current_g2 + SHARED_TRIM_G2_TIE_BREAK_EPS < candidate_g2 {
                Some(false)
            } else {
                None
            }
        }
        (false, false) => None,
    }
}

fn boundary_safe_quality_is_equivalent(
    boundary_safe: SupportCellCompletionRank,
    other: SupportCellCompletionRank,
) -> bool {
    boundary_safe.boundary_safe
        && !other.boundary_safe
        && boundary_safe.quality.worst_g0_gap
            <= other.quality.worst_g0_gap * (1.0 + BOUNDARY_SAFE_G0_PREFERENCE_RELATIVE_BAND)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::construct::kp20_polar_support_cell_shared_trim_quality::{
        SUPPORT_CELL_TOP_G1_SAMPLE_COUNT, SUPPORT_CELL_TOP_G2_SAMPLE_COUNT,
    };

    fn rank(g1: f64, normal_aligned_g2: f64) -> SupportCellCompletionRank {
        SupportCellCompletionRank {
            boundary_safe: true,
            quality: Kp20SupportCellSharedTrimQuality {
                shared_trim_count: 1,
                partial_shared_trim_count: 0,
                pcurve_out_of_unit_count: 0,
                degenerate_pcurve_count: 0,
                max_boundary_row_shared_trim_count: 1,
                worst_abs_uv: 1.0,
                worst_g0_gap: 0.0,
                worst_normal_disagreement: g1,
                worst_g2_residual: normal_aligned_g2,
                worst_normal_aligned_g2_residual: normal_aligned_g2,
                g2_unsupported_sample_count: 0,
                g2_degenerate_sample_count: 0,
                worst_trim: None,
                worst_g1_trim: None,
                top_g1_samples: [None; SUPPORT_CELL_TOP_G1_SAMPLE_COUNT],
                top_g2_samples: [None; SUPPORT_CELL_TOP_G2_SAMPLE_COUNT],
                top_normal_aligned_g2_samples: [None; SUPPORT_CELL_TOP_G2_SAMPLE_COUNT],
            },
        }
    }

    #[test]
    fn g2_breaks_tie_after_g0_and_g1_are_commercially_met() {
        let decision = selection_decision(rank(1.0e-3, 4.0e-2), Some(rank(2.0e-4, 6.0e-2)));

        assert!(decision.candidate_is_better);
        assert_eq!(
            decision.reason,
            "candidate_commercial_g0_g1_g2_tiebreak_preferred"
        );
    }

    #[test]
    fn g1_remains_preferred_when_candidate_exceeds_commercial_g1() {
        let decision = selection_decision(rank(2.0e-3, 1.0e-2), Some(rank(1.0e-3, 4.0e-2)));

        assert!(!decision.candidate_is_better);
        assert_eq!(
            decision.reason,
            "current_commercial_g0_g1_tiebreak_preferred"
        );
    }

    #[test]
    fn g1_remains_preferred_when_both_g2_values_exceed_commercial_target() {
        let decision = selection_decision(rank(1.0e-3, 5.5e-2), Some(rank(2.0e-4, 6.0e-2)));

        assert!(!decision.candidate_is_better);
        assert_eq!(
            decision.reason,
            "current_commercial_g0_g1_tiebreak_preferred"
        );
    }
}
