use super::shared_trim_quality::Kp20SupportCellSharedTrimQuality;

const BOUNDARY_SAFE_G0_PREFERENCE_RELATIVE_BAND: f64 = 5.0e-3;
const SHARED_TRIM_G0_TIE_BREAK_EPS: f64 = 1.0e-12;
const SHARED_TRIM_G1_TIE_BREAK_EPS: f64 = 1.0e-12;

#[derive(Debug, Clone, Copy)]
pub(super) struct SupportCellCompletionRank {
    pub(super) boundary_safe: bool,
    pub(super) quality: Kp20SupportCellSharedTrimQuality,
}

pub(super) fn selection_is_better(
    candidate: SupportCellCompletionRank,
    current: Option<SupportCellCompletionRank>,
) -> bool {
    let Some(current) = current else {
        return true;
    };
    let a = candidate.quality;
    let b = current.quality;
    match (a.pcurve_out_of_unit_count, a.degenerate_pcurve_count)
        .cmp(&(b.pcurve_out_of_unit_count, b.degenerate_pcurve_count))
    {
        std::cmp::Ordering::Less => true,
        std::cmp::Ordering::Greater => false,
        std::cmp::Ordering::Equal => g0_tie_break_prefers_candidate(candidate, current),
    }
}

fn g0_tie_break_prefers_candidate(
    candidate: SupportCellCompletionRank,
    current: SupportCellCompletionRank,
) -> bool {
    if boundary_safe_quality_is_equivalent(candidate, current) {
        return true;
    }
    if boundary_safe_quality_is_equivalent(current, candidate) {
        return false;
    }
    if shared_trim_g0_is_equivalent(candidate, current) {
        return candidate.quality.worst_normal_disagreement + SHARED_TRIM_G1_TIE_BREAK_EPS
            < current.quality.worst_normal_disagreement;
    }
    candidate.quality.worst_g0_gap + SHARED_TRIM_G0_TIE_BREAK_EPS < current.quality.worst_g0_gap
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

fn boundary_safe_quality_is_equivalent(
    boundary_safe: SupportCellCompletionRank,
    other: SupportCellCompletionRank,
) -> bool {
    boundary_safe.boundary_safe
        && !other.boundary_safe
        && boundary_safe.quality.worst_g0_gap
            <= other.quality.worst_g0_gap * (1.0 + BOUNDARY_SAFE_G0_PREFERENCE_RELATIVE_BAND)
}
