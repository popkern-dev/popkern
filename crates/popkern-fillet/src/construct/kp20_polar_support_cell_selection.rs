use super::Kp20SupportCellPrimitiveTrial;

const TOPOLOGY_PATH_TIE_TOL: f64 = 1.0e-12;
const HIGHER_RESOLUTION_PROMOTION_EPS: f64 = 1.0e-9;

pub(super) fn should_replace_selected_support_cell_trial(
    candidate: &Kp20SupportCellPrimitiveTrial,
    selected: &Kp20SupportCellPrimitiveTrial,
) -> bool {
    if candidate.folded_patch_count != selected.folded_patch_count {
        return candidate.folded_patch_count < selected.folded_patch_count;
    }
    let candidate_tier = support_cell_resolution_tier(candidate);
    let selected_tier = support_cell_resolution_tier(selected);
    if candidate_tier != selected_tier {
        return match (candidate_tier, selected_tier) {
            (SupportCellResolutionTier::Unsupported, _) => false,
            (_, SupportCellResolutionTier::Unsupported) => true,
            (SupportCellResolutionTier::Eight, SupportCellResolutionTier::Four) => {
                candidate.topology_path_g0 + HIGHER_RESOLUTION_PROMOTION_EPS
                    < selected.topology_path_g0
            }
            (SupportCellResolutionTier::Four, SupportCellResolutionTier::Eight) => {
                candidate.topology_path_g0
                    <= selected.topology_path_g0 + HIGHER_RESOLUTION_PROMOTION_EPS
            }
            _ => false,
        };
    }
    candidate.topology_path_g0 + TOPOLOGY_PATH_TIE_TOL < selected.topology_path_g0
        || ((candidate.topology_path_g0 - selected.topology_path_g0).abs() <= TOPOLOGY_PATH_TIE_TOL
            && candidate.patch_count < selected.patch_count)
}

pub(super) fn selection_policy_label(trial: &Kp20SupportCellPrimitiveTrial) -> &'static str {
    match support_cell_resolution_tier(trial) {
        SupportCellResolutionTier::Four => "committed_four_cell_floor",
        SupportCellResolutionTier::Eight => "higher_resolution_eight_cell_promoted_by_g0",
        SupportCellResolutionTier::Unsupported => "unsupported_measurement_only_inventory",
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SupportCellResolutionTier {
    Unsupported,
    Four,
    Eight,
}

fn support_cell_resolution_tier(
    trial: &Kp20SupportCellPrimitiveTrial,
) -> SupportCellResolutionTier {
    match trial.cells_u * trial.cells_v {
        4 => SupportCellResolutionTier::Four,
        8 => SupportCellResolutionTier::Eight,
        _ => SupportCellResolutionTier::Unsupported,
    }
}
