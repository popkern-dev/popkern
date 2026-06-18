//! Route-level support-cell selection diagnostics for KP20 dispatch.

use crate::types::ApexConstructionDiagnostic;

pub(super) fn support_cell_primary_macro_bypass_fact(
    n: usize,
    accepted_primary: &'static str,
    accepted_primary_detail: String,
) -> ApexConstructionDiagnostic {
    support_cell_nonselection_fact(
        n,
        format!(
            "evaluated=0;target_count=0;patch_count=0;folded_patch_count=0;cells_u=0;cells_v=0;topology_path_g0=inf;selection_policy=support_cell_not_evaluated;construction_mode=none;target_piece=NearMidpointLeft;committed_floor_patch_count=0;decline_reason=support_cell_primitive_not_evaluated_primary_macro_accepted;accepted_primary={accepted_primary};accepted_primary_detail={accepted_primary_detail}",
        ),
    )
}

pub(super) fn support_cell_candidate_nonselection_fact(
    n: usize,
    lifted: &super::limited_lift::Kp20LimitedLiftCandidate,
) -> Option<ApexConstructionDiagnostic> {
    let topology = lifted.expanded_support_topology.as_ref()?;
    let cells = &topology.support_cell_primitive;
    if cells.evaluated == 0 {
        return None;
    }
    Some(support_cell_nonselection_fact(
        n,
        format!(
            "evaluated={};target_count={};patch_count={};folded_patch_count={};cells_u={};cells_v={};topology_path_g0={:.6e};selection_policy={};construction_mode={};target_piece={:?};committed_floor_patch_count={};decline_reason={}",
            cells.evaluated,
            cells.target_count,
            cells.patch_count,
            cells.folded_patch_count,
            cells.cells_u,
            cells.cells_v,
            cells.topology_path_g0,
            cells.selection_policy,
            cells.construction_mode,
            cells.target_piece,
            cells.committed_floor_patches.len(),
            cells.decline_reason,
        ),
    ))
}

fn support_cell_nonselection_fact(n: usize, detail: String) -> ApexConstructionDiagnostic {
    super::super::kp20_polar_diagnostics::dispatch_fact(
        "support_cell_candidate_nonselection",
        n,
        Some(detail),
        None,
        None,
    )
}
