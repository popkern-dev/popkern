//! Operation-status formatting for chain constructive re-authoring.

use super::super::cell_reauthor::CellReauthorOutcome;

pub(super) fn reauthor_detail(outcome: &CellReauthorOutcome) -> String {
    format!(
        ";derived_row_count={};leg_solve_declined_count={};field_declined_sub_chain_count={};twist_conflict_count={};max_twist_conflict_distance={:.6e};max_derived_span_residual={:.6e};worst_span_residual_row={};repaired_leg_count={};max_leg_out_of_plane_ratio={:.6e};max_twist_conflict_normal_component={:.6e};max_relative_twist_displacement={:.6e};refilled_cell_count={};predicted_worst_trim_g1={:.6e};cross_field_sample_count={};near_tangent_cross_count={};beta_sign_flip_count={};max_tangential_fraction={:.6e};min_beta_magnitude_fraction={:.6e};worst_cross_field_detail=[{}];beta_envelope_sample_count={};beta_envelope_constraint_count={};beta_envelope_repaired_row_count={};max_beta_envelope_control_delta={:.6e};min_beta_envelope_fraction_before={:.6e};min_beta_envelope_fraction_after={:.6e};pair_normal_sample_count={};pair_normal_constraint_count={};pair_normal_repaired_row_count={};pair_normal_limited_row_count={};max_pair_normal_control_delta={:.6e};max_pair_normal_disagreement_before={:.6e};max_pair_normal_disagreement_after={:.6e};worst_pair_normal_detail=[{}];clamped_leg_beta_count={};one_sided_band_count={};dual_multi_demand_side_count={};dual_twist_pair_count={};max_dual_twist_pair_separation={:.6e};predicted_worst_detail=[{}]",
        outcome.derived_row_count,
        outcome.leg_solve_declined_count,
        outcome.field_declined_sub_chain_count,
        outcome.twist_conflict_count,
        outcome.max_twist_conflict_distance,
        outcome.max_derived_span_residual,
        outcome.worst_span_residual_row.as_deref().unwrap_or("none"),
        outcome.repaired_leg_count,
        outcome.max_leg_out_of_plane_ratio,
        outcome.max_twist_conflict_normal_component,
        outcome.max_relative_twist_displacement,
        outcome.refilled_cell_count,
        outcome.predicted_worst_trim_g1,
        outcome.cross_field_sample_count,
        outcome.near_tangent_cross_count,
        outcome.beta_sign_flip_count,
        outcome.max_tangential_fraction,
        outcome.min_beta_magnitude_fraction,
        outcome
            .worst_cross_field_detail
            .as_deref()
            .unwrap_or("none"),
        outcome.beta_envelope_sample_count,
        outcome.beta_envelope_constraint_count,
        outcome.beta_envelope_repaired_row_count,
        outcome.max_beta_envelope_control_delta,
        outcome.min_beta_envelope_fraction_before,
        outcome.min_beta_envelope_fraction_after,
        outcome.pair_normal_sample_count,
        outcome.pair_normal_constraint_count,
        outcome.pair_normal_repaired_row_count,
        outcome.pair_normal_limited_row_count,
        outcome.max_pair_normal_control_delta,
        outcome.max_pair_normal_disagreement_before,
        outcome.max_pair_normal_disagreement_after,
        outcome
            .worst_pair_normal_detail
            .as_deref()
            .unwrap_or("none"),
        outcome.clamped_leg_beta_count,
        outcome.one_sided_band_count,
        outcome.dual_multi_demand_side_count,
        outcome.dual_twist_pair_count,
        outcome.max_dual_twist_pair_separation,
        outcome.predicted_worst_detail.as_deref().unwrap_or("none")
    )
}
