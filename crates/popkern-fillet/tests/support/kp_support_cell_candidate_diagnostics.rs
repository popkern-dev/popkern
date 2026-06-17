pub(crate) fn assert_candidate_diagnostics(
    fixture: &str,
    support_cell_dispatch_details: &[String],
) {
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("repair_status=accepted_non_regressing_final_repair")
                && detail.contains("selected_candidate=")
        }),
        "{fixture} support-cell route should emit accepted transactional shared-trim repair diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details
            .iter()
            .any(|detail| detail.contains("candidate_status=accepted_selected")),
        "{fixture} support-cell selector should emit per-candidate score diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("candidate_status=accepted_selected")
                && detail.contains("rank_reason=selected_best_rank")
        }),
        "{fixture} selected support-cell candidate should report the final rank reason; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("candidate_status=accepted_not_selected")
                && detail.contains("rank_reason=current_")
        }),
        "{fixture} non-selected viable support-cell candidates should report why the selected rank won; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details
            .iter()
            .filter(|detail| detail.contains("candidate_status="))
            .all(|detail| detail.contains("rank_reason=")),
        "{fixture} every support-cell candidate score should carry a structured rank reason; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details
            .iter()
            .filter(|detail| detail.contains("candidate_status="))
            .all(|detail| detail.contains("after_quality_status=")
                && detail.contains("after_shared_trim_count=")
                && detail.contains("after_full_row_owned_trim_count=")
                && detail.contains("after_segment_owned_trim_count=")
                && detail.contains("after_max_boundary_row_shared_trim_count=")),
        "{fixture} every support-cell candidate score should expose quality classification and ownership counts; support_cell_details={support_cell_dispatch_details:?}"
    );
    for candidate in [
        "pcurve_then_boundary_if_g0_safe",
        "pcurve_then_worst_boundary_if_g0_safe",
        "pcurve_then_top3_boundaries_if_g0_safe",
        "pcurve_then_ranked_boundary_prefix_if_g0_safe",
        "pcurve_then_worst_boundary_neighborhood_if_g0_safe",
        "pcurve_then_worst_endpoint_cluster_if_g0_safe",
        "pcurve_then_ranked_boundary_prefix_then_worst_endpoint_cluster_if_g0_safe",
        "pcurve_then_ranked_boundary_endpoint_cluster_then_row_group_if_g0_safe",
        "pcurve_then_ranked_boundary_endpoint_cluster_then_two_row_groups_if_g0_safe",
        "pcurve_then_ranked_boundary_endpoint_cluster_then_three_row_groups_if_g0_safe",
        "pcurve_then_ranked_boundary_endpoint_row_group_then_tangent_if_g0_safe",
        "pcurve_then_ranked_boundary_endpoint_row_group_then_exact_boundary_if_g0_safe",
        "pcurve_then_ranked_boundary_endpoint_row_group_then_partial_boundary_if_g0_safe",
        "pcurve_then_ranked_boundary_endpoint_row_group_then_exact_and_partial_boundary_if_g0_safe",
        "pcurve_then_ranked_boundary_endpoint_row_group_then_exact_partial_tangent_if_g0_safe",
        "pcurve_then_ranked_boundary_endpoint_row_group_then_exact_partial_top_g1_normal_solve_if_g0_safe",
        "pcurve_then_chain_curve_network_if_g0_safe",
        "pcurve_then_chain_network_cross_tangent_ribbon_if_g0_safe",
        "pcurve_then_chain_constructive_cells_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_cells_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_joint_leg_solve_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_beta_envelope_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_pair_normal_solve_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_active_normal_solve_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_beta_envelope_active_normal_solve_if_g0_safe",
        "pcurve_then_exact_partial_owner_row_chain_constructive_dual_twist_beta_envelope_active_normal_solve_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_if_g0_safe",
        "pcurve_then_exact_partial_owner_row_chain_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_if_g0_safe",
        "pcurve_then_exact_partial_owner_row_chain_constructive_dual_twist_beta_envelope_active_normal_all_derived_row_band_curvature_solve_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_worst_g1_tangent_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_worst_g1_normal_solve_if_g0_safe",
        "pcurve_then_chain_constructive_dual_twist_endpoint_adjacent_worst_g1_normal_solve_if_g0_safe",
        "pcurve_then_worst_boundary_row_group_if_g0_safe",
    ] {
        assert!(
            support_cell_dispatch_details
                .iter()
                .any(|detail| detail.contains(&format!("candidate={candidate}"))),
            "{fixture} support-cell selector should score {candidate}; support_cell_details={support_cell_dispatch_details:?}"
        );
    }
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("candidate=pcurve_then_chain_curve_network_if_g0_safe")
                && detail.contains("chain_curve_network")
                && detail.contains("sub_chain_count=")
        }),
        "{fixture} chain-curve network candidate should report its sub-chain inventory; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("candidate=pcurve_then_chain_network_cross_tangent_ribbon_if_g0_safe")
                && detail.contains("chain_ribbon_network")
                && detail.contains("projected_row_count=")
        }),
        "{fixture} chain ribbon candidate should report its projected-row inventory; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("candidate=pcurve_then_chain_constructive_cells_if_g0_safe")
                && detail.contains("chain_constructive_cells")
                && detail.contains("derived_row_count=")
        }),
        "{fixture} constructive-cells candidate should report its derived-row inventory; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("candidate=pcurve_then_chain_constructive_dual_twist_cells_if_g0_safe")
                && detail.contains("chain_constructive_dual_twist_cells")
                && detail.contains("dual_twist_pair_count=")
        }),
        "{fixture} dual-twist constructive candidate should report its pair inventory; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_chain_constructive_dual_twist_beta_envelope_active_normal_solve_if_g0_safe",
            ) && detail.contains("chain_constructive_dual_twist_beta_envelope_active_normal_solve")
                && detail.contains("beta_envelope_constraint_count=")
                && detail.contains("active_normal_constraint_count=")
        }),
        "{fixture} combined beta-envelope active-normal candidate should report both solve inventories; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_exact_partial_owner_row_chain_constructive_dual_twist_beta_envelope_active_normal_solve_if_g0_safe",
            ) && detail.contains("chain_owner_row_constructive_dual_twist_beta_envelope_active_normal_solve")
                && detail.contains("owner_row_source_sub_chain_count=")
                && detail.contains("fit_source_sub_chain_count=")
                && detail.contains("active_normal_constraint_count=")
        }),
        "{fixture} owner-row chain candidate should report source ownership and active-normal diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_chain_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_if_g0_safe",
            ) && detail.contains(
                "chain_constructive_dual_twist_beta_envelope_active_normal_curvature_solve",
            ) && detail.contains("curvature_solve_status=")
                && detail.contains("curvature_weighted_normal_solve_status=")
                && detail.contains("max_curvature_control_delta=")
        }),
        "{fixture} curvature candidate should report constrained G2 solve diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_exact_partial_owner_row_chain_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_if_g0_safe",
            ) && detail.contains(
                "chain_owner_row_constructive_dual_twist_beta_envelope_active_normal_curvature_solve",
            ) && detail.contains("owner_row_source_sub_chain_count=")
                && detail.contains("fit_source_sub_chain_count=")
                && detail.contains("curvature_solve_status=")
                && detail.contains("max_curvature_control_delta=")
        }),
        "{fixture} owner-row curvature candidate should report source ownership and constrained G2 diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_exact_partial_owner_row_chain_constructive_dual_twist_beta_envelope_active_normal_all_derived_row_band_curvature_solve_if_g0_safe",
            ) && detail.contains(
                "chain_owner_row_constructive_dual_twist_beta_envelope_active_normal_all_derived_row_band_curvature_solve",
            ) && detail.contains("owner_row_source_sub_chain_count=")
                && detail.contains("fit_source_sub_chain_count=")
                && detail.contains("curvature_solve_scope=all_derived_row_bands")
                && detail.contains("curvature_solve_budget=expanded_row_band")
                && detail.contains("max_curvature_control_delta=")
        }),
        "{fixture} owner-row all-derived curvature candidate should report source ownership, all-row scope, and constrained G2 diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("top_normal_aligned_g2_samples=[")
                && detail.contains("a_proj=(")
                && detail.contains("a_ii=(")
                && detail.contains("a_terms=(")
                && detail.contains("b_proj=(")
                && detail.contains("b_ii=(")
                && detail.contains("b_terms=(")
                && detail.contains("q_gap=")
                && detail.contains("normal_only=(")
        }),
        "{fixture} support-cell G2 diagnostics should expose projected-direction and second-form decomposition for top normal-aligned samples; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert_candidate_stage_progression(fixture, support_cell_dispatch_details);
}

fn assert_candidate_stage_progression(fixture: &str, support_cell_dispatch_details: &[String]) {
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("candidate=pcurve_then_top3_boundaries_if_g0_safe")
                && detail.contains("candidate_status=accepted")
                && detail.contains("operation_status=boundary_kept")
        }),
        "{fixture} targeted top-3 boundary candidate should remain accepted and boundary-safe; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("candidate=pcurve_then_ranked_boundary_prefix_if_g0_safe")
                && detail.contains("operation_status=boundary_kept")
        }),
        "{fixture} ranked-prefix candidate should remain an accepted measured support-cell repair pass; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains("candidate_status=accepted_selected")
                && detail.contains(
                    "candidate=pcurve_then_chain_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_if_g0_safe",
                )
                && detail.contains(
                    "operation_status=chain_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_authored",
                )
        }),
        "{fixture} should select the combined beta-envelope active-normal curvature dual-twist chain-cell repair; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_ranked_boundary_endpoint_cluster_then_two_row_groups_if_g0_safe",
            ) && detail.contains("operation_status=ranked_prefix_")
        }),
        "{fixture} should score the bounded two-row-group composed candidate; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_ranked_boundary_endpoint_cluster_then_three_row_groups_if_g0_safe",
            ) && detail.contains("operation_status=ranked_prefix_")
        }),
        "{fixture} should score the bounded three-row-group composed candidate; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_ranked_boundary_endpoint_row_group_then_exact_boundary_if_g0_safe",
            ) && detail.contains("operation_status=ranked_prefix_endpoint_")
                && detail.contains("exact_boundary_rows_authored")
        }),
        "{fixture} should score the exact shared-boundary row candidate as a real construction pass; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_ranked_boundary_endpoint_row_group_then_partial_boundary_if_g0_safe",
            ) && detail.contains("segment_topology_already_segmented")
                && detail.contains("bounded_partial_boundary_segments_fit")
        }),
        "{fixture} should score the partial shared-boundary segment candidate as a real construction pass; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_ranked_boundary_endpoint_row_group_then_exact_and_partial_boundary_if_g0_safe",
            ) && detail.contains("exact_and_partial_boundary_")
        }),
        "{fixture} should score the exact-plus-partial shared-boundary candidate as a real construction pass; support_cell_details={support_cell_dispatch_details:?}"
    );
    // The post-partial tangent candidate must stay a scored,
    // accepted measured pass. Its internal stage trajectory
    // (which sub-operations keep or decline) legitimately varies
    // with the surrounding selection — since the constructive
    // chain candidates win, this candidate's tangent stage no
    // longer carries the repair and may decline.
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_ranked_boundary_endpoint_row_group_then_exact_partial_tangent_if_g0_safe",
            ) && detail.contains("candidate_status=accepted")
        }),
        "{fixture} should score the post-partial tangent candidate as an accepted measured pass; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| {
            detail.contains(
                "candidate=pcurve_then_ranked_boundary_endpoint_row_group_then_exact_partial_top_g1_normal_solve_if_g0_safe",
            ) && detail.contains("exact_partial_top_g1_normal_solve")
                && detail.contains("worst_g1_normal_solve_status=")
        }),
        "{fixture} should score the exact/partial constrained normal solve with structured diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
}
