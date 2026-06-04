//! Structured diagnostic labels for KP20 multipatch candidates.

pub(super) fn multipatch_validation_label(
    validation: &super::kp20_polar_multipatch_validate::Kp20MultiPatchValidationReport,
) -> String {
    format!(
        "sectors={},strips={},interfaces={},inners={},finite_cps={},invalid_cps={},boundary_constraints={},adjacent_constraints={},interior_constraints={},exact_endpoint_adjacent={},strip_g0={:.3e},strip_adjacent={:.3e},inner_g0={:.3e},inner_adjacent={:.3e},inner_core={:.3e},decline={}",
        validation.validated_sector_count,
        validation.validated_strip_count,
        validation.validated_interface_count,
        validation.validated_inner_patch_count,
        validation.finite_control_point_count,
        validation.invalid_control_point_count,
        validation.validated_boundary_constraint_count,
        validation.validated_adjacent_constraint_count,
        validation.validated_interior_constraint_count,
        validation.exact_boundary_owned_endpoint_count,
        validation.max_strip_interface_g0_residual,
        validation.max_strip_adjacent_residual,
        validation.max_inner_boundary_g0_residual,
        validation.max_inner_adjacent_residual,
        validation.max_inner_core_residual,
        validation.decline_reason.unwrap_or("none"),
    )
}

pub(super) fn multipatch_evaluation_label(
    evaluation: &super::kp20_polar_multipatch_evaluate::Kp20MultiPatchEvaluationReport,
) -> String {
    format!(
        "attempted_strips={},attempted_inners={},changed_strips={},changed_inners={},changed_pieces={},strip_optimizer_trials={},strip_optimizer_mode={},strip_optimizer_scale={},inner_construction_mode={},fast_evaluated={},exact_evaluated={},fast_folds={},exact_folds={},fast_support_folds={},fast_inner_folds={},exact_support_folds={},exact_inner_folds={},exact_fold_assessment={},fast_piece_folds={},exact_piece_folds={},fast_min_cos={},exact_min_cos={},fast_worst_site={},exact_worst_site={},fast_fold_zones={},exact_fold_zones={},fast_validation={},exact_validation={},best_rejected_exact_guided={},validation_clean={},accepted={},decline={}",
        evaluation.attempted_strip_patch_count,
        evaluation.attempted_inner_patch_count,
        evaluation.changed_strip_piece_count,
        evaluation.changed_inner_piece_count,
        evaluation.changed_piece_count,
        evaluation.strip_optimizer_trial_count,
        evaluation.strip_optimizer_mode.unwrap_or("none"),
        optional_f64_label(evaluation.strip_optimizer_selected_scale),
        evaluation.inner_construction_mode.unwrap_or("none"),
        evaluation.fast_evaluated,
        evaluation.exact_evaluated,
        optional_usize_label(evaluation.fast_fold_count),
        optional_usize_label(evaluation.exact_fold_count),
        optional_usize_label(evaluation.fast_support_fold_count),
        optional_usize_label(evaluation.fast_inner_fold_count),
        optional_usize_label(evaluation.exact_support_fold_count),
        optional_usize_label(evaluation.exact_inner_fold_count),
        exact_fold_assessment(evaluation),
        optional_piece_counts_label(evaluation.fast_piece_fold_counts),
        optional_piece_counts_label(evaluation.exact_piece_fold_counts),
        optional_f64_label(evaluation.fast_min_jacobian_normal_cos),
        optional_f64_label(evaluation.exact_min_jacobian_normal_cos),
        optional_fold_site_label(evaluation.fast_worst_fold_site),
        optional_fold_site_label(evaluation.exact_worst_fold_site),
        optional_fold_zone_counts_label(evaluation.fast_fold_zone_counts),
        optional_fold_zone_counts_label(evaluation.exact_fold_zone_counts),
        optional_validation_label(evaluation.fast_validation),
        optional_validation_label(evaluation.exact_validation),
        optional_rejected_exact_guided_label(evaluation.best_rejected_exact_guided),
        evaluation.validation_clean,
        evaluation.accepted,
        evaluation.decline_reason.unwrap_or("none"),
    )
}

fn exact_fold_assessment(
    evaluation: &super::kp20_polar_multipatch_evaluate::Kp20MultiPatchEvaluationReport,
) -> &'static str {
    match (
        evaluation.exact_fold_count,
        evaluation.exact_support_fold_count,
        evaluation.exact_inner_fold_count,
    ) {
        (None, _, _) => "not_exact_evaluated",
        (Some(0), _, _) => "no_exact_folds",
        (Some(_), Some(support), Some(0)) if support > 0 => "support_exact_folds",
        (Some(_), Some(0), Some(inner)) if inner > 0 => "inner_exact_folds",
        (Some(_), Some(support), Some(inner)) if support > 0 && inner > 0 => "mixed_exact_folds",
        _ => "unknown_exact_fold_distribution",
    }
}

fn optional_rejected_exact_guided_label(
    rejected: Option<super::kp20_polar_multipatch_evaluate::Kp20RejectedExactGuidedReport>,
) -> String {
    rejected.map_or_else(
        || "none".to_string(),
        |rejected| {
            format!(
                "strip:{}|inner:{}|scale:{:.3e}|fast_folds:{}|exact_folds:{}|exact_min_cos:{:.3e}|exact_worst_site:{}|exact_fold_zones:{}|exact_validation:{}",
                rejected.strip_mode,
                rejected.inner_mode,
                rejected.scale,
                rejected.fast_fold_count,
                rejected.exact_fold_count,
                rejected.exact_min_jacobian_normal_cos,
                optional_fold_site_label(rejected.exact_worst_fold_site),
                fold_zone_counts_label(rejected.exact_fold_zone_counts),
                optional_validation_label(Some(rejected.exact_validation)),
            )
        },
    )
}

fn optional_fold_zone_counts_label(
    counts: Option<super::kp20_polar_fold_summary::Kp20FoldZoneCounts>,
) -> String {
    counts.map_or_else(|| "none".to_string(), fold_zone_counts_label)
}

fn fold_zone_counts_label(counts: super::kp20_polar_fold_summary::Kp20FoldZoneCounts) -> String {
    format!(
        "corner:{}|edge:{}|interior:{}",
        counts.origin_corner, counts.edge_shoulder, counts.interior
    )
}

fn optional_fold_site_label(site: Option<super::kp20_polar_fold_summary::Kp20FoldSite>) -> String {
    site.map_or_else(
        || "none".to_string(),
        |site| {
            format!(
                "sector:{}|piece:{:?}|u:{:.3}|v:{:.3}|ratio:{:.3e}|cos:{:.3e}",
                site.sector,
                site.piece,
                site.u,
                site.v,
                site.min_jacobian_ratio,
                site.min_jacobian_normal_cos,
            )
        },
    )
}

fn optional_validation_label(
    residuals: Option<super::kp20_polar_multipatch_evaluate::Kp20CandidateValidationResiduals>,
) -> String {
    residuals.map_or_else(
        || "none".to_string(),
        |residuals| {
            format!(
                "tensor:{:.3e}|c2_pos:{:.3e}|c2_first:{:.3e}|origin_g2:{:.3e}",
                residuals.tensor_border,
                residuals.internal_c2_position,
                residuals.internal_c2_first,
                residuals.origin_g2_strip,
            )
        },
    )
}

pub(super) fn multipatch_inner_label(
    inner: &super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport,
) -> String {
    format!(
        "attempted={},built={},finite_cps={},invalid_cps={},boundary_constraints={},adjacent_constraints={},interior_constraints={},exact_endpoint_adjacent={},max_g0={:.3e},max_adjacent={:.3e},max_inner_core={:.3e},min_adjacent_step={:.3e},patches={},decline={},names={}",
        inner.attempted_inner_patch_count,
        inner.built_inner_patch_count,
        inner.finite_control_point_count,
        inner.invalid_control_point_count,
        inner.boundary_constraint_count,
        inner.adjacent_row_constraint_count,
        inner.interior_constraint_count,
        inner.exact_boundary_owned_endpoint_count,
        inner.max_boundary_g0_residual,
        inner.max_adjacent_row_residual,
        inner.max_interior_residual,
        inner.min_adjacent_row_step,
        inner_patch_anchor_label(inner),
        inner.decline_reason.unwrap_or("none"),
        if inner.inner_patch_names.is_empty() {
            "none".to_string()
        } else {
            inner.inner_patch_names.join(",")
        },
    )
}

pub(super) fn multipatch_interfaces_label(
    interfaces: &super::kp20_polar_multipatch_interface::Kp20MultiPatchInterfaceBuildReport,
) -> String {
    format!(
        "attempted={},built={},finite_cps={},invalid_cps={},anchors={},decline={},names={}",
        interfaces.attempted_interface_count,
        interfaces.built_interface_count,
        interfaces.finite_control_point_count,
        interfaces.invalid_control_point_count,
        interface_anchor_label(interfaces),
        interfaces.decline_reason.unwrap_or("none"),
        if interfaces.interface_names.is_empty() {
            "none".to_string()
        } else {
            interfaces.interface_names.join(",")
        },
    )
}

pub(super) fn multipatch_strips_label(
    strips: &super::kp20_polar_multipatch_strip::Kp20MultiPatchStripBuildReport,
) -> String {
    format!(
        "attempted={},built={},reauthored={},finite_cps={},invalid_cps={},max_interior_displacement={:.3e},min_parameter_step={:.3e},max_parameter_skew={:.3e},max_tangent_step_imbalance={:.3e},min_tangent_alignment={:.3e},row_strategy_columns=uniform:{}|fold_site:{}|tangent_balance:{}|combined:{},row_strategy_assessment={},tangent_row_reauthored={},max_tangent_row_displacement={:.3e},max_row_spacing_ratio={:.3e},max_row_curvature_ratio={:.3e},strip_cps={},surface_count={},surface_sectors={},decline={},names={}",
        strips.attempted_strip_patch_count,
        strips.built_strip_patch_count,
        strips.reauthored_strip_patch_count,
        strips.finite_control_point_count,
        strips.invalid_control_point_count,
        strips.max_reauthored_interior_displacement,
        strips.min_parameter_step,
        strips.max_parameter_skew,
        strips.max_tangent_step_imbalance,
        strips.min_tangent_alignment,
        strips.uniform_row_strategy_columns,
        strips.fold_site_row_strategy_columns,
        strips.tangent_balance_row_strategy_columns,
        strips.combined_row_strategy_columns,
        strips.row_strategy_assessment,
        strips.tangent_row_reauthored_columns,
        strips.max_tangent_row_displacement,
        strips.max_row_spacing_ratio,
        strips.max_row_curvature_ratio,
        strips.strip_control_point_count,
        strips.strip_surfaces.len(),
        strip_surface_sector_label(strips),
        strips.decline_reason.unwrap_or("none"),
        if strips.strip_names.is_empty() {
            "none".to_string()
        } else {
            strips.strip_names.join(",")
        },
    )
}

pub(super) fn multipatch_inventory_label(
    inventory: &super::kp20_polar_multipatch_inventory::Kp20MultiPatchInventoryReport,
) -> String {
    let sectors = inventory
        .sector_reports
        .iter()
        .map(|sector| {
            format!(
                "{}:strips=[{}],inner={},interfaces=[{}]",
                sector.sector_name,
                sector.strip_names.join(","),
                sector.inner_patch_name,
                sector.interface_names.join(","),
            )
        })
        .collect::<Vec<_>>()
        .join("|");
    format!(
        "sectors={},strips={},inners={},interfaces={},strip_cps={},inner_cps={},mismatch={},names={}",
        inventory.attempted_sector_count,
        inventory.attempted_strip_patch_count,
        inventory.attempted_inner_patch_count,
        inventory.attempted_interface_count,
        inventory.strip_control_point_count,
        inventory.inner_control_point_count,
        inventory.inventory_mismatch.unwrap_or("none"),
        if sectors.is_empty() { "none" } else { &sectors },
    )
}

fn inner_patch_anchor_label(
    inner: &super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport,
) -> String {
    let mut patches = inner
        .inner_patches
        .iter()
        .map(|patch| {
            let center = patch.surface.cps[3][3];
            format!(
                "{}@{}:[{}]({:.3e},{:.3e},{:.3e})",
                patch.patch_name,
                patch.sector_index,
                patch.boundary_interface_names.join(","),
                center.x,
                center.y,
                center.z
            )
        })
        .collect::<Vec<_>>();
    patches.sort();
    if patches.is_empty() {
        "none".to_string()
    } else {
        patches.join("|")
    }
}

fn interface_anchor_label(
    interfaces: &super::kp20_polar_multipatch_interface::Kp20MultiPatchInterfaceBuildReport,
) -> String {
    let mut anchors = interfaces
        .interface_curves
        .iter()
        .map(|interface| {
            let start = interface.controls[0];
            let end = interface.controls[6];
            format!(
                "{}@{}:{}({:.3e},{:.3e},{:.3e})->({:.3e},{:.3e},{:.3e})",
                interface.interface_name,
                interface.sector_index,
                interface.strip_name,
                start.x,
                start.y,
                start.z,
                end.x,
                end.y,
                end.z
            )
        })
        .collect::<Vec<_>>();
    anchors.sort();
    if anchors.is_empty() {
        "none".to_string()
    } else {
        anchors.join("|")
    }
}

fn strip_surface_sector_label(
    strips: &super::kp20_polar_multipatch_strip::Kp20MultiPatchStripBuildReport,
) -> String {
    let mut sectors = strips
        .strip_surfaces
        .iter()
        .map(|strip| {
            let anchor = strip.surface.cps[0][0];
            format!(
                "{}@{}:{:?}({:.3e},{:.3e},{:.3e})",
                strip.patch_name,
                strip.sector_index,
                strip.strip_piece,
                anchor.x,
                anchor.y,
                anchor.z
            )
        })
        .collect::<Vec<_>>();
    sectors.sort();
    if sectors.is_empty() {
        "none".to_string()
    } else {
        sectors.join("|")
    }
}

fn optional_usize_label(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "none".to_string())
}

fn optional_f64_label(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.3e}"))
        .unwrap_or_else(|| "none".to_string())
}

fn optional_piece_counts_label(value: Option<[usize; 4]>) -> String {
    value
        .map(|counts| {
            format!(
                "NearOriginLeft:{}|NearOriginRight:{}|NearMidpointRight:{}|NearMidpointLeft:{}",
                counts[0], counts[1], counts[2], counts[3]
            )
        })
        .unwrap_or_else(|| "none".to_string())
}
