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
        "attempted_strips={},attempted_inners={},changed_strips={},changed_inners={},changed_pieces={},fast_evaluated={},exact_evaluated={},fast_folds={},exact_folds={},fast_min_cos={},exact_min_cos={},validation_clean={},accepted={},decline={}",
        evaluation.attempted_strip_patch_count,
        evaluation.attempted_inner_patch_count,
        evaluation.changed_strip_piece_count,
        evaluation.changed_inner_piece_count,
        evaluation.changed_piece_count,
        evaluation.fast_evaluated,
        evaluation.exact_evaluated,
        optional_usize_label(evaluation.fast_fold_count),
        optional_usize_label(evaluation.exact_fold_count),
        optional_f64_label(evaluation.fast_min_jacobian_normal_cos),
        optional_f64_label(evaluation.exact_min_jacobian_normal_cos),
        evaluation.validation_clean,
        evaluation.accepted,
        evaluation.decline_reason.unwrap_or("none"),
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
        "attempted={},built={},reauthored={},finite_cps={},invalid_cps={},max_interior_displacement={:.3e},strip_cps={},surface_count={},surface_sectors={},decline={},names={}",
        strips.attempted_strip_patch_count,
        strips.built_strip_patch_count,
        strips.reauthored_strip_patch_count,
        strips.finite_control_point_count,
        strips.invalid_control_point_count,
        strips.max_reauthored_interior_displacement,
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
