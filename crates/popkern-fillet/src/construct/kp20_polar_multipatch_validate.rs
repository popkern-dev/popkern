//! Candidate-only KP20 multipatch construction validation.

use popkern_core::geometry::Point3;

use super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerBuildReport;
use super::kp20_polar_multipatch_interface::{
    Kp20MultiPatchInterfaceBuildReport, Kp20MultiPatchInterfaceCurve,
};
use super::kp20_polar_multipatch_inventory::Kp20MultiPatchInventoryReport;
use super::kp20_polar_multipatch_strip::{
    Kp20MultiPatchStripBuildReport, Kp20MultiPatchStripSurface,
};

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchValidationReport {
    pub(super) validated_sector_count: usize,
    pub(super) validated_strip_count: usize,
    pub(super) validated_interface_count: usize,
    pub(super) validated_inner_patch_count: usize,
    pub(super) finite_control_point_count: usize,
    pub(super) invalid_control_point_count: usize,
    pub(super) validated_boundary_constraint_count: usize,
    pub(super) validated_adjacent_constraint_count: usize,
    pub(super) exact_boundary_owned_endpoint_count: usize,
    pub(super) max_strip_interface_g0_residual: f64,
    pub(super) max_strip_adjacent_residual: f64,
    pub(super) max_inner_boundary_g0_residual: f64,
    pub(super) max_inner_adjacent_residual: f64,
    pub(super) decline_reason: Option<&'static str>,
}

pub(super) fn validate_multipatch_candidate(
    inventory: &Kp20MultiPatchInventoryReport,
    strips: &Kp20MultiPatchStripBuildReport,
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
    inner: &Kp20MultiPatchInnerBuildReport,
) -> Kp20MultiPatchValidationReport {
    let mut report = Kp20MultiPatchValidationReport {
        validated_sector_count: 0,
        validated_strip_count: 0,
        validated_interface_count: 0,
        validated_inner_patch_count: 0,
        finite_control_point_count: 0,
        invalid_control_point_count: 0,
        validated_boundary_constraint_count: 0,
        validated_adjacent_constraint_count: 0,
        exact_boundary_owned_endpoint_count: 0,
        max_strip_interface_g0_residual: 0.0,
        max_strip_adjacent_residual: 0.0,
        max_inner_boundary_g0_residual: 0.0,
        max_inner_adjacent_residual: 0.0,
        decline_reason: None,
    };

    for sector in &inventory.sector_reports {
        let Some(inner_patch) = inner
            .inner_patches
            .iter()
            .find(|patch| patch.patch_name == sector.inner_patch_name)
        else {
            report.decline_reason = report
                .decline_reason
                .or(Some("missing_validated_inner_patch"));
            continue;
        };
        report.validated_sector_count += 1;
        report.validated_inner_patch_count += 1;
        count_surface_controls(&inner_patch.surface, &mut report);

        for (boundary_index, interface_name) in sector.interface_names.iter().enumerate() {
            let Some(interface) = interface_by_name(interfaces, interface_name) else {
                report.decline_reason = report
                    .decline_reason
                    .or(Some("missing_validated_interface"));
                continue;
            };
            let Some(strip) = strip_by_name(strips, &interface.strip_name) else {
                report.decline_reason = report.decline_reason.or(Some("missing_validated_strip"));
                continue;
            };
            report.validated_interface_count += 1;
            report.validated_strip_count += 1;
            count_strip_once(strips, strip, &mut report);
            report.max_strip_interface_g0_residual = report
                .max_strip_interface_g0_residual
                .max(max_row_residual(&strip.surface.cps[6], &interface.controls));
            report.max_strip_adjacent_residual = report.max_strip_adjacent_residual.max(
                max_row_residual(&strip.surface.cps[5], &interface.support_adjacent_controls),
            );
            report.max_inner_boundary_g0_residual =
                report
                    .max_inner_boundary_g0_residual
                    .max(inner_boundary_residual(
                        inner_patch,
                        boundary_index,
                        interface,
                    ));
            report.max_inner_adjacent_residual =
                report
                    .max_inner_adjacent_residual
                    .max(inner_adjacent_residual(
                        inner_patch,
                        boundary_index,
                        interface,
                        sector.interface_names.len(),
                    ));
            report.validated_boundary_constraint_count += interface.controls.len();
            report.validated_adjacent_constraint_count += ADJACENT_CONTROL_COUNT;
            report.exact_boundary_owned_endpoint_count += exact_boundary_owned_endpoint_count_for(
                boundary_index,
                sector.interface_names.len(),
            );
        }
    }
    report.decline_reason = report
        .decline_reason
        .or_else(|| validation_decline(&report));
    report
}

fn validation_decline(report: &Kp20MultiPatchValidationReport) -> Option<&'static str> {
    if report.invalid_control_point_count > 0 {
        Some("non_finite_multipatch_control_point")
    } else if report.max_strip_interface_g0_residual > EXACT_BOUNDARY_G0_TOLERANCE {
        Some("strip_interface_g0_residual")
    } else if report.max_strip_adjacent_residual > EXACT_BOUNDARY_G0_TOLERANCE {
        Some("strip_adjacent_residual")
    } else if report.max_inner_boundary_g0_residual > EXACT_BOUNDARY_G0_TOLERANCE {
        Some("inner_interface_g0_residual")
    } else if report.max_inner_adjacent_residual > EXACT_BOUNDARY_G0_TOLERANCE {
        Some("inner_adjacent_residual")
    } else {
        None
    }
}

fn interface_by_name<'a>(
    interfaces: &'a Kp20MultiPatchInterfaceBuildReport,
    name: &str,
) -> Option<&'a Kp20MultiPatchInterfaceCurve> {
    interfaces
        .interface_curves
        .iter()
        .find(|interface| interface.interface_name == name)
}

fn strip_by_name<'a>(
    strips: &'a Kp20MultiPatchStripBuildReport,
    name: &str,
) -> Option<&'a Kp20MultiPatchStripSurface> {
    strips
        .strip_surfaces
        .iter()
        .find(|strip| strip.patch_name == name)
}

fn count_strip_once(
    strips: &Kp20MultiPatchStripBuildReport,
    strip: &Kp20MultiPatchStripSurface,
    report: &mut Kp20MultiPatchValidationReport,
) {
    if strips
        .strip_surfaces
        .iter()
        .take_while(|candidate| candidate.patch_name != strip.patch_name)
        .any(|candidate| candidate.patch_name == strip.patch_name)
    {
        return;
    }
    count_surface_controls(&strip.surface, report);
}

fn count_surface_controls(
    surface: &popkern_core::geometry::gregory::Kp20Bi6Patch,
    report: &mut Kp20MultiPatchValidationReport,
) {
    for row in surface.cps {
        for point in row {
            if point.x.is_finite() && point.y.is_finite() && point.z.is_finite() {
                report.finite_control_point_count += 1;
            } else {
                report.invalid_control_point_count += 1;
            }
        }
    }
}

fn inner_boundary_residual(
    patch: &super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerPatch,
    boundary_index: usize,
    interface: &Kp20MultiPatchInterfaceCurve,
) -> f64 {
    let mut residual = 0.0_f64;
    for idx in 0..7 {
        residual = residual.max(
            inner_boundary_control(patch, boundary_index, idx)
                .distance_to(&interface.controls[idx]),
        );
    }
    residual
}

fn inner_adjacent_residual(
    patch: &super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerPatch,
    boundary_index: usize,
    interface: &Kp20MultiPatchInterfaceCurve,
    boundary_count: usize,
) -> f64 {
    let mut residual = 0.0_f64;
    for idx in 0..7 {
        if endpoint_exact_boundary_owner(boundary_index, idx)
            .is_some_and(|owner| owner < boundary_count)
        {
            continue;
        }
        let expected = support_tangent_adjacent_control(
            interface.controls[idx],
            interface.support_adjacent_controls[idx],
        );
        residual =
            residual.max(inner_adjacent_control(patch, boundary_index, idx).distance_to(&expected));
    }
    residual
}

fn inner_boundary_control(
    patch: &super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerPatch,
    boundary_index: usize,
    idx: usize,
) -> Point3 {
    match boundary_index {
        0 => patch.surface.cps[0][idx],
        1 => patch.surface.cps[6][idx],
        2 => patch.surface.cps[idx][0],
        _ => patch.surface.cps[idx][6],
    }
}

fn inner_adjacent_control(
    patch: &super::kp20_polar_multipatch_inner::Kp20MultiPatchInnerPatch,
    boundary_index: usize,
    idx: usize,
) -> Point3 {
    match boundary_index {
        0 => patch.surface.cps[1][idx],
        1 => patch.surface.cps[5][idx],
        2 => patch.surface.cps[idx][1],
        _ => patch.surface.cps[idx][5],
    }
}

fn exact_boundary_owned_endpoint_count_for(boundary_index: usize, boundary_count: usize) -> usize {
    [0_usize, 6]
        .into_iter()
        .filter(|idx| {
            endpoint_exact_boundary_owner(boundary_index, *idx)
                .is_some_and(|owner| owner < boundary_count)
        })
        .count()
}

fn endpoint_exact_boundary_owner(boundary_index: usize, idx: usize) -> Option<usize> {
    match (boundary_index, idx) {
        (0, 0) | (1, 0) => Some(2),
        (0, 6) | (1, 6) => Some(3),
        (2, 0) | (3, 0) => Some(0),
        (2, 6) | (3, 6) => Some(1),
        _ => None,
    }
}

fn support_tangent_adjacent_control(boundary: Point3, support_adjacent: Point3) -> Point3 {
    boundary + (boundary - support_adjacent)
}

fn max_row_residual(a: &[Point3; 7], b: &[Point3; 7]) -> f64 {
    let mut residual = 0.0_f64;
    for idx in 0..7 {
        residual = residual.max(a[idx].distance_to(&b[idx]));
    }
    residual
}

const EXACT_BOUNDARY_G0_TOLERANCE: f64 = 1.0e-12;
const ADJACENT_CONTROL_COUNT: usize = 7;

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

    use super::*;
    use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};
    use crate::construct::kp20_polar_multipatch_inner::build_multipatch_inner_patches;
    use crate::construct::kp20_polar_multipatch_interface::build_multipatch_interfaces;
    use crate::construct::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;
    use crate::construct::kp20_polar_multipatch_strip::build_multipatch_strip_surfaces;

    #[test]
    fn validator_accepts_candidate_object_invariants() {
        let inventory = multipatch_inventory_from_fold_summary(&Kp20SamplerFoldSummary {
            fold_count: 1,
            min_jacobian_normal_cos: -0.1,
            piece_fold_counts: [1, 0, 0, 0],
            fold_details: Vec::new(),
            fold_sites: vec![Kp20FoldSite {
                sector: 0,
                piece: Kp20PieceId::NearOriginLeft,
                u: 0.5,
                v: 0.5,
                min_jacobian_ratio: -0.1,
                min_jacobian_normal_cos: -0.1,
            }],
        });
        let mut macro_patch = Kp20MacroPatch::zero();
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][3] = Point3::new(1.0, 2.0, 3.0);
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[5][3] =
            Point3::new(0.25, 0.5, 0.75);
        let strips = build_multipatch_strip_surfaces(&inventory, &[macro_patch]);
        let interfaces = build_multipatch_interfaces(&strips);
        let inner = build_multipatch_inner_patches(&inventory, &interfaces);

        let report = validate_multipatch_candidate(&inventory, &strips, &interfaces, &inner);

        assert_eq!(report.validated_sector_count, 1);
        assert_eq!(report.validated_strip_count, 1);
        assert_eq!(report.validated_interface_count, 1);
        assert_eq!(report.validated_inner_patch_count, 1);
        assert_eq!(report.invalid_control_point_count, 0);
        assert_eq!(report.validated_boundary_constraint_count, 7);
        assert_eq!(report.validated_adjacent_constraint_count, 7);
        assert_eq!(report.exact_boundary_owned_endpoint_count, 0);
        assert_eq!(report.max_strip_interface_g0_residual, 0.0);
        assert_eq!(report.max_strip_adjacent_residual, 0.0);
        assert_eq!(report.max_inner_boundary_g0_residual, 0.0);
        assert_eq!(report.max_inner_adjacent_residual, 0.0);
        assert_eq!(report.decline_reason, None);
    }
}
