//! Candidate-only KP20 strip/inner interface curve extraction.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::Kp20PieceId;

use super::kp20_polar_multipatch_strip::Kp20MultiPatchStripBuildReport;

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchInterfaceBuildReport {
    pub(super) attempted_interface_count: usize,
    pub(super) built_interface_count: usize,
    pub(super) finite_control_point_count: usize,
    pub(super) invalid_control_point_count: usize,
    pub(super) decline_reason: Option<&'static str>,
    pub(super) interface_names: Vec<String>,
    pub(super) interface_curves: Vec<Kp20MultiPatchInterfaceCurve>,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchInterfaceCurve {
    pub(super) interface_name: String,
    pub(super) strip_name: String,
    pub(super) sector_index: usize,
    pub(super) strip_piece: Kp20PieceId,
    pub(super) controls: [Point3; 7],
    pub(super) support_adjacent_controls: [Point3; 7],
}

pub(super) fn build_multipatch_interfaces(
    strips: &Kp20MultiPatchStripBuildReport,
) -> Kp20MultiPatchInterfaceBuildReport {
    let mut interface_curves = Vec::new();
    let mut finite_control_point_count = 0_usize;
    let mut invalid_control_point_count = 0_usize;
    for strip in &strips.strip_surfaces {
        let controls = strip.surface.cps[6];
        let support_adjacent_controls = strip.surface.cps[5];
        let (finite, invalid) = count_control_points(&controls, &support_adjacent_controls);
        finite_control_point_count += finite;
        invalid_control_point_count += invalid;
        interface_curves.push(Kp20MultiPatchInterfaceCurve {
            interface_name: interface_name_from_strip(&strip.patch_name),
            strip_name: strip.patch_name.clone(),
            sector_index: strip.sector_index,
            strip_piece: strip.strip_piece,
            controls,
            support_adjacent_controls,
        });
    }
    let interface_names = interface_curves
        .iter()
        .map(|interface| interface.interface_name.clone())
        .collect::<Vec<_>>();
    let decline_reason = if interface_curves.len() != strips.built_strip_patch_count {
        Some("missing_interface_curve")
    } else if invalid_control_point_count > 0 {
        Some("non_finite_interface_control_point")
    } else {
        None
    };
    Kp20MultiPatchInterfaceBuildReport {
        attempted_interface_count: strips.built_strip_patch_count,
        built_interface_count: interface_curves.len(),
        finite_control_point_count,
        invalid_control_point_count,
        decline_reason,
        interface_names,
        interface_curves,
    }
}

fn interface_name_from_strip(strip_name: &str) -> String {
    strip_name.replace("_support_strip", "_support_interface")
}

fn count_control_points(
    controls: &[Point3; 7],
    support_adjacent_controls: &[Point3; 7],
) -> (usize, usize) {
    let mut finite = 0_usize;
    let mut invalid = 0_usize;
    for point in controls.iter().chain(support_adjacent_controls) {
        if point.x.is_finite() && point.y.is_finite() && point.z.is_finite() {
            finite += 1;
        } else {
            invalid += 1;
        }
    }
    (finite, invalid)
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

    use super::*;
    use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};
    use crate::construct::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;
    use crate::construct::kp20_polar_multipatch_strip::build_multipatch_strip_surfaces;

    #[test]
    fn interface_builder_extracts_strip_inner_boundary() {
        let summary = Kp20SamplerFoldSummary {
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
        };
        let inventory = multipatch_inventory_from_fold_summary(&summary);
        let mut macro_patch = Kp20MacroPatch::zero();
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][3] = Point3::new(1.0, 2.0, 3.0);
        let strips = build_multipatch_strip_surfaces(&inventory, &[macro_patch], &summary);

        let interfaces = build_multipatch_interfaces(&strips);

        assert_eq!(interfaces.attempted_interface_count, 1);
        assert_eq!(interfaces.built_interface_count, 1);
        assert_eq!(interfaces.finite_control_point_count, 14);
        assert_eq!(interfaces.invalid_control_point_count, 0);
        assert_eq!(interfaces.decline_reason, None);
        assert_eq!(interfaces.interface_curves[0].sector_index, 0);
        assert_eq!(
            interfaces.interface_curves[0].strip_piece,
            Kp20PieceId::NearMidpointLeft
        );
        assert!(
            interfaces.interface_curves[0]
                .interface_name
                .ends_with("_support_interface")
        );
        assert!(
            interfaces.interface_curves[0]
                .strip_name
                .ends_with("_support_strip")
        );
        assert_eq!(
            interfaces.interface_curves[0].controls[3],
            Point3::new(1.0, 2.0, 3.0)
        );
        assert_eq!(
            interfaces.interface_curves[0].support_adjacent_controls[3],
            Point3::new(0.0, 0.0, 0.0)
        );
    }
}
