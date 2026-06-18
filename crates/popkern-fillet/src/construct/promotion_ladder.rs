//! Promotion ladder for the three-quad decomposition.
//!
//! Walks the rung hierarchy:
//! deg5 → biquintic-full → biseptic → two-stripe → KP13 bi-5
//! (env-gated) → KP20 bi-6 macro (env-gated) → bi-6 polar →
//! KP binonic.
//!
//! Extracted from `three_quad_apply.rs` to keep that file under
//! the 500-line structural limit.

use popkern_core::geometry::gregory::GregoryPatch;

use super::three_quad_kp_diagnostics;
use super::three_quad_promotions::{
    promote_decomp_to_biquintic_full_with_spoke_g1, promote_decomp_to_biseptic_with_spoke_g1,
    promote_decomp_to_two_stripe_with_spoke_g1,
};
use crate::blend::gregory_three_quad::NQuadDecomposition;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::ApexConstructionDiagnostic;
use crate::types::CornerSubPatchFamily;

/// Result of running the promotion ladder.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(super) enum PromotionResult {
    /// Caller should use the ordinary N-way poke.
    #[default]
    Standard,
    /// Caller should use the 4N macro-cap poke with this construction family.
    Macro4n(CornerSubPatchFamily),
    /// Caller should use the expanded KP20 shared-trim support topology.
    ExpandedKp20Support,
    /// Caller should use the KP20 support-cell topology.
    Kp20SupportCell,
}

#[derive(Debug, Clone, Default)]
pub(crate) struct KpBinonicPromotionReport {
    pub relaxed_coplanarity: Option<f64>,
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

#[derive(Debug, Clone, Default)]
pub(super) struct PromotionLadderReport {
    pub result: PromotionResult,
    pub kp20_support_cell_floor_patches: Option<Vec<GregoryPatch>>,
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

/// Run the full promotion ladder and return its selected topology
/// family plus structured diagnostics for declined promotion rungs.
///
/// Phase 6c-B: `promote_binonic` callback takes both `sub_patches`
/// and `fillet_surfaces` so the binonic optimiser can use outer-
/// fillet G1 residuals (cap-vs-fillet tangent-plane match along
/// each cap's outer boundary).
#[allow(clippy::too_many_arguments)]
pub(super) fn run_promotion_ladder(
    decomp: &mut NQuadDecomposition,
    fillet_surfaces: &[FilletSurfaceResult],
    relaxed_coplanarity: &mut Option<f64>,
    use_deg5: bool,
    use_biquintic_full: bool,
    use_biseptic: bool,
    use_two_stripe: bool,
    use_kp_binonic: bool,
    use_kp_polar: bool,
    promote_binonic: impl FnOnce(
        &mut [GregoryPatch],
        &[FilletSurfaceResult],
    ) -> KpBinonicPromotionReport,
) -> PromotionLadderReport {
    let n_sides = decomp.sub_patches.len();
    if !(use_biquintic_full && use_deg5) {
        return standard_report_with_prerequisite_decline(
            use_kp_binonic,
            use_kp_polar,
            decomp.sub_patches.len(),
            "pre_biquintic_disabled",
            Some("deg5 or biquintic-full promotion disabled".to_string()),
        );
    }
    promote_decomp_to_biquintic_full_with_spoke_g1(&mut decomp.sub_patches);
    if !use_biseptic || !use_two_stripe {
        return standard_report_with_prerequisite_decline(
            use_kp_binonic,
            use_kp_polar,
            decomp.sub_patches.len(),
            "pre_two_stripe_disabled",
            Some("biseptic or two-stripe promotion disabled".to_string()),
        );
    }
    if !promote_decomp_to_biseptic_with_spoke_g1(&mut decomp.sub_patches) {
        return standard_report_with_prerequisite_decline(
            use_kp_binonic,
            use_kp_polar,
            decomp.sub_patches.len(),
            "biseptic_unavailable",
            None,
        );
    }
    if !promote_decomp_to_two_stripe_with_spoke_g1(&mut decomp.sub_patches) {
        return standard_report_with_prerequisite_decline(
            use_kp_binonic,
            use_kp_polar,
            decomp.sub_patches.len(),
            "two_stripe_unavailable",
            None,
        );
    }
    if use_kp_polar
        && decomp.sub_patches.len() >= 3
        && let Some(patches_4n) =
            super::kp13_bi5_promote::try_kp13_bi5_dispatch(&decomp.sub_patches, fillet_surfaces)
    {
        decomp.sub_patches = patches_4n;
        return PromotionLadderReport {
            result: PromotionResult::Macro4n(CornerSubPatchFamily::Kp13Bi5),
            kp20_support_cell_floor_patches: None,
            diagnostics: Vec::new(),
        };
    }
    let mut diagnostics = Vec::new();
    if use_kp_polar && decomp.sub_patches.len() >= 3 {
        let mut report = super::kp20_polar_dispatch::try_kp20_polar_dispatch_with_diagnostics(
            &decomp.sub_patches,
            fillet_surfaces,
        );
        diagnostics.append(&mut report.diagnostics);
        if let Some(patches_4n) = report.patches {
            decomp.sub_patches = patches_4n;
            let result = if report.expanded_support_topology_accepted {
                PromotionResult::ExpandedKp20Support
            } else if report.support_cell_topology_accepted {
                PromotionResult::Kp20SupportCell
            } else {
                PromotionResult::Macro4n(CornerSubPatchFamily::Kp20Bi6)
            };
            return PromotionLadderReport {
                result,
                kp20_support_cell_floor_patches: report.support_cell_floor_patches,
                diagnostics,
            };
        }
    }
    let (polar_committed, mut polar_diagnostics) = if use_kp_polar && decomp.sub_patches.len() >= 5
    {
        super::kp_polar_bi6_promote::promote_and_update_relaxed_with_diagnostics(
            &mut decomp.sub_patches,
            fillet_surfaces,
            relaxed_coplanarity,
        )
    } else {
        (false, Vec::new())
    };
    diagnostics.append(&mut polar_diagnostics);
    if use_kp_binonic && !polar_committed {
        let report = promote_binonic(&mut decomp.sub_patches, fillet_surfaces);
        *relaxed_coplanarity = report.relaxed_coplanarity;
        diagnostics.extend(report.diagnostics);
    }
    if use_kp_polar && n_sides >= 5 {
        diagnostics.push(
            three_quad_kp_diagnostics::support_cell_candidate_nonselection(
                "promotion_ladder_standard",
                n_sides,
                decomp.sub_patches.len(),
            ),
        );
    }
    PromotionLadderReport {
        result: PromotionResult::Standard,
        kp20_support_cell_floor_patches: None,
        diagnostics,
    }
}

fn standard_report_with_prerequisite_decline(
    use_kp_binonic: bool,
    use_kp_polar: bool,
    n_sides: usize,
    reason_code: &'static str,
    detail: Option<String>,
) -> PromotionLadderReport {
    let mut diagnostics = Vec::new();
    if use_kp_binonic {
        diagnostics.push(three_quad_kp_diagnostics::promotion_declined(
            reason_code,
            n_sides,
            detail,
            None,
            None,
        ));
    }
    if use_kp_polar && n_sides >= 5 {
        diagnostics.push(
            three_quad_kp_diagnostics::support_cell_candidate_nonselection(
                reason_code,
                n_sides,
                n_sides,
            ),
        );
    }
    PromotionLadderReport {
        result: PromotionResult::Standard,
        kp20_support_cell_floor_patches: None,
        diagnostics,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    #[test]
    fn binonic_prerequisite_decline_is_structured_diagnostic() {
        let mut decomp = NQuadDecomposition {
            central_point: Point3::new(0.0, 0.0, 0.0),
            sub_patches: Vec::new(),
        };
        let mut relaxed_coplanarity = None;
        let report = run_promotion_ladder(
            &mut decomp,
            &[],
            &mut relaxed_coplanarity,
            false,
            false,
            true,
            true,
            true,
            false,
            |_, _| KpBinonicPromotionReport::default(),
        );
        assert_eq!(report.result, PromotionResult::Standard);
        assert_eq!(report.diagnostics.len(), 1);
        let crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(cap) =
            &report.diagnostics[0]
        else {
            panic!("expected exact-boundary cap diagnostic");
        };
        assert_eq!(cap.code, "kp_binonic.promotion.declined");
        assert_eq!(cap.stage, "kp_binonic.promotion.pre_biquintic_disabled");
        assert_eq!(cap.sector_name, "kp_binonic_ring_n0");
    }
}
