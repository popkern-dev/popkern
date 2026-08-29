//! Fold-triggered balanced-fan retry for the KP binonic rung
//! (PK-240).
//!
//! The straight-spoke fan is always tried first — every serving
//! class keeps its exact behavior. When the binonic promotion
//! declines with a self-intersection (the guided construction folds
//! against a degenerate sector split — measured 173.6° at 1:20 vs
//! the regular 120° chart its tables assume), the decomposition is
//! rebuilt with balanced spoke tangents and the ladder re-runs. The
//! retry walks the balancing strength weakest-first because full
//! 2π/N regularization can over-distort the smallest sector (its
//! outer arc is tiny while its angular span grows to 2π/N; the
//! initial fan itself folds pre-solve — measured at 1:20/1:30),
//! while the mildest balance that unfolds the guided construction
//! serves best. A retry is accepted ONLY when every sub-patch
//! committed binonic, the ring boundary held, and the candidate is
//! fold-free with margin; otherwise the straight-spoke result
//! stands and the recorded fold decline remains the serving truth.

use popkern_core::geometry::gregory::GregoryPatch;

use super::promotion_ladder::{KpBinonicPromotionReport, PromotionLadderReport};
use crate::blend::boundary::BlendBoundary;
use crate::blend::gregory_three_quad::{NQuadDecomposition, build_n_quad_decomposition_balanced};
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::ApexConstructionDiagnostic;
use crate::vertex_classifier::VertexBlendClassification;

/// Balancing strengths tried in order. Weakest first — measured on
/// the 1:20–1:100 ladder, the mildest balance that unfolds the
/// guided construction serves best (least geometric distortion of
/// the fan: 0.33 at 1:50 measured ring G1 2.7e-3 while 0.66 at 1:20
/// measured 1.7e-2), and stronger balances increasingly risk the
/// initial-fan fold screen and the outer-boundary residual gate.
const RETRY_STRENGTHS: [f64; 5] = [0.33, 0.5, 0.66, 0.85, 1.0];

/// Fold margin for retry acceptance: the ladder's own gate and
/// `validate_all` both test at 1e-6 arc length, but post-assignment
/// support solves nudge the committed boundary rows (measured
/// 5.6e-8 at 1:20), which flipped a margin-less 0.66 candidate into
/// a served `TwistedSurface`. Accept only candidates fold-free at
/// 10x the validation sensitivity so serving-time nudges cannot
/// cross the line.
const RETRY_FOLD_MARGIN_TOLERANCE: f64 = 1.0e-7;

/// Ladder flags forwarded unchanged to the retry runs.
#[derive(Clone, Copy)]
pub(super) struct BalancedRetryLadderFlags {
    pub use_deg5: bool,
    pub use_biquintic_full: bool,
    pub use_biseptic: bool,
    pub use_two_stripe: bool,
    pub use_kp_binonic: bool,
    pub use_kp_polar: bool,
}

/// Accepted balanced retry: the decomposition and ladder report to
/// serve in place of the straight-spoke attempt.
pub(super) struct BalancedRetryOutcome {
    pub decomp: NQuadDecomposition,
    pub report: PromotionLadderReport,
    pub relaxed_coplanarity: Option<f64>,
}

/// What the retry ladder did: an accepted outcome, or the
/// per-strength decline reasons (empty when the natural split never
/// engaged balancing at all).
pub(super) enum BalancedRetryResult {
    Accepted(Box<BalancedRetryOutcome>),
    Declined(Vec<String>),
}

/// Driver: when the straight-spoke ladder fold-declined the binonic,
/// run the retry and either serve its outcome (mutating the caller's
/// decomposition/report/relaxed state, preserving the straight
/// attempt's diagnostics as provenance) or record the per-strength
/// decline reasons.
#[allow(clippy::too_many_arguments)]
pub(super) fn apply_fold_triggered_balanced_retry(
    decomp: &mut NQuadDecomposition,
    promotion_report: &mut PromotionLadderReport,
    relaxed_coplanarity: &mut Option<f64>,
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    boundary: &BlendBoundary,
    brep: &popkern_core::topology::BRep,
    fillet_surfaces: &[FilletSurfaceResult],
    classification: VertexBlendClassification,
    flags: BalancedRetryLadderFlags,
    promote_binonic: impl Fn(&mut [GregoryPatch], &[FilletSurfaceResult]) -> KpBinonicPromotionReport,
) {
    if !flags.use_kp_binonic || !binonic_fold_declined(&promotion_report.diagnostics) {
        return;
    }
    match retry_binonic_with_balanced_fan(
        boundary,
        brep,
        fillet_surfaces,
        classification,
        flags,
        promote_binonic,
    ) {
        BalancedRetryResult::Accepted(outcome) => {
            // Keep the straight attempt's diagnostics (the fold
            // decline is real provenance), then serve the balanced
            // result.
            diagnostics.append(&mut promotion_report.diagnostics);
            *decomp = outcome.decomp;
            *promotion_report = outcome.report;
            *relaxed_coplanarity = outcome.relaxed_coplanarity;
        }
        BalancedRetryResult::Declined(reasons) if !reasons.is_empty() => {
            diagnostics.push(
                super::three_quad_kp_diagnostics::balanced_spoke_fan_declined(&reasons.join(";")),
            );
        }
        BalancedRetryResult::Declined(_) => {}
    }
}

/// The `kp_binonic.promotion.boundary_candidates` fact's detail
/// string from a ladder report, when the promotion recorded one.
fn boundary_candidates_detail(diagnostics: &[ApexConstructionDiagnostic]) -> Option<String> {
    diagnostics.iter().find_map(|diagnostic| {
        let ApexConstructionDiagnostic::ExactBoundaryCap(d) = diagnostic else {
            return None;
        };
        (d.stage == "kp_binonic.promotion.boundary_candidates")
            .then(|| d.detail_code.clone())
            .flatten()
    })
}

/// Whether the ladder report carries the binonic self-intersection
/// decline that arms the retry.
fn binonic_fold_declined(diagnostics: &[ApexConstructionDiagnostic]) -> bool {
    diagnostics.iter().any(|diagnostic| {
        let ApexConstructionDiagnostic::ExactBoundaryCap(d) = diagnostic else {
            return false;
        };
        d.code == "kp_binonic.promotion.declined"
            && d.stage == "kp_binonic.promotion.self_intersection"
    })
}

/// Run the balanced-fan retry ladder. Returns `Declined` (with
/// per-strength reasons) when no strength both engages and serves —
/// the caller keeps the straight-spoke result and records the
/// reasons.
#[allow(clippy::too_many_arguments)]
pub(super) fn retry_binonic_with_balanced_fan(
    boundary: &BlendBoundary,
    brep: &popkern_core::topology::BRep,
    fillet_surfaces: &[FilletSurfaceResult],
    classification: VertexBlendClassification,
    flags: BalancedRetryLadderFlags,
    promote_binonic: impl Fn(&mut [GregoryPatch], &[FilletSurfaceResult]) -> KpBinonicPromotionReport,
) -> BalancedRetryResult {
    let mut declined: Vec<String> = Vec::new();
    let mut banded_acceptance: Option<Box<BalancedRetryOutcome>> = None;
    for strength in RETRY_STRENGTHS {
        let Some(mut retry_decomp) = build_n_quad_decomposition_balanced(
            boundary,
            fillet_surfaces,
            flags.use_deg5,
            classification,
            strength,
        ) else {
            declined.push(format!("s{strength:.2}=decomposition_rejected"));
            continue;
        };
        let Some((natural_span, balanced_span)) = retry_decomp.balanced_spokes else {
            // The natural split is not degenerate enough to engage —
            // no strength will change that; the fold has another
            // cause and the recorded decline stands.
            return BalancedRetryResult::Declined(Vec::new());
        };
        let mut relaxed_coplanarity = None;
        let mut report = super::promotion_ladder::run_promotion_ladder(
            &mut retry_decomp,
            fillet_surfaces,
            Some((boundary, brep)),
            &mut relaxed_coplanarity,
            flags.use_deg5,
            flags.use_biquintic_full,
            flags.use_biseptic,
            flags.use_two_stripe,
            flags.use_kp_binonic,
            flags.use_kp_polar,
            &promote_binonic,
        );
        let committed = !retry_decomp.sub_patches.is_empty()
            && retry_decomp
                .sub_patches
                .iter()
                .all(|p| matches!(p, GregoryPatch::RectangularBiNonic(_)));
        // The balanced binonic must also HOLD the ring boundary: at
        // 1:50 the strength-0.66 candidate committed in-ladder with
        // an outer-boundary residual of 1.3e-1, and the downstream
        // residual gate then reverted the corner to the PRE-ladder
        // fan (worse than the straight serve). Gate here so a
        // boundary-breaking retry falls back to the straight
        // POST-ladder serve instead.
        const RETRY_BOUNDARY_RESIDUAL_GATE: f64 = 5.0e-3;
        let measured_residual = committed
            .then(|| {
                super::three_quad_apply::three_quad_apply_helpers::max_outer_boundary_residual(
                    &retry_decomp,
                    boundary,
                )
            })
            .flatten();
        let holds_boundary =
            committed && measured_residual.is_none_or(|r| r <= RETRY_BOUNDARY_RESIDUAL_GATE);
        let fold_margined = committed
            && !retry_decomp
                .sub_patches
                .iter()
                .any(|p| p.self_intersects_with_tolerance(RETRY_FOLD_MARGIN_TOLERANCE));
        if committed && holds_boundary && fold_margined {
            // Chart-health preference (PK-249): the accepted
            // candidate's charts can carry handedness-reversal bands
            // hugging the distorted sector's spoke sides (the
            // Jacobian crosses zero in bands the mid-face twist
            // screen and the arc-length fold gate cannot see; the
            // geometry stays curvature-continuous). A chart-honest
            // candidate (bands 0, possibly via the promotion's
            // orientation repair) is served immediately; a banded
            // one is stashed and the ladder keeps walking — stronger
            // balancing measurably yields chart-clean candidates
            // (1:100 serves clean at 0.85 while 0.33–0.66 band).
            // When no strength reaches a clean chart, the weakest
            // banded acceptance serves exactly as before — never
            // trade a 1e-3 G1 serve for the 1.9e-1 creased fallback.
            let reversal_bands =
                super::three_quad_apply::three_quad_chart_orientation::chart_reversal_band_count(
                    &retry_decomp.sub_patches,
                );
            report.diagnostics.push(
                super::three_quad_kp_diagnostics::balanced_spoke_fan_engaged(
                    natural_span,
                    balanced_span,
                    strength,
                    reversal_bands,
                ),
            );
            let outcome = Box::new(BalancedRetryOutcome {
                decomp: retry_decomp,
                report,
                relaxed_coplanarity,
            });
            if reversal_bands == 0 {
                return BalancedRetryResult::Accepted(outcome);
            }
            declined.push(format!(
                "s{strength:.2}=accepted_with_bands({reversal_bands})"
            ));
            if banded_acceptance.is_none() {
                banded_acceptance = Some(outcome);
            }
            continue;
        }
        // Carry the promotion's candidate-selection anatomy (which
        // boundary-fidelity path served and why the stronger ones
        // fell short) into the per-strength decline reason — the
        // declined runs' reports are dropped, and without this the
        // extreme-asymmetry class (PK-249) records only the final
        // residual, not which candidate produced it.
        let selection = boundary_candidates_detail(&report.diagnostics)
            .map(|detail| format!("[{detail}]"))
            .unwrap_or_default();
        declined.push(if !committed {
            format!("s{strength:.2}=binonic_not_committed{selection}")
        } else if !holds_boundary {
            format!(
                "s{strength:.2}=boundary_residual_exceeded({:.3e}){selection}",
                measured_residual.unwrap_or(f64::NAN)
            )
        } else {
            format!("s{strength:.2}=fold_margin_exceeded{selection}")
        });
    }
    if let Some(mut outcome) = banded_acceptance {
        // The full strength walk stays attributable on the served
        // report: which stronger candidates were tried toward a
        // chart-clean serve and why each fell short.
        outcome
            .report
            .diagnostics
            .push(super::three_quad_kp_diagnostics::balanced_spoke_fan_walk(
                &declined.join(";"),
            ));
        return BalancedRetryResult::Accepted(outcome);
    }
    BalancedRetryResult::Declined(declined)
}
