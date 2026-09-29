//! Evidence collection for the replay report: what one feature
//! execution reports, and the classification of kernel warnings into
//! the report's three record kinds — warning, fallback, skip.
//!
//! The kernel subsystems each report non-fatal findings as one flat
//! warning enum. The report channel separates the ones that mean "a
//! construction served through its explicit weaker path" (fallbacks)
//! and "work was declined while the rest served" (skips) from
//! informational notes, because a downgraded route is a first-class
//! fact a GUI shows next to the geometry, not a line in a log. The
//! mapping is by variant, here, so the kernel enums stay the single
//! source of truth for what happened and this file is the single
//! source of truth for what it means to a document.

use std::collections::BTreeSet;

use popkern_blend::{BlendResult, BlendWarning, ContinuityLevel};
use popkern_bool::merge_fallback::MergeFallbackEvent;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};
use popkern_draft::DraftWarning;
use popkern_shell::ShellWarning;
use popkern_sweep::SweepWarning;

use super::{
    BlendEvidenceRecord, EVIDENCE_SAMPLES_PER_SEAM, FallbackRecord, KernelWarning, RecordSource,
    SkipRecord,
};
use crate::pattern::{InstanceOutcome, PatternResult};

/// What one feature execution collected, before the replay loop
/// attaches identity, outcome, and timing.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct FeatureEvidence {
    /// Non-fatal warnings that are neither fallbacks nor skips.
    pub warnings: Vec<KernelWarning>,
    /// Blend evidence, for fillet and chamfer executions.
    pub blend: Option<BlendEvidenceRecord>,
    /// Fallback firings.
    pub fallbacks: Vec<FallbackRecord>,
    /// Skip records.
    pub skipped: Vec<SkipRecord>,
}

impl FeatureEvidence {
    /// Nothing collected.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Evidence for a served blend: the kernel's per-edge warnings
    /// classified, and — for a fillet — the continuity measured over
    /// the faces the blend created (`pre_faces` is the face census
    /// before the blend). `certification` says whether gated seams are
    /// certified now or deferred (see [`crate::deferred_certification`]).
    pub fn blend(
        result: &BlendResult,
        pre_faces: &BTreeSet<FaceHandle>,
        edges_requested: usize,
        requested: Option<ContinuityLevel>,
        certification: popkern_blend::ContinuityCertification,
    ) -> Self {
        let mut evidence = Self::empty();
        for (edge, warning) in &result.edge_warnings {
            classify_blend_warning(&mut evidence, *edge, warning);
        }
        let continuity = requested.map(|level| {
            let new_faces: Vec<FaceHandle> = result
                .brep
                .faces
                .keys()
                .filter(|f| !pre_faces.contains(f))
                .collect();
            popkern_blend::measure_blend_continuity_evidence_with(
                &result.brep,
                &new_faces,
                level,
                &popkern_blend::ContinuityGates::commercial(),
                EVIDENCE_SAMPLES_PER_SEAM,
                certification,
            )
        });
        evidence.blend = Some(BlendEvidenceRecord {
            edges_requested,
            corner_patches: result.corner_patches.len(),
            requested_continuity: requested.unwrap_or(ContinuityLevel::G0),
            continuity,
            corners: super::replay_report_corners::corner_evidence(result),
        });
        evidence
    }

    /// Evidence for a served shell/thicken/offset.
    pub fn shell(warnings: &[popkern_shell::ShellWarning]) -> Self {
        let mut evidence = Self::empty();
        for warning in warnings {
            classify_shell_warning(&mut evidence, warning);
        }
        evidence
    }

    /// Evidence for a served draft.
    pub fn draft(warnings: &[popkern_draft::DraftWarning]) -> Self {
        let mut evidence = Self::empty();
        for warning in warnings {
            classify_draft_warning(&mut evidence, warning);
        }
        evidence
    }

    /// Evidence for a served sweep/loft/pipe/revolve.
    pub fn sweep(warnings: &[popkern_sweep::SweepWarning]) -> Self {
        let mut evidence = Self::empty();
        for warning in warnings {
            classify_sweep_warning(&mut evidence, warning);
        }
        evidence
    }

    /// Record the boolean pipeline's cellular→graph-merge fallback
    /// firings drained from [`popkern_bool::merge_fallback`].
    pub fn record_boolean_fallbacks(&mut self, events: Vec<MergeFallbackEvent>) {
        for event in events {
            self.fallbacks.push(FallbackRecord {
                source: RecordSource::Boolean,
                mechanism: format!(
                    "{:?}: cellular selection -> graph-merge selection",
                    event.operation
                ),
                reason: format!("{:?}", event.reason),
            });
        }
    }

    /// Evidence for a served pattern: transform warnings, per-instance
    /// geometry fallbacks, suppressed instances, and instances that
    /// could not be produced.
    pub fn pattern(result: &PatternResult) -> Self {
        let mut evidence = Self::empty();
        let (warnings, outcomes, failed): (&[String], Vec<(usize, &InstanceOutcome)>, &[_]) =
            match result {
                PatternResult::AllSucceeded { warnings, .. } => (warnings, Vec::new(), &[]),
                PatternResult::PartialSuccess {
                    succeeded,
                    failed,
                    warnings,
                } => (
                    warnings,
                    succeeded.iter().map(|(i, o)| (*i, o)).collect(),
                    failed,
                ),
                PatternResult::AllFailed { errors } => (&[], Vec::new(), errors),
            };
        for message in warnings {
            evidence.warnings.push(KernelWarning {
                source: RecordSource::Pattern,
                message: message.clone(),
            });
        }
        for (index, outcome) in outcomes {
            classify_instance_outcome(&mut evidence, index, outcome);
        }
        for (index, error) in failed {
            evidence.skipped.push(SkipRecord {
                source: RecordSource::Pattern,
                what: format!("instance {index}"),
                reason: format!("could not be produced: {error}"),
            });
        }
        evidence
    }
}

/// Face census of `brep` before an operation, for blend evidence.
pub(crate) fn face_census(brep: &BRep) -> BTreeSet<FaceHandle> {
    brep.faces.keys().collect()
}

/// Classify one per-edge blend warning.
pub(super) fn classify_blend_warning(
    evidence: &mut FeatureEvidence,
    edge: EdgeHandle,
    warning: &BlendWarning,
) {
    let source = RecordSource::Blend;
    match warning {
        BlendWarning::CollapsedEdgeSkipped { .. } => evidence.skipped.push(SkipRecord {
            source,
            what: format!("{edge:?}"),
            reason: format!("{warning:?}"),
        }),
        BlendWarning::FilletPairDropped {
            fillet_a,
            fillet_b,
            reason,
        } => evidence.skipped.push(SkipRecord {
            source,
            what: format!("fillet specs {fillet_a} and {fillet_b}"),
            reason: format!("strip interference could not be completed: {reason}"),
        }),
        BlendWarning::VertexModeFallback { .. }
        | BlendWarning::ContinuityDowngraded { .. }
        | BlendWarning::CornerContinuityDowngraded { .. }
        | BlendWarning::CornerContinuityRelaxed { .. }
        | BlendWarning::SurfaceBridgeFallback { .. }
        | BlendWarning::CornerRemediationOverlap { .. }
        | BlendWarning::PolarKpBinonicToBicubic { .. }
        | BlendWarning::PolarKpToSphereFit { .. }
        | BlendWarning::PolarKpToPolarNurbs { .. }
        | BlendWarning::PolarKpToMitredPlane { .. }
        | BlendWarning::DirectApexFellBackToVertexBlend { .. }
        | BlendWarning::CapInsertionDeclined { .. } => evidence.fallbacks.push(FallbackRecord {
            source,
            mechanism: format!("{edge:?}: {}", variant_name(warning)),
            reason: format!("{warning:?}"),
        }),
        _ => evidence.warnings.push(KernelWarning {
            source,
            message: format!("{edge:?}: {warning:?}"),
        }),
    }
}

/// Classify one shell/thicken/offset warning.
pub(super) fn classify_shell_warning(evidence: &mut FeatureEvidence, warning: &ShellWarning) {
    let source = RecordSource::Shell;
    match warning {
        ShellWarning::FaceAutoRemoved { face, reason } => evidence.skipped.push(SkipRecord {
            source,
            what: format!("{face:?}"),
            reason: reason.clone(),
        }),
        ShellWarning::JoinStyleDowngraded { .. }
        | ShellWarning::JoinStyleDowngradedAtEdge { .. } => {
            evidence.fallbacks.push(FallbackRecord {
                source,
                mechanism: "join style Arc -> Intersection".to_string(),
                reason: warning.to_string(),
            })
        }
        ShellWarning::FaceOffsetApproximated { .. } => evidence.fallbacks.push(FallbackRecord {
            source,
            mechanism: "exact offset surface -> fitted approximation".to_string(),
            reason: warning.to_string(),
        }),
        _ => evidence.warnings.push(KernelWarning {
            source,
            message: warning.to_string(),
        }),
    }
}

/// Classify one draft warning.
pub(super) fn classify_draft_warning(evidence: &mut FeatureEvidence, warning: &DraftWarning) {
    let source = RecordSource::Draft;
    match warning {
        DraftWarning::AdjacentFallbackToTransition { .. } => {
            evidence.fallbacks.push(FallbackRecord {
                source,
                mechanism: "adjacent-face draft -> transition surface".to_string(),
                reason: format!("{warning:?}"),
            })
        }
        _ => evidence.warnings.push(KernelWarning {
            source,
            message: format!("{warning:?}"),
        }),
    }
}

/// Classify one sweep/loft/pipe/revolve warning.
pub(super) fn classify_sweep_warning(evidence: &mut FeatureEvidence, warning: &SweepWarning) {
    let source = RecordSource::Sweep;
    match warning {
        SweepWarning::FrenetFallback { .. } => evidence.fallbacks.push(FallbackRecord {
            source,
            mechanism: "Frenet frame -> Bishop frame".to_string(),
            reason: format!("{warning:?}"),
        }),
        SweepWarning::HelixApproximated => evidence.fallbacks.push(FallbackRecord {
            source,
            mechanism: "exact helix -> NURBS approximation".to_string(),
            reason: format!("{warning:?}"),
        }),
        _ => evidence.warnings.push(KernelWarning {
            source,
            message: format!("{warning:?}"),
        }),
    }
}

/// Classify one pattern instance outcome.
pub(super) fn classify_instance_outcome(
    evidence: &mut FeatureEvidence,
    index: usize,
    outcome: &InstanceOutcome,
) {
    let source = RecordSource::Pattern;
    match outcome {
        InstanceOutcome::Replayed => {}
        InstanceOutcome::GeometryFallback => evidence.fallbacks.push(FallbackRecord {
            source,
            mechanism: format!("instance {index}: parametric replay -> geometry copy-transform"),
            reason: "the instance's parametric replay failed".to_string(),
        }),
        InstanceOutcome::Skipped => evidence.skipped.push(SkipRecord {
            source,
            what: format!("instance {index}"),
            reason: "instance suppressed".to_string(),
        }),
        InstanceOutcome::Failed { error } => evidence.skipped.push(SkipRecord {
            source,
            what: format!("instance {index}"),
            reason: format!("could not be produced: {error}"),
        }),
    }
}

/// The variant name of a warning's `Debug` form (`RadiusClamped { .. }`
/// → `RadiusClamped`), for a mechanism label.
fn variant_name(warning: &BlendWarning) -> String {
    let debug = format!("{warning:?}");
    debug
        .split(|c: char| !c.is_alphanumeric())
        .next()
        .unwrap_or("")
        .to_string()
}
