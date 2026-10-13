//! Replay report channel — the observability spine of a feature-tree
//! replay (scripted-modeling design §6.2, §8).
//!
//! Every feature execution returns a [`FeatureEvidence`] beside its
//! BRep mutation: the kernel's non-fatal warnings, the measured blend
//! continuity record, explicit fallback firings (a construction that
//! served through its weaker path), and skip records (work the
//! operation declined to do while still serving). The replay loop
//! wraps each one in a [`FeatureReport`] with the feature's identity,
//! outcome, and wall time, and aggregates them in timeline order into
//! a [`ReplayReport`] returned in [`ReplayResult::report`] and kept on
//! the tree ([`crate::FeatureTree::last_report`]).
//!
//! No print-trace artifacts: everything a script, GUI, daemon, or MCP
//! client learns about a replay arrives here as structured data. The
//! records are serde-derived so the same payload crosses the
//! daemon/MCP boundary unchanged.
//!
//! # Honesty rules
//!
//! - A blend feature ALWAYS carries measured continuity evidence
//!   (fillets: `Some`; chamfers: `None` — their boundary is a designed
//!   crease and nothing smooth was requested). The evidence is a record
//!   of the measurement at blend time; its handles are provenance, not
//!   live references.
//! - Fallback firings are recorded once, as [`FallbackRecord`]s, not as
//!   warnings — a downgraded route is a first-class fact, not a note.
//! - The boolean cellular→graph-merge registry
//!   ([`popkern_bool::merge_fallback`]) is process-global; the replay
//!   loop drains it immediately before and after each boolean so a
//!   record is attributed to the feature whose boolean fired it. A
//!   concurrent boolean in another thread of the same process during
//!   that window would be attributed here too — the registry is a
//!   process census, and this channel reports it as observed.

use std::collections::HashMap;
use std::time::Duration;

use popkern_blend::{BlendContinuityEvidence, ContinuityLevel};
use serde::{Deserialize, Serialize};

use crate::error::FeatureResult;
use crate::node::FeatureId;
use crate::parameter::{EvaluationWarning, ParameterId};

#[path = "replay_report_evidence.rs"]
mod evidence;
#[path = "replay_report_corners.rs"]
mod replay_report_corners;
pub(crate) use evidence::{FeatureEvidence, face_census};

/// Seam sampling density for blend-time continuity evidence: the
/// commercial gate default (the apex fixtures use 100), and the same
/// density `popkern-part` records so the two channels agree.
pub const EVIDENCE_SAMPLES_PER_SEAM: usize = 33;

/// The result of replaying a feature tree.
///
/// Contains per-feature outcomes, the replay report, and the total
/// elapsed time.
#[derive(Debug)]
pub struct ReplayResult {
    /// Per-feature outcomes, keyed by [`FeatureId`].
    pub feature_results: HashMap<FeatureId, FeatureResult>,
    /// Total replay wall time, or `None` on platforms without a
    /// monotonic clock (see `popkern_core::time::Stopwatch`).
    pub elapsed: Option<Duration>,
    /// Warnings from expression evaluation (math errors, bounds clamping).
    pub eval_warnings: Vec<(ParameterId, EvaluationWarning)>,
    /// Per-feature evidence in timeline order. A whole-document
    /// [`replay`](crate::FeatureTree::replay) reports every feature
    /// (steps served from the step cache republish their recorded
    /// rows); [`replay_from`](crate::FeatureTree::replay_from) reports
    /// only the steps it executed, and the tree's
    /// [`last_report`](crate::FeatureTree::last_report) keeps the rest.
    pub report: ReplayReport,
    /// The tolerance this replay built under.
    pub tolerance: popkern_core::tolerance::Tolerance,
    /// How many leading timeline steps were served from the step cache
    /// (their inputs were unchanged) instead of executing; the replay
    /// executed every step after them.
    pub cached_steps: usize,
}

/// The final document, per-feature outcomes, and the report from a
/// replay pass.
pub(crate) struct ReplayRangeOutput {
    /// The document after the last executed step (bodies shared with
    /// the step cache; the published bodies are its handles).
    pub state: crate::document_state::DocumentState,
    /// Per-feature outcomes.
    pub results: HashMap<FeatureId, FeatureResult>,
    /// Per-feature evidence, timeline order.
    pub report: ReplayReport,
}

/// Aggregated per-feature evidence for one replay, in timeline order.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ReplayReport {
    /// One report per executed (or suppressed) feature, timeline order.
    pub features: Vec<FeatureReport>,
}

impl ReplayReport {
    /// The report for `feature`, if this replay covered it.
    pub fn feature(&self, feature: FeatureId) -> Option<&FeatureReport> {
        self.features.iter().find(|r| r.feature == feature)
    }

    /// Every kernel warning with the feature that produced it.
    pub fn warnings(&self) -> impl Iterator<Item = (&FeatureReport, &KernelWarning)> {
        self.features
            .iter()
            .flat_map(|r| r.warnings.iter().map(move |w| (r, w)))
    }

    /// Every fallback firing with the feature that produced it.
    pub fn fallbacks(&self) -> impl Iterator<Item = (&FeatureReport, &FallbackRecord)> {
        self.features
            .iter()
            .flat_map(|r| r.fallbacks.iter().map(move |f| (r, f)))
    }

    /// Every skip record with the feature that produced it.
    pub fn skips(&self) -> impl Iterator<Item = (&FeatureReport, &SkipRecord)> {
        self.features
            .iter()
            .flat_map(|r| r.skipped.iter().map(move |s| (r, s)))
    }

    /// Blend features whose measured continuity is below what was
    /// requested — the silent-downgrade census.
    pub fn downgraded_blends(&self) -> impl Iterator<Item = &FeatureReport> {
        self.features.iter().filter(|r| {
            r.blend
                .as_ref()
                .is_some_and(BlendEvidenceRecord::downgraded)
        })
    }

    /// Features that failed to execute.
    pub fn failures(&self) -> impl Iterator<Item = &FeatureReport> {
        self.features
            .iter()
            .filter(|r| matches!(r.outcome, ReplayOutcome::Failed { .. }))
    }

    /// `true` when every feature served with nothing to report: no
    /// failures, warnings, fallbacks, skips, or continuity downgrades.
    /// Suppressed features are not findings.
    pub fn is_clean(&self) -> bool {
        self.features.iter().all(FeatureReport::is_clean)
    }

    /// Merge `newer` into this report: features it covers replace their
    /// previous record, everything else is kept, and the result is
    /// ordered by `timeline`. Records for features no longer in the
    /// timeline are dropped.
    pub(crate) fn absorb(&mut self, newer: ReplayReport, timeline: &[FeatureId]) {
        let mut by_id: HashMap<FeatureId, FeatureReport> =
            self.features.drain(..).map(|r| (r.feature, r)).collect();
        for report in newer.features {
            by_id.insert(report.feature, report);
        }
        self.features = timeline.iter().filter_map(|id| by_id.remove(id)).collect();
    }
}

/// Everything one feature's execution reported.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FeatureReport {
    /// The feature this report describes.
    pub feature: FeatureId,
    /// The feature's user-visible label.
    pub label: String,
    /// The operation kind (`"fillet"`, `"boolean"`, `"linear_pattern"`,
    /// …) — the same label the feature's structured errors carry.
    pub operation: String,
    /// Whether the feature served, failed, or was suppressed.
    pub outcome: ReplayOutcome,
    /// Wall time spent executing the feature (kernel operation plus
    /// any auto-merge), or `None` when it was not executed or the
    /// platform has no monotonic clock.
    pub elapsed: Option<Duration>,
    /// Non-fatal kernel warnings that are neither fallbacks nor skips.
    pub warnings: Vec<KernelWarning>,
    /// Measured blend evidence — `Some` for every served fillet or
    /// chamfer, `None` for every other feature.
    pub blend: Option<BlendEvidenceRecord>,
    /// Constructions that served through their explicit weaker path.
    pub fallbacks: Vec<FallbackRecord>,
    /// Work the operation declined while still serving.
    pub skipped: Vec<SkipRecord>,
}

impl FeatureReport {
    /// A report for a feature that was not executed (suppressed, or
    /// executed with nothing collected yet).
    pub(crate) fn new(
        feature: FeatureId,
        label: impl Into<String>,
        operation: impl Into<String>,
        outcome: ReplayOutcome,
    ) -> Self {
        Self {
            feature,
            label: label.into(),
            operation: operation.into(),
            outcome,
            elapsed: None,
            warnings: Vec::new(),
            blend: None,
            fallbacks: Vec::new(),
            skipped: Vec::new(),
        }
    }

    /// Attach the evidence an execution collected.
    pub(crate) fn with_evidence(mut self, evidence: FeatureEvidence) -> Self {
        self.warnings.extend(evidence.warnings);
        if evidence.blend.is_some() {
            self.blend = evidence.blend;
        }
        self.fallbacks.extend(evidence.fallbacks);
        self.skipped.extend(evidence.skipped);
        self
    }

    /// `true` when the feature served (or was suppressed) with nothing
    /// to report.
    pub fn is_clean(&self) -> bool {
        !matches!(self.outcome, ReplayOutcome::Failed { .. })
            && self.warnings.is_empty()
            && self.fallbacks.is_empty()
            && self.skipped.is_empty()
            && !self
                .blend
                .as_ref()
                .is_some_and(BlendEvidenceRecord::downgraded)
    }
}

/// How a feature's execution ended.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ReplayOutcome {
    /// The feature executed and its body was committed.
    Served,
    /// The feature failed; the body was restored to its pre-feature
    /// state (transactional per feature).
    Failed {
        /// The structured failure payload.
        error: crate::failure::FeatureFailure,
    },
    /// The feature was not executed.
    Suppressed {
        /// Why it was skipped.
        reason: SuppressionReason,
    },
}

/// Why a feature was not executed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SuppressionReason {
    /// The node's own `suppressed` flag.
    Explicit,
    /// The active configuration overrides the node's suppression.
    Configuration,
    /// The node's `suppress_if` expression evaluated non-zero.
    Conditional,
    /// The feature depends (by entity or parameter) on a Reference
    /// feature that failed to evaluate and has no last-known value.
    DependsOnFailedReference,
}

/// Which kernel subsystem produced a record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RecordSource {
    /// Fillet / chamfer (`popkern-blend`).
    Blend,
    /// Shell / thicken / offset (`popkern-shell`).
    Shell,
    /// Draft (`popkern-draft`).
    Draft,
    /// Sweep / loft / pipe / revolve (`popkern-sweep`).
    Sweep,
    /// Boolean features and auto-merge (`popkern-bool`).
    Boolean,
    /// Pattern and mirror instancing (this crate).
    Pattern,
    /// Reference resolution: a recorded name rebound by provenance
    /// (this crate, PK-677).
    Naming,
}

/// A non-fatal kernel warning that is neither a fallback nor a skip.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KernelWarning {
    /// The subsystem that emitted it.
    pub source: RecordSource,
    /// The kernel warning's own rendering — `Display` where the
    /// subsystem defines one, its structured `Debug` form otherwise —
    /// with the entity handles it names.
    pub message: String,
}

/// A construction that served through its explicit weaker path.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FallbackRecord {
    /// The subsystem that fell back.
    pub source: RecordSource,
    /// The boundary that was crossed, in mechanism terms
    /// (`"cellular selection -> graph-merge selection"`).
    pub mechanism: String,
    /// Why the stronger path did not apply.
    pub reason: String,
}

/// Work an operation declined while still serving the rest.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SkipRecord {
    /// The subsystem that skipped.
    pub source: RecordSource,
    /// What was skipped (an edge, a face, a pattern instance).
    pub what: String,
    /// Why.
    pub reason: String,
}

/// The measured record of one served blend feature.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlendEvidenceRecord {
    /// Edges the feature asked to blend.
    pub edges_requested: usize,
    /// Corner patches (vertex blends) the kernel inserted.
    pub corner_patches: usize,
    /// Continuity the feature asked for; chamfers record `G0`.
    pub requested_continuity: ContinuityLevel,
    /// Measured continuity over the blend's live seams: `Some` for a
    /// fillet, `None` for a chamfer (nothing smooth was requested).
    pub continuity: Option<BlendContinuityEvidence>,
    /// What closed each filleted corner (PK-603), in corner order.
    #[serde(default)]
    pub corners: Vec<CornerEvidence>,
}

/// What closed one filleted corner (PK-603).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CornerEvidence {
    /// Filleted edges meeting at the corner.
    pub sides: usize,
    /// The construction that served: `setback` (the KP20 cap or its G1
    /// downgrade), `fitted` (the variational fit), or
    /// `classic:<method>` (`classic:gregory`, `classic:rolling_ball_sphere`,
    /// `classic:exact_boundary`, `classic:mitre`, …).
    pub construction: String,
    /// Continuity the corner's construction reports achieving.
    pub achieved: ContinuityLevel,
    /// Faces the corner construction authored (zero for a mitre).
    pub patches: usize,
    /// Mean of the corner faces' boundary vertices, for locating the
    /// corner; `None` when no face was authored.
    pub centroid: Option<[f64; 3]>,
}

impl BlendEvidenceRecord {
    /// `true` when the measured continuity is below the request.
    pub fn downgraded(&self) -> bool {
        self.continuity
            .as_ref()
            .is_some_and(BlendContinuityEvidence::downgraded)
    }
}

#[cfg(test)]
#[path = "replay_report_tests.rs"]
mod tests;
