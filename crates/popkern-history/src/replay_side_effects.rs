//! What a feature execution reports back to the replay loop besides
//! its body: consumed bodies and the evidence record.

use crate::node::BodyId;
use crate::replay_report::FeatureEvidence;

/// Returned from [`execute_feature`] so that the caller can apply deferred
/// mutations (e.g. removing consumed tool bodies) only on success.
pub(crate) struct ExecuteSideEffects {
    /// Body IDs that were consumed by this feature (e.g. Boolean tool body).
    pub consumed_bodies: Vec<BodyId>,
    /// What the operation reported: warnings, blend evidence, fallback
    /// firings, skip records (the replay report channel).
    pub evidence: FeatureEvidence,
}

impl ExecuteSideEffects {
    /// No side effects and nothing reported.
    pub fn empty() -> Self {
        Self {
            consumed_bodies: Vec::new(),
            evidence: FeatureEvidence::empty(),
        }
    }

    /// No consumed bodies; the operation's evidence.
    pub fn reported(evidence: FeatureEvidence) -> Self {
        Self {
            consumed_bodies: Vec::new(),
            evidence,
        }
    }
}
