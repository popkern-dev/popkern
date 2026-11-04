//! What a feature execution reports back to the replay loop besides
//! its body: consumed bodies and the evidence record.

use popkern_core::topology::BRep;

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
    /// Other bodies this feature rewrote (the source of separated
    /// pieces), stored in place of their BReps on success.
    pub written_bodies: Vec<(BodyId, BRep)>,
}

impl ExecuteSideEffects {
    /// No side effects and nothing reported.
    pub fn empty() -> Self {
        Self {
            consumed_bodies: Vec::new(),
            evidence: FeatureEvidence::empty(),
            written_bodies: Vec::new(),
        }
    }

    /// No consumed bodies; the operation's evidence.
    pub fn reported(evidence: FeatureEvidence) -> Self {
        Self {
            consumed_bodies: Vec::new(),
            evidence,
            written_bodies: Vec::new(),
        }
    }
}
