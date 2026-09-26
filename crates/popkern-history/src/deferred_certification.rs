//! Deferred seam certification for blend continuity evidence.
//!
//! A fillet's continuity evidence has two parts: sampled residuals,
//! cheap and always measured, and interval certificates that bound the
//! residuals over each whole seam — often 95% of a G2 fillet feature's
//! time. Construction never depends on the certificate, so a host that
//! must show the part first (the browser playground) sets
//! [`ContinuityCertification::Deferred`]: replay serves on sampled
//! evidence, marks the gated seams deferred, and keeps each such
//! feature's exact output BRep here. The host then certifies one seam
//! at a time with [`FeatureTree::certify_next_deferred_seam`] until
//! none remain, and the report upgrades in place.
//!
//! # Invariants
//!
//! - A feature has a kept BRep exactly when its latest report row
//!   carries deferred seams: every replayed feature drops its entry
//!   before executing and re-adds it only from that execution's
//!   evidence, the same per-feature granularity the report merges at.
//! - The kept BRep is the one the evidence was measured on, so its
//!   edge handles are valid there. The cache is never used for this:
//!   its strategy may evict or skip snapshots.
//! - The policy is session-local (never serialized); a loaded document
//!   certifies inline until a host opts in.

use std::collections::HashMap;

pub use popkern_blend::ContinuityCertification;
use popkern_core::topology::{BRep, EdgeHandle};

use crate::FeatureTree;
use crate::node::FeatureId;

/// Session-local deferred-certification state of a [`FeatureTree`].
#[derive(Debug, Clone, Default)]
pub struct DeferredCertification {
    pub(crate) policy: ContinuityCertification,
    /// Output BRep of each feature whose evidence has deferred seams.
    pub(crate) bodies: HashMap<FeatureId, BRep>,
}

/// One deferred seam certified by
/// [`FeatureTree::certify_next_deferred_seam`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeferredSeamProgress {
    /// The feature whose evidence gained the certificate.
    pub feature: FeatureId,
    /// The seam edge, on that feature's output BRep.
    pub edge: EdgeHandle,
    /// Deferred seams left across the whole last report.
    pub remaining: usize,
}

/// Why a deferred seam could not be certified.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum DeferredCertificationError {
    /// The report row has deferred seams but the BRep they were
    /// measured on is not held — the invariant above is broken.
    #[error("feature '{label}' has deferred seams but its measured BRep is not held")]
    BodyUnavailable {
        /// The feature.
        feature: FeatureId,
        /// Its label.
        label: String,
    },
}

impl FeatureTree {
    /// When blend seams are certified during replay (default
    /// [`ContinuityCertification::Inline`]).
    pub fn continuity_certification(&self) -> ContinuityCertification {
        self.deferred.policy
    }

    /// Set when blend seams are certified during replay. Takes effect
    /// at the next replay; already-deferred seams stay certifiable.
    pub fn set_continuity_certification(&mut self, policy: ContinuityCertification) {
        self.deferred.policy = policy;
    }

    /// Deferred seams across the last report.
    pub fn deferred_seam_count(&self) -> usize {
        self.last_report
            .features
            .iter()
            .filter_map(|row| row.blend.as_ref()?.continuity.as_ref())
            .map(|evidence| evidence.deferred_seams)
            .sum()
    }

    /// Certify the first deferred seam (timeline order) against the
    /// BRep its evidence was measured on, updating the last report in
    /// place. `Ok(None)` when nothing is deferred.
    ///
    /// # Errors
    ///
    /// [`DeferredCertificationError::BodyUnavailable`] when a row's
    /// measured BRep is missing (an invariant violation, reported
    /// rather than certified against a different BRep).
    pub fn certify_next_deferred_seam(
        &mut self,
    ) -> Result<Option<DeferredSeamProgress>, DeferredCertificationError> {
        let Some(row) = self.last_report.features.iter_mut().find(|row| {
            row.blend
                .as_ref()
                .and_then(|b| b.continuity.as_ref())
                .is_some_and(|e| e.deferred_seams > 0)
        }) else {
            return Ok(None);
        };
        let feature = row.feature;
        let Some(brep) = self.deferred.bodies.get(&feature) else {
            return Err(DeferredCertificationError::BodyUnavailable {
                feature,
                label: row.label.clone(),
            });
        };
        let certified = row
            .blend
            .as_mut()
            .and_then(|b| b.continuity.as_mut())
            .and_then(|evidence| evidence.certify_next_deferred(brep));
        let Some(certified) = certified else {
            return Ok(None);
        };
        if certified.remaining == 0 {
            self.deferred.bodies.remove(&feature);
        }
        Ok(Some(DeferredSeamProgress {
            feature,
            edge: certified.edge,
            remaining: self.deferred_seam_count(),
        }))
    }

    /// Replay bookkeeping: forget `feature`'s kept BRep before it
    /// executes again.
    pub(crate) fn forget_deferred_body(&mut self, feature: FeatureId) {
        self.deferred.bodies.remove(&feature);
    }

    /// Replay bookkeeping: keep `brep` when `feature`'s fresh evidence
    /// deferred any seam.
    pub(crate) fn keep_deferred_body(
        &mut self,
        feature: FeatureId,
        evidence: &crate::replay_report::FeatureEvidence,
        brep: &BRep,
    ) {
        let deferred = evidence
            .blend
            .as_ref()
            .and_then(|b| b.continuity.as_ref())
            .is_some_and(|e| e.deferred_seams > 0);
        if deferred {
            self.deferred.bodies.insert(feature, brep.clone());
        }
    }
}

#[cfg(test)]
#[path = "deferred_certification_tests.rs"]
mod tests;
