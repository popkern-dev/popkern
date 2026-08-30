//! Document-level attributes: the build tolerance and the naming scheme.
//!
//! Geometric tolerance is always an explicit parameter (design invariant
//! 5). A document records the tolerance it was built with so that
//! opening it anywhere rebuilds the same body: [`FeatureTree::replay`]
//! still takes the tolerance explicitly (scripts override per ephemeral
//! build), and [`FeatureTree::replay_document`] passes the stored one.
//! Every [`ReplayResult`](crate::ReplayResult) states the tolerance
//! that produced it.

use popkern_core::tolerance::Tolerance;

use crate::FeatureTree;
use crate::error::ReplayError;
use crate::replay_report::ReplayResult;
use crate::undo::TreeCommand;

impl FeatureTree {
    /// The document's build tolerance, if one has been set.
    pub fn build_tolerance(&self) -> Option<Tolerance> {
        self.build_tolerance
    }

    /// The tolerance [`replay_document`](Self::replay_document) uses:
    /// the stored attribute, or [`Tolerance::default`] when unset.
    pub fn effective_build_tolerance(&self) -> Tolerance {
        self.build_tolerance.unwrap_or_default()
    }

    /// Set or clear the document build tolerance. Records an undo
    /// command. Cached replay snapshots are evicted, since they were
    /// built under the previous tolerance.
    pub fn set_build_tolerance(&mut self, tolerance: Option<Tolerance>) {
        let old = self.build_tolerance;
        if old == tolerance {
            return;
        }
        self.build_tolerance = tolerance;
        self.cache.clear();
        self.record_undo(TreeCommand::SetBuildTolerance {
            old,
            new: tolerance,
        });
    }

    /// Replay the whole tree under the document's own build tolerance
    /// (see [`effective_build_tolerance`](Self::effective_build_tolerance)).
    ///
    /// # Errors
    ///
    /// As [`replay`](Self::replay).
    pub fn replay_document(&mut self) -> Result<ReplayResult, ReplayError> {
        let tolerance = self.effective_build_tolerance();
        self.replay(tolerance)
    }
}
