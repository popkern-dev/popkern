//! Replay entry points and timeline scrubbing on [`FeatureTree`].
//!
//! Every replay resumes after the longest prefix of steps whose cached
//! records match the current tree (content-addressed, so an edit can
//! never be served stale geometry; see `cache` and `step_key`).
//! Scrubbing reads the [`DocumentState`] after any step: every alive
//! body, not only the active one.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::FeatureTree;
use crate::document_state::DocumentState;
use crate::error::{HistoryError, ReplayError};
use crate::node::FeatureId;
use crate::replay_report::ReplayResult;
use crate::undo::TreeCommand;

impl FeatureTree {
    // -- Replay ---------------------------------------------------------------

    /// Replay the feature tree up to the rollback bar (or the end of the
    /// timeline).
    ///
    /// Suppressed features are skipped. Failed features record an error
    /// but do not abort the replay — the bodies continue from the last
    /// good state. Steps whose inputs are unchanged since an earlier
    /// replay are not executed again: the replay resumes from the
    /// cached state after the longest unchanged prefix and republishes
    /// those steps' results and the tree's report rows, so the outcome
    /// equals a replay from scratch, report included.
    ///
    /// # Errors
    ///
    /// Returns [`ReplayError`] only for fatal infrastructure failures.
    /// Individual feature failures are recorded in
    /// [`ReplayResult::feature_results`].
    pub fn replay(&mut self, tolerance: Tolerance) -> Result<ReplayResult, ReplayError> {
        let end = self.rollback_end();
        self.replay_resuming(end, end, true, tolerance)
    }

    /// Replay re-executing every step from `id` on: at most the steps
    /// before `id` are served from the cache (from the nearest cached
    /// state at or before it). The returned report carries the executed
    /// steps' rows; [`FeatureTree::last_report`] keeps every feature's.
    ///
    /// # Errors
    ///
    /// Returns [`ReplayError`] if `id` is not in the timeline.
    pub fn replay_from(
        &mut self,
        id: FeatureId,
        tolerance: Tolerance,
    ) -> Result<ReplayResult, ReplayError> {
        let pos = self
            .timeline
            .iter()
            .position(|&fid| fid == id)
            .ok_or(HistoryError::FeatureNotFound)?;
        let end = self.rollback_end();
        self.replay_resuming(end, pos, false, tolerance)
    }

    /// Replay executing every step (only a frozen prefix is skipped).
    /// The reference a cached replay is checked against.
    ///
    /// # Errors
    ///
    /// As [`FeatureTree::replay`].
    pub fn replay_from_scratch(
        &mut self,
        tolerance: Tolerance,
    ) -> Result<ReplayResult, ReplayError> {
        let end = self.rollback_end();
        self.replay_resuming(end, 0, true, tolerance)
    }

    fn rollback_end(&self) -> usize {
        self.rollback_index
            .unwrap_or(self.timeline.len())
            .min(self.timeline.len())
    }

    // -- Rollback / scrubbing -------------------------------------------------

    /// Move the rollback bar to `index`.
    ///
    /// After this call the timeline is logically "rolled back" to just before
    /// the feature at `index`.  The next [`replay`](Self::replay) will stop at
    /// that point.
    ///
    /// Returns the body the feature just before the bar produced (see
    /// [`FeatureTree::brep_at`]) when the last replay's state there is
    /// cached.
    pub fn rollback_to(&mut self, index: usize) -> Option<BRep> {
        let old = self.rollback_index;
        self.rollback_index = Some(index);
        self.record_undo(TreeCommand::MoveRollback {
            old,
            new: Some(index),
        });
        let state = self.cached_state_at(index)?.clone();
        self.step_body(&state, index).cloned()
    }

    /// Remove the rollback bar, making all features in the timeline active.
    pub fn rollback_to_end(&mut self) {
        let old = self.rollback_index;
        self.rollback_index = None;
        self.record_undo(TreeCommand::MoveRollback { old, new: None });
    }

    /// Return the current rollback position, or `None` if no bar is set.
    pub fn rollback_position(&self) -> Option<usize> {
        self.rollback_index
    }

    /// The document after the first `index` timeline features of the
    /// LAST replay, when its state is cached (no execution; `index = 0`
    /// is the empty document). What a history scrubber reads.
    pub fn cached_state_at(&self, index: usize) -> Option<&DocumentState> {
        static EMPTY: std::sync::OnceLock<DocumentState> = std::sync::OnceLock::new();
        if index == 0 {
            return Some(EMPTY.get_or_init(DocumentState::empty));
        }
        self.cache.chain_state(index - 1)
    }

    /// The document after the first `index` timeline features of the
    /// CURRENT tree (`index = 0` is the empty document), ignoring the
    /// rollback bar.
    ///
    /// Served from the step cache when the steps are unchanged since
    /// they were cached; otherwise the missing steps execute from the
    /// nearest cached state (recorded in the cache, not published: the
    /// tree's current bodies and report keep the last replay's).
    ///
    /// # Errors
    ///
    /// [`HistoryError::InvalidPosition`] past the end of the timeline;
    /// [`ReplayError`] for fatal replay failures.
    pub fn state_at(
        &mut self,
        index: usize,
        tolerance: Tolerance,
    ) -> Result<DocumentState, ReplayError> {
        if index > self.timeline.len() {
            return Err(ReplayError::HistoryError(HistoryError::InvalidPosition(
                index,
            )));
        }
        if index == 0 {
            return Ok(DocumentState::empty());
        }
        self.evaluate_parameters_for_replay();
        let keys = self.step_keys(index, &tolerance);
        if let Some(state) = keys[index - 1]
            .and_then(|key| self.cache.get(key))
            .and_then(|record| record.state.clone())
        {
            return Ok(state);
        }
        let (from, start) = self.resume_point(&keys, index);
        Ok(self
            .replay_range(start, from, index, &keys, &tolerance)?
            .state)
    }

    /// The BRep after the feature at `index` in the timeline (one-based:
    /// `brep_at(1)` is after the first feature; `0` is an empty BRep):
    /// the body that feature built or modified — or, when it merged into
    /// a host body, that host — otherwise the document's active body.
    ///
    /// # Errors
    ///
    /// As [`FeatureTree::state_at`], and [`ReplayError::FeatureFailed`]
    /// when no body is alive after the step.
    pub fn brep_at(&mut self, index: usize, tolerance: Tolerance) -> Result<BRep, ReplayError> {
        if index == 0 {
            return Ok(BRep::new(tolerance));
        }
        let state = self.state_at(index, tolerance)?;
        self.step_body(&state, index)
            .cloned()
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: "brep_at".into(),
                error: "no body is alive after this step".into(),
            })
    }

    /// The body of `state` that the feature at one-based `index` built,
    /// modified, or merged into; otherwise the active body.
    fn step_body<'s>(&self, state: &'s DocumentState, index: usize) -> Option<&'s BRep> {
        let id = *self.timeline.get(index.checked_sub(1)?)?;
        let body = self.features.get(id)?.body;
        state
            .body(body)
            .or_else(|| self.merged_into.get(&id).and_then(|h| state.body(*h)))
            .or_else(|| state.active_body())
    }

    /// Plant hand-authored states as the last replay's chain, keyed as
    /// the current tree's steps (what the in-place route reads).
    #[cfg(test)]
    pub(crate) fn plant_test_chain(&mut self, steps: &[(FeatureId, BRep)]) {
        let keys = self.step_keys(steps.len(), &Tolerance::default());
        let planted = steps
            .iter()
            .zip(keys)
            .filter_map(|((id, brep), key)| {
                let body = self.features.get(*id)?.body;
                let state = DocumentState::empty().with_body(body, brep.clone(), *id);
                Some((*id, key?, state))
            })
            .collect();
        self.cache.plant_chain(planted);
    }

    /// The cached state after `id` under the CURRENT tree's keys.
    #[cfg(test)]
    pub(crate) fn current_cached_state(&self, id: FeatureId) -> Option<DocumentState> {
        let pos = self.timeline.iter().position(|&f| f == id)?;
        let key = self.step_keys(pos + 1, &Tolerance::default())[pos]?;
        self.cache.peek(key)?.state.clone()
    }
}
