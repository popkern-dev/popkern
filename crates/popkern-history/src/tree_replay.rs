//! Private replay helpers for [`FeatureTree`].
//!
//! Extracted from `lib.rs` to stay within the 500-line structural limit.

use std::collections::HashMap;

#[cfg(test)]
#[path = "replay_integration_tests.rs"]
mod integration_tests;

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::{FeatureResult, ReplayError};
use crate::node::BodyId;
use crate::{FeatureTree, ReplayRangeOutput, replay};

impl FeatureTree {
    /// Persist replay results: update body_manager and current_brep.
    pub(crate) fn store_body_breps(&mut self, body_breps: &HashMap<BodyId, BRep>) {
        for (body_id, brep) in body_breps {
            self.body_manager.set_brep(*body_id, brep.clone());
        }
        self.current_brep = body_breps
            .get(&BodyId(0))
            .cloned()
            .or_else(|| body_breps.values().next().cloned());
    }

    /// Execute features `from..to` using per-body BRep maps, caching results
    /// and returning the final body BReps together with per-feature outcomes.
    ///
    /// This is the shared core used by both [`replay`](FeatureTree::replay) and
    /// [`replay_from`](FeatureTree::replay_from).
    ///
    /// Each feature's body BRep is temporarily removed from the map before
    /// calling [`replay::execute_feature`] so that Boolean features can freely
    /// access (and consume) the tool body's BRep without borrow conflicts.
    pub(crate) fn replay_range(
        &mut self,
        body_breps: HashMap<BodyId, BRep>,
        from: usize,
        to: usize,
        tolerance: &Tolerance,
    ) -> Result<ReplayRangeOutput, ReplayError> {
        let mut body_breps = body_breps;
        let mut results = HashMap::new();

        for i in from..to {
            let id = self.timeline[i];
            let node = &self.features[id];

            if node.suppressed {
                results.insert(id, FeatureResult::Suppressed);
                continue;
            }

            // Evaluate conditional suppression. Clone suppress_if first to
            // release the immutable borrow on `node` before we need `self`
            // mutably below.
            let suppress_if = node.suppress_if.clone();
            let auto_suppressed = suppress_if
                .as_ref()
                .and_then(|expr| crate::expression::evaluate(expr, &self.parameter_table).ok())
                .map(|v| v != 0.0)
                .unwrap_or(false);

            if auto_suppressed {
                results.insert(id, FeatureResult::Suppressed);
                continue;
            }

            let mut ctx = TrackingContext::new(node.operation_id);
            let feature = node.feature.clone();
            let body = node.body;

            // Take this body's BRep out of the map so we can pass both
            // `&mut brep` and `&mut body_breps` to execute_feature.
            let mut brep = body_breps
                .remove(&body)
                .unwrap_or_else(|| BRep::new(*tolerance));

            // Snapshot BRep before execution for atomic rollback on failure.
            let brep_backup = brep.clone();

            let exec_result = replay::execute_feature(
                &mut brep,
                &feature,
                &mut ctx,
                tolerance,
                &mut body_breps,
                &self.parameter_table,
            );

            match exec_result {
                Ok(side_effects) => {
                    // Apply deferred side effects only on success.
                    for consumed in &side_effects.consumed_bodies {
                        body_breps.remove(consumed);
                    }
                    let history = ctx.finish();
                    results.insert(
                        id,
                        FeatureResult::Success {
                            history: history.clone(),
                        },
                    );
                    self.features[id].last_result = Some(FeatureResult::Success { history });
                    self.cache.store(id, brep.clone(), i);
                }
                Err(e) => {
                    // Atomically restore the pre-execution BRep so the body
                    // remains in a consistent state for subsequent features.
                    // Side effects NOT applied — tool body remains in the map.
                    brep = brep_backup;
                    let error_msg = e.to_string();
                    results.insert(
                        id,
                        FeatureResult::Failed {
                            error: error_msg.clone(),
                        },
                    );
                    self.features[id].last_result =
                        Some(FeatureResult::Failed { error: error_msg });
                }
            }

            // Put the BRep back (restored to backup on failure, or updated on success).
            body_breps.insert(body, brep);
        }

        Ok((body_breps, results))
    }
}
