//! Multi-body tracking for the feature tree.
//!
//! [`BodyManager`] maintains the mapping from [`BodyId`]s to their associated
//! features and current BRep snapshots.  Each body represents a discrete solid
//! in a multi-body workflow (e.g. the target body and a boolean tool body).

use std::collections::HashMap;

use popkern_core::topology::BRep;

use crate::node::{BodyId, FeatureId};

/// Tracks multiple bodies within the feature tree.
///
/// Each body has an ordered list of features and a current BRep snapshot.
#[derive(Debug, Clone)]
pub(crate) struct BodyManager {
    /// Features assigned to each body.
    bodies: HashMap<BodyId, Vec<FeatureId>>,
    /// Current BRep for each body (from most recent replay).
    body_breps: HashMap<BodyId, BRep>,
    /// Next body id to assign.
    next_id: u32,
}

impl Default for BodyManager {
    fn default() -> Self {
        Self::new()
    }
}

impl BodyManager {
    /// Create a new, empty body manager.
    pub fn new() -> Self {
        Self {
            bodies: HashMap::new(),
            body_breps: HashMap::new(),
            next_id: 0,
        }
    }

    /// Create a new body and return its [`BodyId`].
    ///
    /// Body ids are assigned sequentially starting from 0.
    pub fn create_body(&mut self) -> BodyId {
        let id = BodyId(self.next_id);
        self.next_id += 1;
        self.bodies.insert(id, Vec::new());
        id
    }

    /// Assign a feature to a body.
    ///
    /// The feature is appended to the body's feature list.  If the body does
    /// not yet exist it is created on demand.
    pub fn assign_feature(&mut self, body: BodyId, feature: FeatureId) {
        self.bodies.entry(body).or_default().push(feature);
    }

    /// Store the current BRep for a body (called during replay).
    #[allow(dead_code)]
    pub fn set_brep(&mut self, body: BodyId, brep: BRep) {
        self.body_breps.insert(body, brep);
    }

    /// Return the current BRep for a body, if one has been stored.
    pub fn get_brep(&self, body: BodyId) -> Option<&BRep> {
        self.body_breps.get(&body)
    }

    /// Return all known body ids.
    pub fn bodies(&self) -> Vec<BodyId> {
        // Ascending id = creation order, and deterministic: the map is
        // a HashMap, and its iteration order once leaked into every
        // multi-body consumer — the whole-model STEP and thumbnail
        // concatenated bodies in a per-process random order, and
        // auto-merge host detection scanned candidates in it (PK-525).
        let mut ids: Vec<BodyId> = self.bodies.keys().copied().collect();
        ids.sort_unstable_by_key(|id| id.0);
        ids
    }

    /// Return the features assigned to a body, in insertion order.
    ///
    /// Returns an empty [`Vec`] if the body does not exist.
    pub fn body_features(&self, body: BodyId) -> Vec<FeatureId> {
        self.bodies.get(&body).cloned().unwrap_or_default()
    }

    /// Mark a body as consumed (e.g. the tool body after a boolean operation).
    ///
    /// Removes the body's BRep snapshot and feature list.  The [`BodyId`] will
    /// not be reused.
    #[allow(dead_code)]
    pub fn consume_body(&mut self, body: BodyId) {
        self.bodies.remove(&body);
        self.body_breps.remove(&body);
    }

    /// Clear all stored BRep snapshots (used when resetting for replay).
    pub fn clear_breps(&mut self) {
        self.body_breps.clear();
    }

    /// Clone every alive body's stored BRep snapshot — the complete
    /// end-of-timeline state of the last completed replay (consumed
    /// bodies are absent by [`FeatureTree::store_body_breps`]'s
    /// clear-then-store contract). The seed of the append fast path.
    pub fn all_breps(&self) -> std::collections::HashMap<BodyId, BRep> {
        self.body_breps.clone()
    }

    /// Clear all feature assignments (used when resetting for replay).
    ///
    /// Body ids are preserved; only the per-body feature lists are cleared.
    #[allow(dead_code)]
    pub fn clear_features(&mut self) {
        for features in self.bodies.values_mut() {
            features.clear();
        }
    }

    /// Resolve the body that a feature belongs to, by searching all bodies.
    ///
    /// Returns `None` if the feature is not assigned to any body.
    #[allow(dead_code)]
    pub fn body_of(&self, feature: FeatureId) -> Option<BodyId> {
        for (body_id, features) in &self.bodies {
            if features.contains(&feature) {
                return Some(*body_id);
            }
        }
        None
    }
}

#[cfg(test)]
#[path = "body_tests.rs"]
mod tests;
