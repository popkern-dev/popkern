//! Body accessors on [`FeatureTree`]: the live bodies, their BReps,
//! their features, and the host a Boss/Cut feature merged into.

use popkern_core::topology::BRep;

use crate::FeatureTree;
use crate::auto_merge;
use crate::node::{BodyId, FeatureId};

impl FeatureTree {
    // -- Body accessors -------------------------------------------------------

    /// Return all body IDs currently tracked by this tree.
    pub fn bodies(&self) -> Vec<BodyId> {
        self.body_manager.bodies()
    }

    /// Return the current BRep snapshot for a body.
    ///
    /// Returns `None` if no snapshot has been stored yet (i.e. the tree has
    /// not been replayed, or the body does not exist).
    pub fn body_brep(&self, body: BodyId) -> Option<&BRep> {
        self.body_manager.get_brep(body)
    }

    /// Return the features assigned to a body, in insertion order.
    ///
    /// Returns an empty [`Vec`] if the body does not exist.
    pub fn body_features(&self, body: BodyId) -> Vec<FeatureId> {
        self.body_manager.body_features(body)
    }

    /// The body a Boss/Cut feature merges into: the body whose BRep
    /// carries the face its sketch was placed on. `None` for a feature
    /// with no on-face sketch in its body (a new body, or a base-plane
    /// tool that is merged by an explicit boolean).
    pub fn merge_host(&self, feature: FeatureId) -> Option<BodyId> {
        if let Some(host) = self.merged_into.get(&feature) {
            return Some(*host);
        }
        let bodies: Vec<(BodyId, &BRep)> = self
            .body_manager
            .bodies()
            .into_iter()
            .filter_map(|b| self.body_manager.get_brep(b).map(|brep| (b, brep)))
            .collect();
        auto_merge::detect_host_body_among(feature, &self.features, bodies)
    }
}
