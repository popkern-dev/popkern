//! Serialization and deserialization for the parametric feature tree.
//!
//! The [`FeatureTree`] is the document format: BRep state is NOT saved
//! because it can be deterministically rebuilt by calling
//! [`FeatureTree::replay`] after loading.  Only the structural data
//! (features, timeline, rollback index, operation-id counter) is persisted.
//!
//! # Format
//!
//! Documents are stored as **JSON** (`.pkn` files).  The format is
//! human-readable and git-friendly.  Future versions may offer a compact
//! binary encoding via the `bincode` dev-dependency.
//!
//! # Round-trip
//!
//! ```no_run
//! # use popkern_history::FeatureTree;
//! # use popkern_core::tolerance::Tolerance;
//! let mut tree = FeatureTree::new();
//! // … add features …
//! let bytes = tree.save().expect("save failed");
//! let mut loaded = FeatureTree::load(&bytes).expect("load failed");
//! loaded.replay(Tolerance::default()).expect("replay failed");
//! ```

use crate::FeatureTree;
use crate::error::HistoryError;

impl FeatureTree {
    /// Serialize the feature tree to bytes (JSON format).
    ///
    /// Session-local state (cache, undo log, current BRep, body BReps) is
    /// **not** serialized.  After loading, call [`replay`](Self::replay) to
    /// rebuild all runtime state.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::SerializeError`] if JSON serialization fails.
    pub fn save(&self) -> Result<Vec<u8>, HistoryError> {
        serde_json::to_vec(self).map_err(|e| HistoryError::SerializeError(e.to_string()))
    }

    /// Deserialize a feature tree from bytes.
    ///
    /// The returned tree has no BRep state — call [`replay`](Self::replay) to
    /// rebuild it.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::SerializeError`] if JSON deserialization fails
    /// or the data is structurally invalid.
    pub fn load(data: &[u8]) -> Result<Self, HistoryError> {
        serde_json::from_slice(data).map_err(|e| HistoryError::SerializeError(e.to_string()))
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "serialize_tests.rs"]
mod tests;
