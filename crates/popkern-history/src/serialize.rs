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

use serde::{Deserialize, Serialize};

use crate::FeatureTree;
use crate::error::HistoryError;

/// Current document format version.
const CURRENT_VERSION: u32 = 3;

/// Versioned envelope used for serialization (writing).
///
/// The `tree` field is the full [`FeatureTree`] payload.
#[derive(Serialize)]
struct PknDocumentOut<'a> {
    /// Format version — currently always [`CURRENT_VERSION`].
    version: u32,
    /// The feature tree payload.
    tree: &'a FeatureTree,
}

/// Versioned envelope used for deserialization (reading).
#[derive(Deserialize)]
struct PknDocumentIn {
    /// Format version found in the file.
    version: u32,
    /// The feature tree payload.
    tree: FeatureTree,
}

impl FeatureTree {
    /// Serialize the feature tree to bytes (JSON format).
    ///
    /// The output is wrapped in a versioned document envelope at
    /// [`CURRENT_VERSION`].  Session-local state (cache, undo log, current
    /// BRep, body BReps) is **not** serialized.  After loading, call
    /// [`replay`](Self::replay) to rebuild all runtime state.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::SerializeError`] if JSON serialization fails.
    pub fn save(&self) -> Result<Vec<u8>, HistoryError> {
        let doc = PknDocumentOut {
            version: CURRENT_VERSION,
            tree: self,
        };
        serde_json::to_vec(&doc).map_err(|e| HistoryError::SerializeError(e.to_string()))
    }

    /// Deserialize a feature tree from bytes.
    ///
    /// Only version 3 is supported. Returns
    /// [`HistoryError::UnsupportedVersion`] for any other version.
    ///
    /// The returned tree has no BRep state — call [`replay`](Self::replay) to
    /// rebuild it.
    ///
    /// # Errors
    ///
    /// - [`HistoryError::SerializeError`] — malformed JSON or structural mismatch.
    /// - [`HistoryError::UnsupportedVersion`] — file was written by a different build.
    pub fn load(data: &[u8]) -> Result<Self, HistoryError> {
        let doc: PknDocumentIn = serde_json::from_slice(data)
            .map_err(|e| HistoryError::SerializeError(e.to_string()))?;

        if doc.version != CURRENT_VERSION {
            return Err(HistoryError::UnsupportedVersion {
                file_version: doc.version,
                max_supported: CURRENT_VERSION,
            });
        }

        let mut tree = doc.tree;
        tree.rebuild_indices();
        Ok(tree)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "serialize_tests.rs"]
mod tests;
