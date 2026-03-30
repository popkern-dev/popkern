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
const CURRENT_VERSION: u32 = 1;

/// Versioned envelope wrapping the [`FeatureTree`] on disk.
///
/// Adding a version field up-front makes it possible to introduce
/// backward-incompatible schema changes in the future without breaking
/// existing `.pkn` files.
#[derive(Serialize, Deserialize)]
struct PknDocument {
    /// Format version — currently always [`CURRENT_VERSION`].
    version: u32,
    /// The feature tree payload.
    tree: FeatureTree,
}

impl FeatureTree {
    /// Serialize the feature tree to bytes (JSON format).
    ///
    /// The output is wrapped in a versioned [`PknDocument`] envelope.
    /// Session-local state (cache, undo log, current BRep, body BReps) is
    /// **not** serialized.  After loading, call [`replay`](Self::replay) to
    /// rebuild all runtime state.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::SerializeError`] if JSON serialization fails.
    pub fn save(&self) -> Result<Vec<u8>, HistoryError> {
        let doc = PknDocument {
            version: CURRENT_VERSION,
            tree: self.clone(),
        };
        serde_json::to_vec(&doc).map_err(|e| HistoryError::SerializeError(e.to_string()))
    }

    /// Deserialize a feature tree from bytes.
    ///
    /// Expects the versioned [`PknDocument`] envelope.  The returned tree has
    /// no BRep state — call [`replay`](Self::replay) to rebuild it.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::SerializeError`] if JSON deserialization fails,
    /// the data is structurally invalid, or the document version is
    /// unsupported.
    pub fn load(data: &[u8]) -> Result<Self, HistoryError> {
        let doc: PknDocument = serde_json::from_slice(data)
            .map_err(|e| HistoryError::SerializeError(e.to_string()))?;
        if doc.version != CURRENT_VERSION {
            return Err(HistoryError::SerializeError(format!(
                "unsupported document version {} (expected {})",
                doc.version, CURRENT_VERSION
            )));
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
