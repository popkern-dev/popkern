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
//! # Versions
//!
//! | Version | Description |
//! |---------|-------------|
//! | 1       | Feature dimensions stored as bare `f64` fields. |
//! | 2       | Feature dimensions stored as [`crate::parameter::ParameterId`] handles; full [`crate::parameter::ParameterTable`] persisted. |
//! | 3       | `Extrude` uses `end_condition: ExtrudeEnd` and `Revolve` uses `axis: AxisDesc` + `end_condition: RevolveEnd`. |
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

#[path = "serialize_v1.rs"]
mod serialize_v1;

#[path = "serialize_v2.rs"]
mod serialize_v2;

use serde::{Deserialize, Serialize};
use serde_json::Value;

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
///
/// The `tree` field is kept as a raw [`Value`] so that [`load`] can dispatch
/// to a version-specific decoder before attempting structural parsing.
///
/// [`load`]: FeatureTree::load
#[derive(Deserialize)]
struct PknDocumentIn {
    /// Format version found in the file.
    version: u32,
    /// Raw JSON of the tree payload.
    tree: Value,
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
    /// Dispatches to a version-specific loader:
    ///
    /// - **v1** — migrates bare `f64` dimensions to a new
    ///   [`crate::parameter::ParameterTable`].  See [`serialize_v1`].
    /// - **v2** — migrates old flat `Extrude`/`Revolve` fields to the v3
    ///   `ExtrudeEnd`/`AxisDesc`/`RevolveEnd` wrappers.  See [`serialize_v2`].
    /// - **v3** — deserializes directly; rebuilds transient indices.
    /// - **other** — returns [`HistoryError::UnsupportedVersion`].
    ///
    /// The returned tree has no BRep state — call [`replay`](Self::replay) to
    /// rebuild it.
    ///
    /// # Errors
    ///
    /// - [`HistoryError::SerializeError`] — malformed JSON or structural mismatch.
    /// - [`HistoryError::UnsupportedVersion`] — file was written by a newer build.
    pub fn load(data: &[u8]) -> Result<Self, HistoryError> {
        let doc: PknDocumentIn = serde_json::from_slice(data)
            .map_err(|e| HistoryError::SerializeError(e.to_string()))?;

        match doc.version {
            1 => serialize_v1::migrate_v1(doc.tree),
            2 => serialize_v2::migrate_v2(doc.tree),
            3 => {
                let mut tree: FeatureTree = serde_json::from_value(doc.tree)
                    .map_err(|e| HistoryError::SerializeError(e.to_string()))?;
                tree.rebuild_indices();
                Ok(tree)
            }
            v => Err(HistoryError::UnsupportedVersion {
                file_version: v,
                max_supported: CURRENT_VERSION,
            }),
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "serialize_tests.rs"]
mod tests;
