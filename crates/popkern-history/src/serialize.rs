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

use popkern_core::naming::NAMING_SCHEME_VERSION;

use crate::FeatureTree;
use crate::error::HistoryError;

/// Current document format version: what [`FeatureTree::save`] writes.
///
/// Version history (the compatibility policy is
/// `docs/decisions/2026-08-30-pkn-format-compatibility-policy.md`):
/// - **4** (2026-08-30): `naming_scheme` envelope field; label table,
///   auto-label counters and build tolerance on the tree; `RawBody`,
///   `Transform`, analytic sketches, draft literal planes, full fillet
///   recipes. Every field defaults, so version-3 documents load.
/// - **3**: the format before the scripted-modeling campaign.
pub const CURRENT_VERSION: u32 = 4;

/// Oldest document format version [`FeatureTree::load`] reads.
pub const MIN_SUPPORTED_VERSION: u32 = 3;

/// The naming scheme every document without a `naming_scheme` field was
/// written under: the field arrived while scheme 1 was the only scheme.
const UNRECORDED_NAMING_SCHEME: u32 = 1;

/// Versioned envelope used for serialization (writing).
///
/// The `tree` field is the full [`FeatureTree`] payload.
#[derive(Serialize)]
struct PknDocumentOut<'a> {
    /// Format version — currently always [`CURRENT_VERSION`].
    version: u32,
    /// Entity-naming scheme the stored references were minted under —
    /// always [`NAMING_SCHEME_VERSION`] when writing.
    naming_scheme: u32,
    /// The feature tree payload.
    tree: &'a FeatureTree,
}

/// Versioned envelope used for deserialization (reading).
#[derive(Deserialize)]
struct PknDocumentIn {
    /// Format version found in the file.
    version: u32,
    /// Naming scheme found in the file ([`UNRECORDED_NAMING_SCHEME`]
    /// when absent).
    #[serde(default = "unrecorded_naming_scheme")]
    naming_scheme: u32,
    /// The feature tree payload.
    tree: FeatureTree,
}

fn unrecorded_naming_scheme() -> u32 {
    UNRECORDED_NAMING_SCHEME
}

/// Check that the file's naming scheme is one this build reads.
///
/// Only scheme 1 has ever been minted, so there is no migration to run;
/// when a scheme 2 arrives this is where the scheme-1 → scheme-2 name
/// migration dispatches.
fn check_naming_scheme(file_scheme: u32) -> Result<(), HistoryError> {
    if file_scheme == NAMING_SCHEME_VERSION {
        Ok(())
    } else {
        Err(HistoryError::UnsupportedNamingScheme {
            file_scheme,
            supported: NAMING_SCHEME_VERSION,
        })
    }
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
            naming_scheme: NAMING_SCHEME_VERSION,
            tree: self,
        };
        serde_json::to_vec(&doc).map_err(|e| HistoryError::SerializeError(e.to_string()))
    }

    /// Deserialize a feature tree from bytes.
    ///
    /// Versions [`MIN_SUPPORTED_VERSION`]..=[`CURRENT_VERSION`] are read;
    /// older ones load through their field defaults (the additions
    /// since version 3 all default), and the tree is written back at
    /// [`CURRENT_VERSION`]. Any other version returns
    /// [`HistoryError::UnsupportedVersion`].
    ///
    /// The returned tree has no BRep state — call [`replay`](Self::replay) to
    /// rebuild it.
    ///
    /// # Errors
    ///
    /// - [`HistoryError::SerializeError`] — malformed JSON or structural mismatch.
    /// - [`HistoryError::UnsupportedVersion`] — file was written by a different build.
    /// - [`HistoryError::UnsupportedNamingScheme`] — the file's entity names
    ///   were minted under a scheme this build cannot read.
    pub fn load(data: &[u8]) -> Result<Self, HistoryError> {
        let doc: PknDocumentIn = serde_json::from_slice(data)
            .map_err(|e| HistoryError::SerializeError(e.to_string()))?;

        if !(MIN_SUPPORTED_VERSION..=CURRENT_VERSION).contains(&doc.version) {
            return Err(HistoryError::UnsupportedVersion {
                file_version: doc.version,
                min_supported: MIN_SUPPORTED_VERSION,
                max_supported: CURRENT_VERSION,
            });
        }

        check_naming_scheme(doc.naming_scheme)?;

        let mut tree = doc.tree;
        tree.rebuild_indices();
        tree.enforce_unique_feature_labels();
        Ok(tree)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "serialize_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "serialize_roundtrip_tests.rs"]
mod roundtrip_tests;

impl FeatureTree {
    /// Rebuild the `OperationId` -> `FeatureId` index from the features
    /// SlotMap.
    ///
    /// Called after deserialization to restore the transient `op_index`
    /// field.
    pub(crate) fn rebuild_indices(&mut self) {
        self.op_index.clear();
        for (fid, node) in &self.features {
            self.op_index.insert(node.operation_id, fid);
        }
    }
}
