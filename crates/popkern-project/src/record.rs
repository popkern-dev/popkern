//! Applied-migration records: the `project` field of a `.pkn`
//! envelope (live-projects-and-history.md §5.3) — which migrations
//! built the document, by content hash, and which features each one
//! recorded or edited — and how a saved document reconciles with the
//! folder it came from (§6.2).

use popkern_assembly::Document;
use popkern_history::FeatureTree;
use popkern_history::node::FeatureId;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::error::ProjectError;
use crate::manifest::PROJECT_FORMAT;
use crate::project::Project;

/// A feature a migration recorded or edited.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FeatureRef {
    /// The part definition holding the feature (assemblies only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub part: Option<String>,
    /// The feature's stable id.
    pub id: FeatureId,
    /// Its label when the migration ran (for people reading the record).
    pub label: String,
}

/// One applied migration.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MigrationRecord {
    /// The migration's name (file stem).
    pub name: String,
    /// SHA-256 of its normalised source.
    pub sha256: String,
    /// Features it recorded.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub features: Vec<FeatureRef>,
    /// Earlier features whose definition, suppression or timeline
    /// position it changed.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub edited: Vec<FeatureRef>,
}

/// The base document a project started from.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BaseRecord {
    /// Its project-relative path.
    pub path: String,
    /// SHA-256 of its bytes.
    pub sha256: String,
}

/// The `project` field: the document's applied migrations.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProjectRecord {
    /// The project-format version.
    pub format: u64,
    /// The base, when the project has one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub base: Option<BaseRecord>,
    /// Applied migrations in order.
    pub migrations: Vec<MigrationRecord>,
}

impl ProjectRecord {
    /// An empty record for `project` (its base, no migrations yet).
    pub fn new(project: &Project) -> ProjectRecord {
        ProjectRecord {
            format: PROJECT_FORMAT,
            base: project.base.as_ref().map(|b| BaseRecord {
                path: b.path.clone(),
                sha256: b.sha256.clone(),
            }),
            migrations: Vec::new(),
        }
    }
}

/// The `.pkn` bytes of `document` with `record` as its `project` field.
///
/// # Errors
///
/// [`ProjectError::Document`] when the document cannot be serialized.
pub fn to_pkn_bytes(document: &Document, record: &ProjectRecord) -> Result<Vec<u8>, ProjectError> {
    let bytes = document
        .to_pkn_bytes()
        .map_err(|e| ProjectError::Document {
            detail: e.to_string(),
        })?;
    let mut envelope: Value =
        serde_json::from_slice(&bytes).map_err(|e| ProjectError::Document {
            detail: e.to_string(),
        })?;
    let record = serde_json::to_value(record).map_err(|e| ProjectError::Document {
        detail: e.to_string(),
    })?;
    if let Value::Object(map) = &mut envelope {
        map.insert("project".into(), record);
    }
    serde_json::to_vec(&envelope).map_err(|e| ProjectError::Document {
        detail: e.to_string(),
    })
}

/// The `project` field of `.pkn` bytes (`None`: the document carries no
/// applied-migration record).
///
/// # Errors
///
/// [`ProjectError::Document`] for malformed JSON or a malformed record.
pub fn read_record(bytes: &[u8]) -> Result<Option<ProjectRecord>, ProjectError> {
    let envelope: Value = serde_json::from_slice(bytes).map_err(|e| ProjectError::Document {
        detail: e.to_string(),
    })?;
    match envelope.get("project") {
        None | Some(Value::Null) => Ok(None),
        Some(value) => serde_json::from_value(value.clone())
            .map(Some)
            .map_err(|e| ProjectError::Document {
                detail: format!("project record: {e}"),
            }),
    }
}

/// How a saved document relates to the folder.
#[derive(Debug, Clone)]
pub enum Reconcile {
    /// Every folder migration is applied.
    UpToDate,
    /// The record is a prefix of the folder: apply migrations from this
    /// index on.
    Pending {
        /// The first unapplied migration's index.
        from: usize,
    },
    /// The record's history differs from the folder's: rebuild.
    Diverged(ProjectError),
    /// The record applied a migration the folder lacks.
    Missing(ProjectError),
    /// The document carries no record: rebuild from the folder.
    Unrecorded,
}

/// Compare a saved document's record with the folder (§6.2).
pub fn reconcile(record: Option<&ProjectRecord>, project: &Project) -> Reconcile {
    let Some(record) = record else {
        return Reconcile::Unrecorded;
    };
    let base_now = project.base.as_ref().map(|b| b.sha256.as_str());
    let base_then = record.base.as_ref().map(|b| b.sha256.as_str());
    if base_now != base_then {
        return Reconcile::Diverged(ProjectError::ProvenanceDiverged {
            migration: "base".into(),
            recorded: base_then.unwrap_or("none").into(),
            current: base_now.unwrap_or("none").into(),
        });
    }
    for (i, applied) in record.migrations.iter().enumerate() {
        let Some(current) = project.migrations.get(i) else {
            return Reconcile::Missing(ProjectError::MissingMigration {
                name: applied.name.clone(),
            });
        };
        if current.sha256 != applied.sha256 {
            // Missing: the folder has neither the applied migration's
            // name nor its content anywhere (it was deleted). An edited
            // migration keeps its name; a moved one keeps its content —
            // both are a diverged history.
            let known = project
                .migrations
                .iter()
                .any(|m| m.name == applied.name || m.sha256 == applied.sha256);
            if !known {
                return Reconcile::Missing(ProjectError::MissingMigration {
                    name: applied.name.clone(),
                });
            }
            return Reconcile::Diverged(ProjectError::ProvenanceDiverged {
                migration: current.name.clone(),
                recorded: applied.sha256.clone(),
                current: current.sha256.clone(),
            });
        }
    }
    if record.migrations.len() == project.migrations.len() {
        Reconcile::UpToDate
    } else {
        Reconcile::Pending {
            from: record.migrations.len(),
        }
    }
}

/// Every part tree of `document`, with its definition name.
fn trees(document: &Document) -> Vec<(Option<String>, &FeatureTree)> {
    match document {
        Document::Part(model) => vec![(None, model.tree())],
        Document::Assembly(assembly) => assembly
            .definitions()
            .iter()
            .filter_map(|d| d.model().map(|m| (Some(d.name.clone()), m.tree())))
            .collect(),
    }
}

/// Each timeline feature's id, label, and what it is (definition,
/// suppression), in timeline order.
fn fingerprints(tree: &FeatureTree) -> Vec<(FeatureId, String, String)> {
    tree.timeline()
        .iter()
        .filter_map(|&id| {
            let node = tree.feature(id)?;
            let definition = serde_json::to_string(&node.feature).ok()?;
            // A parameter edit leaves the definition (which holds
            // parameter ids) unchanged: the values are part of what the
            // feature is.
            let table = tree.parameters();
            let params: Vec<String> = table
                .feature_params(id)
                .into_iter()
                .filter_map(|pid| table.get(pid))
                .map(|p| format!("{}={:?}/{:?}", p.value, p.unit, p.expression))
                .collect();
            let text = format!(
                "{definition}|{}|{:?}|{}",
                node.suppressed,
                node.suppress_if,
                params.join(",")
            );
            Some((id, node.label.clone(), text))
        })
        .collect()
}

/// The features `after` added and the earlier ones it changed, against
/// `before` (the document the migration started from). A feature is
/// edited when its definition or suppression changed, or when its
/// order among the features both documents share changed (a reorder;
/// an insertion before it is not an edit of it).
pub fn attribute(
    before: Option<&Document>,
    after: &Document,
) -> (Vec<FeatureRef>, Vec<FeatureRef>) {
    let mut added = Vec::new();
    let mut edited = Vec::new();
    let earlier = before.map(trees).unwrap_or_default();
    for (part, tree) in trees(after) {
        let then = earlier
            .iter()
            .find(|(p, _)| *p == part)
            .map(|(_, t)| fingerprints(t))
            .unwrap_or_default();
        let now = fingerprints(tree);
        let shared_then: Vec<FeatureId> = then
            .iter()
            .map(|(id, _, _)| *id)
            .filter(|id| now.iter().any(|(n, _, _)| n == id))
            .collect();
        let shared_now: Vec<FeatureId> = now
            .iter()
            .map(|(id, _, _)| *id)
            .filter(|id| shared_then.contains(id))
            .collect();
        for (id, label, text) in now {
            let reference = FeatureRef {
                part: part.clone(),
                id,
                label,
            };
            match then.iter().find(|(old, _, _)| *old == id) {
                None => added.push(reference),
                Some((_, _, old_text)) => {
                    let rank = |order: &[FeatureId]| order.iter().position(|f| *f == id);
                    if *old_text != text || rank(&shared_then) != rank(&shared_now) {
                        edited.push(reference);
                    }
                }
            }
        }
    }
    (added, edited)
}

#[cfg(test)]
#[path = "record_tests.rs"]
mod tests;
