//! The document label table: human-readable handles for features and
//! entities that are never the persistent key.
//!
//! Documents reference entities only by canonical
//! [`EntityName`](popkern_core::naming::EntityName) (provenance-derived)
//! and features only by [`FeatureId`]; labels are an indirection on
//! top. Renaming edits this table and nothing else — zero stored
//! references change, so a rename can never dangle a feature.
//!
//! Feature labels are document-unique. Every feature receives a label
//! when it is recorded: the requested one when it is free, otherwise an
//! auto-label (`extrude-1`, `fillet-3`) minted from a persisted per-
//! operation counter, so migration scripts and GUIs always have a
//! handle without the author having named anything. Counters never
//! reuse a number after a deletion. Entity labels live in a separate
//! namespace (`label("mounting_boss", selection)` and
//! `feature("fillet-3")` do not collide).

use std::collections::BTreeMap;

use popkern_core::naming::EntityName;
use serde::{Deserialize, Serialize};

use crate::FeatureTree;
use crate::error::HistoryError;
use crate::feature::Feature;
use crate::node::BodyId;
use crate::node::FeatureId;
use crate::undo::TreeCommand;

/// Persistent label state: entity labels and the auto-label counters.
///
/// Feature labels themselves live on each
/// [`FeatureNode`](crate::FeatureNode) (`label`); the tree keeps them
/// unique across the timeline.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct LabelTable {
    /// `label → EntityName` for user-labelled entities.
    #[serde(default)]
    pub entities: BTreeMap<String, EntityName>,
    /// Next auto-label number per operation label (`"fillet" → 4` means
    /// the next auto-labelled fillet is `fillet-4`).
    #[serde(default)]
    pub auto_counters: BTreeMap<String, u32>,
    /// `label → BodyId` for user-labelled bodies (tool models). A label
    /// names history, so it survives the body being consumed.
    #[serde(default)]
    pub bodies: BTreeMap<String, BodyId>,
}

impl LabelTable {
    /// Bind `label` to `entity` (`None` removes the label), without
    /// recording undo. Used by the undo log to replay label edits.
    pub(crate) fn bind_entity(&mut self, label: &str, entity: Option<EntityName>) {
        match entity {
            Some(name) => {
                self.entities.insert(label.to_string(), name);
            }
            None => {
                self.entities.remove(label);
            }
        }
    }
}

/// The qualified-parameter-name prefix a feature label produces
/// (`"Extrude 1"` → `extrude_1`, so `extrude_1.distance`).
pub fn qualified_prefix(label: &str) -> String {
    label.to_lowercase().replace(' ', "_")
}

impl FeatureTree {
    // -- Feature labels ------------------------------------------------------

    /// The feature carrying `label`, if any.
    pub fn feature_by_label(&self, label: &str) -> Option<FeatureId> {
        self.timeline
            .iter()
            .copied()
            .find(|&id| self.features.get(id).is_some_and(|n| n.label == label))
    }

    /// Every feature label in timeline order, paired with its feature.
    pub fn feature_labels(&self) -> Vec<(String, FeatureId)> {
        self.timeline
            .iter()
            .filter_map(|&id| self.features.get(id).map(|n| (n.label.clone(), id)))
            .collect()
    }

    /// Whether `label` is free for a feature other than `except`.
    fn feature_label_free(&self, label: &str, except: Option<FeatureId>) -> bool {
        !self.timeline.iter().any(|&id| {
            Some(id) != except && self.features.get(id).is_some_and(|n| n.label == label)
        })
    }

    /// The label a newly recorded feature receives: `requested` when it
    /// is non-empty and unused, otherwise a fresh auto-label for the
    /// feature's operation. Advances the auto-label counter only when
    /// an auto-label is minted.
    pub(crate) fn mint_feature_label(&mut self, feature: &Feature, requested: &str) -> String {
        if !requested.is_empty() && self.feature_label_free(requested, None) {
            return requested.to_string();
        }
        let operation = feature.operation_label();
        loop {
            let counter = self
                .labels
                .auto_counters
                .entry(operation.to_string())
                .or_insert(1);
            let candidate = format!("{operation}-{counter}");
            *counter += 1;
            if self.feature_label_free(&candidate, None) {
                return candidate;
            }
        }
    }

    /// Make every timeline feature's label unique, in timeline order:
    /// a feature whose label repeats an earlier one receives an
    /// auto-label. Run after loading documents written before labels
    /// were unique (format version 3). Deterministic for a given
    /// document.
    pub(crate) fn enforce_unique_feature_labels(&mut self) {
        let ids: Vec<FeatureId> = self.timeline.clone();
        for (pos, id) in ids.iter().enumerate() {
            let Some(label) = self.features.get(*id).map(|n| n.label.clone()) else {
                continue;
            };
            let repeats = label.is_empty()
                || ids[..pos]
                    .iter()
                    .any(|&earlier| self.features.get(earlier).is_some_and(|n| n.label == label));
            if repeats {
                let feature = self.features[*id].feature.clone();
                let fresh = self.mint_feature_label(&feature, "");
                self.relabel_feature(*id, &label, &fresh);
            }
        }
    }

    /// Set a feature's label and re-key its scoped parameters, without
    /// recording undo.
    pub(crate) fn relabel_feature(&mut self, id: FeatureId, old_label: &str, new_label: &str) {
        if let Some(node) = self.features.get_mut(id) {
            node.label = new_label.to_string();
            self.parameter_table.rename_feature_params(
                id,
                &qualified_prefix(old_label),
                &qualified_prefix(new_label),
            );
        }
    }

    /// Rename a feature.
    ///
    /// Edits the label only: the feature's [`FeatureId`], its
    /// [`OperationId`](popkern_core::naming::OperationId), every stored
    /// reference to it, and its geometry are untouched. Parameters
    /// scoped to the feature re-key from `old.field` to `new.field`.
    /// Records an undo command.
    ///
    /// # Errors
    ///
    /// [`HistoryError::FeatureNotFound`], [`HistoryError::EmptyLabel`],
    /// or [`HistoryError::DuplicateLabel`] when another feature already
    /// carries `label`.
    pub fn rename_feature(&mut self, id: FeatureId, label: &str) -> Result<(), HistoryError> {
        let old_label = self
            .features
            .get(id)
            .map(|n| n.label.clone())
            .ok_or(HistoryError::FeatureNotFound)?;
        if label.is_empty() {
            return Err(HistoryError::EmptyLabel);
        }
        if old_label == label {
            return Ok(());
        }
        if !self.feature_label_free(label, Some(id)) {
            return Err(HistoryError::DuplicateLabel(label.to_string()));
        }
        self.relabel_feature(id, &old_label, label);
        self.record_undo(TreeCommand::RenameFeature {
            id,
            old_label,
            new_label: label.to_string(),
        });
        Ok(())
    }

    // -- Entity labels -------------------------------------------------------

    /// The label table.
    pub fn labels(&self) -> &LabelTable {
        &self.labels
    }

    /// Bind `label` to an entity's canonical name. Records an undo
    /// command. Rebinding a label the same entity already carries is a
    /// no-op.
    ///
    /// # Errors
    ///
    /// [`HistoryError::EmptyLabel`]; [`HistoryError::DuplicateLabel`]
    /// when `label` already names a different entity.
    /// Bind a body label (`DuplicateLabel` if the label names another body).
    pub fn label_body(&mut self, label: &str, body: BodyId) -> Result<(), HistoryError> {
        if label.is_empty() {
            return Err(HistoryError::EmptyLabel);
        }
        match self.labels.bodies.get(label) {
            Some(existing) if *existing == body => return Ok(()),
            Some(_) => return Err(HistoryError::DuplicateLabel(label.to_string())),
            None => {}
        }
        self.labels.bodies.insert(label.to_string(), body);
        self.record_undo(TreeCommand::SetBodyLabel {
            label: label.to_string(),
            old: None,
            new: Some(body),
        });
        Ok(())
    }

    /// Remove a body label.
    pub fn unlabel_body(&mut self, label: &str) -> Result<BodyId, HistoryError> {
        let body = self
            .labels
            .bodies
            .remove(label)
            .ok_or_else(|| HistoryError::LabelNotFound(label.to_string()))?;
        self.record_undo(TreeCommand::SetBodyLabel {
            label: label.to_string(),
            old: Some(body),
            new: None,
        });
        Ok(body)
    }

    /// The body a label names.
    pub fn body_by_label(&self, label: &str) -> Option<BodyId> {
        self.labels.bodies.get(label).copied()
    }

    /// The label of a body, if any.
    pub fn body_label(&self, body: BodyId) -> Option<&str> {
        self.labels
            .bodies
            .iter()
            .find(|(_, b)| **b == body)
            .map(|(l, _)| l.as_str())
    }

    /// Bind an entity label (`DuplicateLabel` if the label names another entity).
    pub fn label_entity(&mut self, label: &str, entity: EntityName) -> Result<(), HistoryError> {
        if label.is_empty() {
            return Err(HistoryError::EmptyLabel);
        }
        match self.labels.entities.get(label) {
            Some(existing) if *existing == entity => return Ok(()),
            Some(_) => return Err(HistoryError::DuplicateLabel(label.to_string())),
            None => {}
        }
        self.labels.bind_entity(label, Some(entity.clone()));
        self.record_undo(TreeCommand::SetEntityLabel {
            label: label.to_string(),
            old: None,
            new: Some(entity),
        });
        Ok(())
    }

    /// Remove an entity label, returning the name it carried. Records
    /// an undo command.
    ///
    /// # Errors
    ///
    /// [`HistoryError::LabelNotFound`].
    pub fn unlabel_entity(&mut self, label: &str) -> Result<EntityName, HistoryError> {
        let entity = self
            .labels
            .entities
            .remove(label)
            .ok_or_else(|| HistoryError::LabelNotFound(label.to_string()))?;
        self.record_undo(TreeCommand::SetEntityLabel {
            label: label.to_string(),
            old: Some(entity.clone()),
            new: None,
        });
        Ok(entity)
    }

    /// Rename an entity label. The entity's canonical name is untouched.
    /// Records one undoable step.
    ///
    /// # Errors
    ///
    /// [`HistoryError::LabelNotFound`], [`HistoryError::EmptyLabel`],
    /// [`HistoryError::DuplicateLabel`].
    pub fn rename_entity_label(&mut self, old: &str, new: &str) -> Result<(), HistoryError> {
        if new.is_empty() {
            return Err(HistoryError::EmptyLabel);
        }
        if old == new {
            return self
                .labels
                .entities
                .contains_key(old)
                .then_some(())
                .ok_or_else(|| HistoryError::LabelNotFound(old.to_string()));
        }
        if self.labels.entities.contains_key(new) {
            return Err(HistoryError::DuplicateLabel(new.to_string()));
        }
        let entity = self
            .labels
            .entities
            .remove(old)
            .ok_or_else(|| HistoryError::LabelNotFound(old.to_string()))?;
        self.labels.bind_entity(new, Some(entity.clone()));
        self.record_undo(TreeCommand::Compound(vec![
            TreeCommand::SetEntityLabel {
                label: old.to_string(),
                old: Some(entity.clone()),
                new: None,
            },
            TreeCommand::SetEntityLabel {
                label: new.to_string(),
                old: None,
                new: Some(entity),
            },
        ]));
        Ok(())
    }

    /// The entity `label` names.
    pub fn entity_by_label(&self, label: &str) -> Option<&EntityName> {
        self.labels.entities.get(label)
    }

    /// The label bound to `entity`, if any.
    pub fn entity_label(&self, entity: &EntityName) -> Option<&str> {
        self.labels
            .entities
            .iter()
            .find(|(_, name)| *name == entity)
            .map(|(label, _)| label.as_str())
    }
}

#[cfg(test)]
#[path = "labels_tests.rs"]
mod tests;
