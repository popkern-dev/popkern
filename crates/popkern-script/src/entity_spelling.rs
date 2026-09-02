//! The canonical text form of an `EntityName` (spec §6.1 `named(...)`):
//! `{feature label}/{role category}/{index}`, with provenance
//! ancestors appended as `@n,n` when present. Labels make the form
//! readable and rename-safe; the emitter uses it wherever a document
//! references an entity without a recorded query.

use popkern_core::naming::{EntityName, EntityRole, NameId, RoleCategory};

use crate::error::ScriptError;
use crate::model::Model;
use crate::selection::EntityRef;

fn category_name(category: RoleCategory) -> &'static str {
    match category {
        RoleCategory::Modified => "Modified",
        RoleCategory::Generated => "Generated",
        RoleCategory::Lateral => "Lateral",
        RoleCategory::Cap => "Cap",
        RoleCategory::Intersection => "Intersection",
        RoleCategory::Junction => "Junction",
    }
}

fn parse_category(text: &str) -> Option<RoleCategory> {
    Some(match text {
        "Modified" => RoleCategory::Modified,
        "Generated" => RoleCategory::Generated,
        "Lateral" => RoleCategory::Lateral,
        "Cap" => RoleCategory::Cap,
        "Intersection" => RoleCategory::Intersection,
        "Junction" => RoleCategory::Junction,
        _ => return None,
    })
}

fn bad(text: &str, reason: &str) -> ScriptError {
    ScriptError::Argument {
        argument: "entity".into(),
        reason: format!("entity spelling '{text}': {reason}"),
    }
}

impl Model {
    /// Spell an entity name through its feature's label
    /// (`rounds/Generated/0`, `boolean-1/Modified/3@2,7`).
    pub fn spell_entity(&self, name: &EntityName) -> Result<String, ScriptError> {
        let node = self
            .tree
            .timeline()
            .iter()
            .filter_map(|&id| self.tree.feature(id))
            .find(|n| n.operation_id == name.operation)
            .ok_or_else(|| ScriptError::Argument {
                argument: "entity".into(),
                reason: format!(
                    "entity operation {:?} belongs to no feature of this document",
                    name.operation
                ),
            })?;
        let mut text = format!(
            "{}/{}/{}",
            node.label,
            category_name(name.role.category),
            name.role.index
        );
        if !name.ancestors.is_empty() {
            text.push('@');
            text.push_str(
                &name
                    .ancestors
                    .iter()
                    .map(|a| a.0.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
            );
        }
        Ok(text)
    }

    /// Parse the spelling produced by [`Model::spell_entity`].
    pub fn parse_entity(&self, text: &str) -> Result<EntityName, ScriptError> {
        let (head, ancestors) = match text.split_once('@') {
            Some((head, tail)) => (head, tail),
            None => (text, ""),
        };
        let mut parts = head.rsplitn(3, '/');
        let index = parts.next().ok_or_else(|| bad(text, "missing index"))?;
        let category = parts.next().ok_or_else(|| bad(text, "missing role"))?;
        let label = parts
            .next()
            .ok_or_else(|| bad(text, "expected label/Category/index"))?;
        let index: u32 = index
            .parse()
            .map_err(|_| bad(text, "index is not a number"))?;
        let category =
            parse_category(category).ok_or_else(|| bad(text, "unknown role category"))?;
        let id = self
            .tree
            .feature_by_label(label)
            .ok_or_else(|| popkern_history::HistoryError::LabelNotFound(label.to_string()))?;
        let operation = self
            .tree
            .feature(id)
            .map(|n| n.operation_id)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
        let ancestors = if ancestors.is_empty() {
            Vec::new()
        } else {
            ancestors
                .split(',')
                .map(|a| {
                    a.parse::<u32>()
                        .map(NameId)
                        .map_err(|_| bad(text, "ancestor is not a number"))
                })
                .collect::<Result<Vec<_>, _>>()?
        };
        Ok(EntityName {
            operation,
            role: EntityRole { category, index },
            ancestors,
        })
    }

    /// Resolve a label or a spelled name to an entity.
    pub(crate) fn resolve_entity_ref(
        &self,
        reference: EntityRef<'_>,
    ) -> Result<EntityName, ScriptError> {
        match reference {
            EntityRef::Label(label) => self.entity_by_label(label),
            EntityRef::Spelled(text) => self.parse_entity(text),
        }
    }

    /// Bind an entity label to a spelled entity (`label_entity(name,
    /// "rounds/Generated/0")`).
    pub fn label_spelled(&mut self, label: &str, text: &str) -> Result<EntityName, ScriptError> {
        let entity = self.parse_entity(text)?;
        self.tree.label_entity(label, entity.clone())?;
        Ok(entity)
    }
}
