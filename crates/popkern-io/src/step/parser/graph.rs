//! Entity graph — indexed access to parsed STEP entities.
//!
//! Builds a [`HashMap`]-backed lookup from entity ID to [`EntityInstance`],
//! validating that all entity references resolve to existing entries.

use std::collections::HashMap;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::{EntityId, EntityInstance, Parameter};

/// An indexed graph of STEP entity instances.
///
/// Provides O(1) lookup by [`EntityId`] and validates that all inter-entity
/// references resolve at build time.
#[derive(Debug)]
pub struct EntityGraph {
    entities: HashMap<EntityId, EntityInstance>,
}

impl EntityGraph {
    /// Build from parsed entities.
    ///
    /// Validates that all entity references point to existing entries.
    /// Duplicate IDs and broken references are logged via the provided
    /// [`DiagnosticLog`].
    pub fn build(entities: Vec<EntityInstance>, log: &mut DiagnosticLog) -> Self {
        let mut map = HashMap::with_capacity(entities.len());
        for entity in entities {
            let id = entity.id;
            if map.contains_key(&id) {
                log.warn(
                    Phase::Resolve,
                    Some(id),
                    format!("duplicate entity #{id}, overwriting"),
                );
            }
            map.insert(id, entity);
        }

        // Validate references
        let ids: Vec<EntityId> = map.keys().copied().collect();
        for &id in &ids {
            if let Some(entity) = map.get(&id) {
                let refs = collect_entity_refs(&entity.params);
                for r in refs {
                    if !map.contains_key(&r) {
                        log.error(
                            Phase::Resolve,
                            Some(id),
                            format!("broken reference to #{r}"),
                        );
                    }
                }
            }
        }

        Self { entities: map }
    }

    /// Look up an entity by ID.
    pub fn get(&self, id: EntityId) -> Option<&EntityInstance> {
        self.entities.get(&id)
    }

    /// Number of entities in the graph.
    pub fn len(&self) -> usize {
        self.entities.len()
    }

    /// True if the graph contains no entities.
    pub fn is_empty(&self) -> bool {
        self.entities.is_empty()
    }

    /// Iterate over all `(id, entity)` pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&EntityId, &EntityInstance)> {
        self.entities.iter()
    }
}

/// Recursively collect all entity reference IDs from a parameter list.
fn collect_entity_refs(params: &[Parameter]) -> Vec<EntityId> {
    let mut refs = Vec::new();
    for p in params {
        match p {
            Parameter::EntityRef(id) => refs.push(*id),
            Parameter::List(items) => refs.extend(collect_entity_refs(items)),
            Parameter::Typed(_, inner) => match inner.as_ref() {
                Parameter::List(items) => refs.extend(collect_entity_refs(items)),
                Parameter::EntityRef(id) => refs.push(*id),
                _ => {}
            },
            _ => {}
        }
    }
    refs
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_entity(id: EntityId, type_name: &str, params: Vec<Parameter>) -> EntityInstance {
        EntityInstance {
            id,
            type_name: type_name.into(),
            params,
            complex_parts: vec![],
        }
    }

    #[test]
    fn graph_lookup_by_id() {
        let mut log = DiagnosticLog::new();
        let entities = vec![
            make_entity(1, "CARTESIAN_POINT", vec![Parameter::String("a".into())]),
            make_entity(2, "DIRECTION", vec![Parameter::String("b".into())]),
        ];
        let graph = EntityGraph::build(entities, &mut log);
        assert!(!log.has_errors());
        assert_eq!(
            graph.get(1).map(|e| e.type_name.as_str()),
            Some("CARTESIAN_POINT")
        );
        assert_eq!(
            graph.get(2).map(|e| e.type_name.as_str()),
            Some("DIRECTION")
        );
        assert!(graph.get(99).is_none());
    }

    #[test]
    fn graph_detects_broken_refs() {
        let mut log = DiagnosticLog::new();
        let entities = vec![make_entity(
            1,
            "AXIS2_PLACEMENT_3D",
            vec![Parameter::EntityRef(99)],
        )];
        let graph = EntityGraph::build(entities, &mut log);
        assert!(log.has_errors());
        assert_eq!(log.errors().count(), 1);
        assert_eq!(graph.len(), 1);
    }

    #[test]
    fn graph_duplicate_id_logs_warning() {
        let mut log = DiagnosticLog::new();
        let entities = vec![
            make_entity(1, "CARTESIAN_POINT", vec![]),
            make_entity(1, "DIRECTION", vec![]),
        ];
        let graph = EntityGraph::build(entities, &mut log);
        assert_eq!(log.warnings().count(), 1);
        assert_eq!(graph.len(), 1);
        // Second one overwrites
        assert_eq!(
            graph.get(1).map(|e| e.type_name.as_str()),
            Some("DIRECTION")
        );
    }

    #[test]
    fn graph_len() {
        let mut log = DiagnosticLog::new();
        let entities = vec![
            make_entity(1, "A", vec![]),
            make_entity(2, "B", vec![]),
            make_entity(3, "C", vec![]),
        ];
        let graph = EntityGraph::build(entities, &mut log);
        assert_eq!(graph.len(), 3);
        assert!(!graph.is_empty());
    }
}
