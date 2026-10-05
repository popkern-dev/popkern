//! Name provenance by value (PK-677): a persistent name spelled
//! without pool ids, so a reference survives edits upstream of it.
//!
//! An [`EntityName`] holds its ancestors as [`NameId`]s, indices into
//! the pool of the body it was minted in. A replay that runs the same
//! features rebuilds the same pool, so recorded names resolve. A
//! feature inserted upstream (history-playhead inserts, PK-666) adds
//! names first: every later pool id shifts, and the role indices an
//! operation hands out renumber too. Measured on the `labelled_rim`
//! fixture, a chamfer before the bore subtraction moved the rim edge
//! from `boolean/Modified/14` with ancestor id 42 to `Modified/22` with
//! ancestor id 66, while the ancestor itself (the bore extrude's cap)
//! did not change. An exact lookup then either misses or, worse,
//! lands on an unrelated entity that happens to carry the same numbers.
//!
//! A [`NameValue`] spells a name by value: its operation and role, and
//! its ancestors spelled the same way, to [`PROVENANCE_DEPTH`] levels.
//! [`resolve_by_provenance`] finds the entity a snapshot describes,
//! relaxing step by step (the exact value, then ignoring the entity's
//! own index, then ignoring every ancestor's index). It binds only
//! when a step leaves exactly one candidate, and reports an ambiguity
//! instead of picking one by index.

use serde::{Deserialize, Serialize};

use crate::naming::{EntityName, EntityRole, NameId, NamePool, OperationId};
use crate::naming_resolve::ResolvedEntity;
use crate::topology::BRep;

/// How many ancestor levels a snapshot keeps.
pub const PROVENANCE_DEPTH: usize = 3;

/// A bound on ancestry depth when importing (ancestors are older names,
/// so a real chain is at most the number of operations; this only
/// stops a corrupt, cyclic pool).
const MAX_IMPORT_DEPTH: usize = 4096;

/// A name spelled by value: no pool ids anywhere.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct NameValue {
    /// The operation that produced the entity.
    pub operation: OperationId,
    /// Its role in that operation.
    pub role: EntityRole,
    /// Its ancestors, by value.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ancestors: Vec<NameValue>,
}

/// The kind of entity a reference names.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EntityKind {
    /// A vertex.
    Vertex,
    /// An edge.
    Edge,
    /// A face.
    Face,
}

/// How far [`resolve_by_provenance`] had to relax to find its match.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Relaxation {
    /// The whole value matched (only the pool ids had moved).
    Exact,
    /// The entity's own role index was renumbered.
    OwnIndex,
    /// Role indices were renumbered in the ancestry as well.
    AncestorIndices,
}

/// The outcome of resolving a snapshot.
#[derive(Debug, Clone, PartialEq)]
pub enum ProvenanceMatch {
    /// One entity matches: it and its current name.
    Unique {
        /// The entity.
        entity: ResolvedEntity,
        /// Its name in the BRep's pool now.
        name: EntityName,
        /// How the match was reached.
        relaxation: Relaxation,
    },
    /// Several entities match at the first relaxation that matches any.
    Ambiguous {
        /// How many.
        candidates: usize,
        /// At which relaxation.
        relaxation: Relaxation,
    },
    /// Nothing matches, even with every index ignored.
    NotFound,
}

impl NamePool {
    /// Intern name `id` of `source` into this pool by value: its
    /// ancestors first, so every ancestor id the copy holds is this
    /// pool's; criteria and lineage come along. The copy of an identical
    /// name is the existing entry. `None` for an unknown id (or an
    /// ancestry deeper than any operation chain, which only a corrupt
    /// pool has).
    pub fn import(&mut self, source: &NamePool, id: NameId) -> Option<NameId> {
        let mut memo = std::collections::HashMap::new();
        self.import_memo(source, id, &mut memo, 0)
    }

    fn import_memo(
        &mut self,
        source: &NamePool,
        id: NameId,
        memo: &mut std::collections::HashMap<NameId, NameId>,
        depth: usize,
    ) -> Option<NameId> {
        if let Some(&copied) = memo.get(&id) {
            return Some(copied);
        }
        if depth > MAX_IMPORT_DEPTH {
            return None;
        }
        let name = source.get(id)?.clone();
        let ancestors = name
            .ancestors
            .iter()
            .filter_map(|a| self.import_memo(source, *a, memo, depth + 1))
            .collect();
        let copy = EntityName {
            operation: name.operation,
            role: name.role,
            ancestors,
        };
        let criteria = source.get_criteria(id).cloned();
        let copied = match self.find(&copy) {
            Some(existing) => existing,
            None => self.intern(copy, criteria),
        };
        self.add_lineage(copied, source.lineage(id).iter().cloned());
        memo.insert(id, copied);
        Some(copied)
    }

    /// The name `id` spelled by value (`None` for an unknown id).
    pub fn value_of(&self, id: NameId) -> Option<NameValue> {
        self.value_at(id, PROVENANCE_DEPTH)
    }

    fn value_at(&self, id: NameId, depth: usize) -> Option<NameValue> {
        let name = self.get(id)?;
        let ancestors = if depth == 0 {
            Vec::new()
        } else {
            name.ancestors
                .iter()
                .filter_map(|a| self.value_at(*a, depth - 1))
                .collect()
        };
        Some(NameValue {
            operation: name.operation,
            role: name.role,
            ancestors,
        })
    }
}

/// Equality of two values under a relaxation: `own` compares the role
/// index at this level, `deep` the ancestors' indices.
fn same(a: &NameValue, b: &NameValue, own: bool, deep: bool) -> bool {
    a.operation == b.operation
        && a.role.category == b.role.category
        && (!own || a.role.index == b.role.index)
        && same_ancestors(&a.ancestors, &b.ancestors, deep)
}

/// Ancestors equal as multisets (order is not part of identity).
fn same_ancestors(a: &[NameValue], b: &[NameValue], deep: bool) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut used = vec![false; b.len()];
    a.iter().all(|x| {
        let found = (0..b.len()).find(|&j| !used[j] && same(x, &b[j], deep, deep));
        if let Some(j) = found {
            used[j] = true;
        }
        found.is_some()
    })
}

/// The named entities of `kind` in `brep`, with their name ids.
fn named(brep: &BRep, kind: EntityKind) -> Vec<(ResolvedEntity, NameId)> {
    match kind {
        EntityKind::Vertex => brep
            .vertices
            .iter()
            .filter_map(|(h, v)| Some((ResolvedEntity::Vertex(h), v.name?)))
            .collect(),
        EntityKind::Edge => brep
            .edges
            .iter()
            .filter_map(|(h, e)| Some((ResolvedEntity::Edge(h), e.name?)))
            .collect(),
        EntityKind::Face => brep
            .faces
            .iter()
            .filter_map(|(h, f)| Some((ResolvedEntity::Face(h), f.name?)))
            .collect(),
    }
}

/// Find the entity of `kind` in `brep` that `snapshot` describes.
pub fn resolve_by_provenance(
    snapshot: &NameValue,
    kind: EntityKind,
    brep: &BRep,
) -> ProvenanceMatch {
    let candidates: Vec<(ResolvedEntity, NameId, NameValue)> = named(brep, kind)
        .into_iter()
        .filter_map(|(entity, id)| Some((entity, id, brep.name_pool.value_of(id)?)))
        .filter(|(_, _, value)| same(value, snapshot, false, false))
        .collect();
    let steps = [
        (Relaxation::Exact, true, true),
        (Relaxation::OwnIndex, false, true),
        (Relaxation::AncestorIndices, false, false),
    ];
    for (relaxation, own, deep) in steps {
        let matches: Vec<&(ResolvedEntity, NameId, NameValue)> = candidates
            .iter()
            .filter(|(_, _, value)| same(value, snapshot, own, deep))
            .collect();
        match matches.as_slice() {
            [] => continue,
            [(entity, id, _)] => {
                return match brep.name_pool.get(*id) {
                    Some(name) => ProvenanceMatch::Unique {
                        entity: *entity,
                        name: name.clone(),
                        relaxation,
                    },
                    None => ProvenanceMatch::NotFound,
                };
            }
            many => {
                return ProvenanceMatch::Ambiguous {
                    candidates: many.len(),
                    relaxation,
                };
            }
        }
    }
    ProvenanceMatch::NotFound
}

/// The kind of a resolved entity.
pub fn kind_of(entity: ResolvedEntity) -> EntityKind {
    match entity {
        ResolvedEntity::Vertex(_) => EntityKind::Vertex,
        ResolvedEntity::Edge(_) => EntityKind::Edge,
        ResolvedEntity::Face(_) => EntityKind::Face,
    }
}

#[cfg(test)]
#[path = "naming_provenance_tests.rs"]
mod tests;
