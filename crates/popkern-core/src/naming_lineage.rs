//! Name lineage (PK-633): what each named entity was cut, trimmed or
//! carried from, kept beside its stable name.
//!
//! A boolean keeps the stable names of the faces it merely moves
//! (`Modified`, with the input as ancestor) but mints `Generated`
//! names for faces whose shape it changed past the centroid match. The
//! lineage records, for every output entity, the input entities it is a
//! piece of — found geometrically (same surface, inside the input's
//! trim) — as name VALUES with their own ancestry flattened in, so a
//! label bound to an input entity resolves to its pieces after any
//! number of later operations, without changing any stable name (no
//! naming-scheme change; saved references keep resolving).

use std::collections::HashSet;
use std::sync::Arc;

use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::naming::{EntityName, NameId, NamePool};
use crate::topology::BRep;

/// One name's recorded lineage, shared by reference: the pool's chunks
/// are copy-on-write, and a chunk copy (any write into a chunk a cloned
/// pool still shares) deep-copied 256 lineage lists — every boolean's
/// prune and lineage record paid it (PK-803). Serialized as the plain
/// list it holds.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct LineageEntry(Arc<Vec<EntityName>>);

impl LineageEntry {
    fn as_slice(&self) -> &[EntityName] {
        self.0.as_slice()
    }

    fn push(&mut self, name: EntityName) {
        Arc::make_mut(&mut self.0).push(name);
    }
}

impl Serialize for LineageEntry {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.0.as_ref().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for LineageEntry {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        Vec::<EntityName>::deserialize(deserializer).map(|v| Self(Arc::new(v)))
    }
}

impl NamePool {
    /// The recorded lineage of `id` (empty when none).
    pub fn lineage(&self, id: NameId) -> &[EntityName] {
        self.lineage
            .get(id.0 as usize)
            .map(LineageEntry::as_slice)
            .unwrap_or(&[])
    }

    /// Record that `id` descends from `from` (deduplicated).
    pub fn add_lineage(&mut self, id: NameId, from: impl IntoIterator<Item = EntityName>) {
        let i = id.0 as usize;
        let Some(own) = self.get(id).cloned() else {
            return;
        };
        let fresh: Vec<EntityName> = from
            .into_iter()
            .filter(|name| *name != own && !self.lineage(id).contains(name))
            .collect();
        if fresh.is_empty() {
            return;
        }
        if let Some(entry) = self.lineage.get_mut(i) {
            for name in fresh {
                if !entry.as_slice().contains(&name) {
                    entry.push(name);
                }
            }
        }
    }

    /// Everything `id` descends from, as values: its recorded lineage
    /// and its DIRECT ancestors. Ancestor ids are only meaningful one
    /// level deep — an ancestor's value was interned here from another
    /// operation's pool, so ITS ancestor ids index that pool, not this
    /// one — which is why deeper history travels as lineage values
    /// (each operation records its pieces' full history, read soundly
    /// in the pool it came from) rather than by walking ids.
    pub fn ancestry(&self, id: NameId) -> Vec<EntityName> {
        let mut out: Vec<EntityName> = self.lineage(id).to_vec();
        if let Some(name) = self.get(id) {
            for ancestor in &name.ancestors {
                if let Some(value) = self.get(*ancestor)
                    && !out.contains(value)
                {
                    out.push(value.clone());
                }
                for value in self.lineage(*ancestor) {
                    if !out.contains(value) {
                        out.push(value.clone());
                    }
                }
            }
        }
        out
    }

    /// Drop the recorded lineage of every name outside `keep` (the
    /// names themselves stay: ids and stable names never change).
    pub fn retain_lineage(&mut self, keep: &HashSet<NameId>) {
        self.lineage.edit_where(
            |i, entry| !entry.as_slice().is_empty() && !keep.contains(&NameId(i as u32)),
            |_, entry| *entry = LineageEntry::default(),
        );
    }

    /// Whether `id` is `target` or descends from it (ancestors or
    /// lineage, by value).
    pub fn descends(&self, id: NameId, target: &EntityName) -> bool {
        self.get(id) == Some(target) || self.ancestry(id).contains(target)
    }
}

/// Drop the lineage no reader of `brep` can reach (PK-731).
///
/// Lineage is read only through [`NamePool::ancestry`], which reads the
/// lineage of the queried name and of its DIRECT ancestors; the queried
/// names are live entities' (`resolve_lineage`, and the lineage a later
/// boolean records reads its live inputs'), and every name a later
/// operation mints takes live input entities as ancestors. So the
/// lineage of any name that is neither a live entity's nor a live
/// entity's direct ancestor is never read again. Each boolean records
/// a piece's whole ancestry by value, so keeping it on every historical
/// name made a pool's lineage grow with the cube of the operation
/// count — and every BRep clone copied it. Resolution is unchanged by
/// construction: every value any [`NamePool::ancestry`] of a live name
/// returns is kept.
pub fn prune_unreadable_lineage(brep: &mut BRep) {
    let live = brep
        .vertices
        .values()
        .filter_map(|v| v.name)
        .chain(brep.edges.values().filter_map(|e| e.name))
        .chain(brep.faces.values().filter_map(|f| f.name));
    let mut keep: HashSet<NameId> = HashSet::new();
    for id in live {
        keep.insert(id);
        if let Some(name) = brep.name_pool.get(id) {
            keep.extend(name.ancestors.iter().copied());
        }
    }
    brep.name_pool.retain_lineage(&keep);
}

#[cfg(test)]
mod tests {
    use crate::naming::{EntityRole, OperationId, RoleCategory};

    use super::*;

    fn name(op: u64, cat: RoleCategory, idx: u32) -> EntityName {
        EntityName {
            operation: OperationId(op),
            role: EntityRole {
                category: cat,
                index: idx,
            },
            ancestors: Vec::new(),
        }
    }

    #[test]
    fn lineage_survives_merge_and_compact_and_is_transitive() {
        let mut first = NamePool::new();
        let original = name(1, RoleCategory::Lateral, 0);
        let piece = first.intern(name(2, RoleCategory::Generated, 0), None);
        first.add_lineage(piece, [original.clone()]);
        assert!(first.descends(piece, &original));
        // A later operation records the piece's whole ancestry.
        let mut second = NamePool::new();
        let again = second.intern(name(3, RoleCategory::Generated, 0), None);
        second.add_lineage(
            again,
            first
                .get(piece)
                .cloned()
                .into_iter()
                .chain(first.ancestry(piece)),
        );
        assert!(second.descends(again, &original), "two operations later");
        let mut merged = NamePool::new();
        merged.intern(name(9, RoleCategory::Cap, 0), None);
        let remap = merged.merge(&second);
        assert!(merged.descends(remap[&again], &original));
        let keep: std::collections::HashSet<_> = [remap[&again]].into_iter().collect();
        let compacted = merged.compact(&keep);
        assert!(merged.descends(compacted[&remap[&again]], &original));
    }
}
