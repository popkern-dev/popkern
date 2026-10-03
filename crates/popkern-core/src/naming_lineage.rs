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

use crate::naming::{EntityName, NameId, NamePool};

impl NamePool {
    /// The recorded lineage of `id` (empty when none).
    pub fn lineage(&self, id: NameId) -> &[EntityName] {
        self.lineage
            .get(id.0 as usize)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    /// Record that `id` descends from `from` (deduplicated).
    pub fn add_lineage(&mut self, id: NameId, from: impl IntoIterator<Item = EntityName>) {
        let i = id.0 as usize;
        if i >= self.names.len() {
            return;
        }
        if self.lineage.len() < self.names.len() {
            self.lineage.resize(self.names.len(), Vec::new());
        }
        for name in from {
            if self.names[i] != name && !self.lineage[i].contains(&name) {
                self.lineage[i].push(name);
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

    /// Whether `id` is `target` or descends from it (ancestors or
    /// lineage, by value).
    pub fn descends(&self, id: NameId, target: &EntityName) -> bool {
        self.get(id) == Some(target) || self.ancestry(id).contains(target)
    }
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
            std::iter::once(first.names[piece.0 as usize].clone()).chain(first.ancestry(piece)),
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
