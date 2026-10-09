//! Built part bodies kept across interference checks (PK-723).
//!
//! An interference check needs every definition's bodies and their
//! content keys (`mesh_cache`). Building them meant cloning each body
//! and validating it, and the content key formats the whole body: about
//! 0.5 s per check on the six-axis arm, paid at every step of a swing
//! whose definitions never change.
//!
//! An entry is dropped whenever its definition's model is borrowed
//! mutably (`Assembly::part_model_mut` / `part_models_mut`: drives,
//! edits, appearance), and it carries a signature — the definition's
//! name and parameter values — so an entry can never answer for a
//! different definition at a reused index.

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use popkern_part::Part;

use crate::assembly::{Assembly, DefinitionId};
use crate::error::AssemblyError;

/// One definition's built bodies and their content keys.
#[derive(Debug)]
pub(crate) struct CachedParts {
    signature: u64,
    pub(crate) parts: Vec<Part>,
    pub(crate) contents: Vec<u64>,
}

/// Built bodies by definition.
#[derive(Debug, Default)]
pub(crate) struct PartsCache {
    entries: Mutex<HashMap<DefinitionId, Arc<CachedParts>>>,
}

impl Clone for PartsCache {
    fn clone(&self) -> Self {
        PartsCache {
            entries: Mutex::new(self.entries.lock().map(|m| m.clone()).unwrap_or_default()),
        }
    }
}

impl PartsCache {
    /// Drop `id`'s entry (its model may change).
    pub(crate) fn forget(&self, id: DefinitionId) {
        if let Ok(mut entries) = self.entries.lock() {
            entries.remove(&id);
        }
    }

    /// Drop every entry.
    pub(crate) fn clear(&self) {
        if let Ok(mut entries) = self.entries.lock() {
            entries.clear();
        }
    }
}

impl Assembly {
    /// The built bodies of every part definition in `ids`, from the
    /// cache when its signature still holds.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] when a definition fails to build.
    pub(crate) fn cached_parts(
        &mut self,
        ids: &std::collections::HashSet<DefinitionId>,
    ) -> Result<HashMap<DefinitionId, Arc<CachedParts>>, AssemblyError> {
        let mut out = HashMap::new();
        let mut sorted: Vec<DefinitionId> = ids.iter().copied().collect();
        sorted.sort();
        for id in sorted {
            let Some(definition) = self.definitions.get(id.0 as usize) else {
                continue;
            };
            let Some(model) = definition.model() else {
                continue;
            };
            let signature = {
                use std::hash::{Hash, Hasher};
                let mut h = std::hash::DefaultHasher::new();
                definition.name.hash(&mut h);
                format!("{:?}", model.params()).hash(&mut h);
                h.finish()
            };
            let cached = self
                .parts_cache
                .entries
                .lock()
                .ok()
                .and_then(|e| e.get(&id).cloned())
                .filter(|c| c.signature == signature);
            let entry = match cached {
                Some(c) => c,
                None => {
                    let name = definition.name.clone();
                    let parts = match &mut self.definitions[id.0 as usize].content {
                        crate::assembly::DefinitionContent::Part(model) => model
                            .to_parts_all()
                            .map_err(|source| AssemblyError::Script {
                                definition: name,
                                source,
                            })?,
                        crate::assembly::DefinitionContent::Assembly(_) => continue,
                    };
                    let contents = parts
                        .iter()
                        .map(|p| crate::mesh_cache::fingerprint(p.brep()))
                        .collect();
                    let entry = Arc::new(CachedParts {
                        signature,
                        parts,
                        contents,
                    });
                    if let Ok(mut e) = self.parts_cache.entries.lock() {
                        e.insert(id, entry.clone());
                    }
                    entry
                }
            };
            out.insert(id, entry);
        }
        Ok(out)
    }
}
