//! Wire edges and acorns are held by shells (PK-1102 S3c): every wire
//! edge by exactly one shell, every held entry a live wire edge, every
//! acorn a live vertex nothing else uses.

use std::collections::{HashMap, HashSet};

use crate::topology::{BRep, EdgeHandle};

use super::ValidationError;

impl BRep {
    /// The edges some shell holds as wire edges.
    pub(super) fn held_wire_edges(&self) -> HashSet<EdgeHandle> {
        self.shells
            .values()
            .flat_map(|s| s.wire_edges.iter().copied())
            .collect()
    }

    /// Every wire edge held by exactly one shell, every held edge a
    /// live wire edge, every acorn's vertex live and reached by no edge.
    pub(super) fn validate_wire_shells(&self) -> Result<(), ValidationError> {
        let mut holder: HashMap<EdgeHandle, usize> = HashMap::new();
        for (sh, shell) in &self.shells {
            for &e in &shell.wire_edges {
                if !self.is_wire_edge(e) {
                    return Err(ValidationError::DanglingReference {
                        what: format!("shell {sh:?} holds {e:?}, which is not a wire edge"),
                    });
                }
                *holder.entry(e).or_default() += 1;
            }
            if let Some(v) = shell.acorn {
                if !self.vertices.contains_key(v) {
                    return Err(ValidationError::DanglingReference {
                        what: format!("acorn shell {sh:?} holds missing vertex {v:?}"),
                    });
                }
                if self.half_edges.values().any(|f| f.origin == v) {
                    return Err(ValidationError::DanglingReference {
                        what: format!("acorn shell {sh:?} holds {v:?}, which an edge reaches"),
                    });
                }
            }
        }
        for e in self.wire_edges() {
            match holder.get(&e).copied().unwrap_or(0) {
                1 => {}
                n => {
                    return Err(ValidationError::DanglingReference {
                        what: format!("wire edge {e:?} is held by {n} shells, not one"),
                    });
                }
            }
        }
        Ok(())
    }
}
