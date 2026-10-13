//! The document after one timeline feature: every alive body.
//!
//! A [`DocumentState`] is what the history scrubber shows at a step and
//! what a replay resumes from. Bodies are shared ([`Arc`]) between
//! neighbouring states, so a step that changes one body of a
//! multi-body document stores one new body, not a copy of the whole
//! document; [`DocumentState::changed_by`] names the step that last
//! produced each body, which is also a stable key for a viewer's mesh
//! of it.

use std::collections::{BTreeMap, HashSet};
use std::sync::Arc;

use popkern_core::topology::BRep;

use crate::node::{BodyId, FeatureId};

/// Every alive body after one timeline step (see the module docs).
#[derive(Debug, Clone, Default)]
pub struct DocumentState {
    bodies: BTreeMap<BodyId, Arc<BRep>>,
    changed_by: BTreeMap<BodyId, FeatureId>,
}

impl DocumentState {
    /// The empty document (before the first feature).
    pub fn empty() -> DocumentState {
        DocumentState::default()
    }

    /// The BRep of `body`, when it is alive in this state.
    pub fn body(&self, body: BodyId) -> Option<&BRep> {
        self.bodies.get(&body).map(Arc::as_ref)
    }

    /// The alive bodies in id order, as the shared handles this state
    /// holds.
    pub(crate) fn shared_bodies(&self) -> impl Iterator<Item = (BodyId, &Arc<BRep>)> {
        self.bodies.iter().map(|(id, brep)| (*id, brep))
    }

    /// The alive bodies in id order.
    pub fn bodies(&self) -> impl Iterator<Item = (BodyId, &BRep)> {
        self.bodies.iter().map(|(id, brep)| (*id, brep.as_ref()))
    }

    /// The alive body ids in id order.
    pub fn body_ids(&self) -> Vec<BodyId> {
        self.bodies.keys().copied().collect()
    }

    /// The number of alive bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// The document's active body: [`BodyId`] 0 when alive, otherwise
    /// the lowest alive id (the rule [`crate::FeatureTree::current_brep`]
    /// follows; never hash order).
    pub fn active_body(&self) -> Option<&BRep> {
        self.bodies
            .get(&BodyId(0))
            .or_else(|| self.bodies.values().next())
            .map(Arc::as_ref)
    }

    /// The timeline feature whose step last produced `body` (built it,
    /// modified it, or merged into it). Two states whose `changed_by`
    /// agree for a body hold the identical BRep for it.
    pub fn changed_by(&self, body: BodyId) -> Option<FeatureId> {
        self.changed_by.get(&body).copied()
    }

    /// Whether `body` is the same shared BRep in `self` and `other`
    /// (unchanged between the two states).
    pub fn shares_body(&self, other: &DocumentState, body: BodyId) -> bool {
        match (self.bodies.get(&body), other.bodies.get(&body)) {
            (Some(a), Some(b)) => Arc::ptr_eq(a, b),
            _ => false,
        }
    }

    /// This state with `body` replaced by `brep`, attributed to `step`.
    pub(crate) fn with_body(&self, body: BodyId, brep: BRep, step: FeatureId) -> DocumentState {
        let mut next = self.clone();
        next.bodies.insert(body, Arc::new(brep));
        next.changed_by.insert(body, step);
        next
    }

    /// The bodies as a working map sharing this state's handles; a
    /// step copies a body only when it takes it out to change it, so
    /// the cached state stays immutable (PK-812).
    pub(crate) fn to_working(&self) -> crate::working_bodies::WorkingBodies {
        crate::working_bodies::WorkingBodies::shared(
            self.bodies
                .iter()
                .map(|(id, brep)| (*id, brep.clone()))
                .collect(),
        )
    }

    /// The state after `step`: bodies in `touched` (or new to the map)
    /// take fresh copies from `working`; every other alive body shares
    /// `self`'s BRep. Bodies absent from `working` are gone (consumed by
    /// a boolean or merge).
    pub(crate) fn advanced(
        &self,
        working: &crate::working_bodies::WorkingBodies,
        touched: Touched<'_>,
        step: FeatureId,
    ) -> DocumentState {
        let mut next = DocumentState::default();
        for (&id, brep) in working.as_map() {
            let shared = match touched {
                Touched::All => None,
                Touched::Bodies(set) if set.contains(&id) => None,
                Touched::Bodies(_) => self.bodies.get(&id).cloned(),
            };
            match shared {
                Some(arc) => {
                    next.bodies.insert(id, arc);
                    if let Some(origin) = self.changed_by.get(&id) {
                        next.changed_by.insert(id, *origin);
                    }
                }
                None => {
                    // Shared with the working map: the next step that
                    // changes this body copies it (`WorkingBodies::remove`).
                    next.bodies.insert(id, brep.clone());
                    next.changed_by.insert(id, step);
                }
            }
        }
        next
    }

    /// Estimated bytes of the bodies this state holds that no state in
    /// `seen` already holds (shared bodies count once across a cache).
    pub(crate) fn unique_bytes(&self, seen: &mut HashSet<*const BRep>) -> usize {
        self.bodies
            .values()
            .filter(|brep| seen.insert(Arc::as_ptr(brep)))
            .map(|brep| estimate_brep_bytes(brep))
            .sum()
    }
}

/// Which bodies one step may have changed.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Touched<'a> {
    /// Every body (pattern features, whose bodies are not tracked
    /// individually).
    All,
    /// Exactly these bodies (plus any body new to the working map).
    Bodies(&'a HashSet<BodyId>),
}

/// Estimate the heap memory of a single BRep, in bytes: vertices × 64 +
/// edges × 96 + faces × 256 (a budget heuristic, not an exact count).
pub(crate) fn estimate_brep_bytes(brep: &BRep) -> usize {
    brep.vertices.len() * 64 + brep.edges.len() * 96 + brep.faces.len() * 256
}

#[cfg(test)]
#[path = "document_state_tests.rs"]
mod tests;
