//! The replay's working body map, copy-on-write (PK-812).
//!
//! A replay resumes from a cached [`crate::document_state::DocumentState`]
//! whose bodies are shared `Arc<BRep>` handles. The working map used to
//! deep-copy every alive body on resume, and the next state deep-copied
//! every touched body again, so each recorded feature paid two full
//! copies of every body in the document (measured on a 64-slot disc
//! script: the 255-face disc copied per sketch, per tool extrude and per
//! boolean). Here the working map shares the state's handles: a step
//! that takes a body out to change it ([`WorkingBodies::remove`]) copies
//! it only when a state still shares it, untouched bodies are never
//! copied, and the next state shares the working handles back.

use std::collections::HashMap;
use std::sync::Arc;

use popkern_core::topology::BRep;

use crate::node::BodyId;

/// Alive bodies during a replay, shared with the document states the
/// replay reads and writes (see the module docs).
#[derive(Debug, Clone, Default)]
pub(crate) struct WorkingBodies {
    bodies: HashMap<BodyId, Arc<BRep>>,
}

impl WorkingBodies {
    /// A working map over `bodies`, sharing every handle.
    pub(crate) fn shared(bodies: HashMap<BodyId, Arc<BRep>>) -> Self {
        Self { bodies }
    }

    /// The BRep of `body`, if alive.
    pub(crate) fn get(&self, body: &BodyId) -> Option<&BRep> {
        self.bodies.get(body).map(Arc::as_ref)
    }

    /// The shared handle of `body`, if alive.
    #[cfg(test)]
    pub(crate) fn get_shared(&self, body: &BodyId) -> Option<&Arc<BRep>> {
        self.bodies.get(body)
    }

    /// Take `body` out to change it. Copies the BRep only when another
    /// holder (a cached state) still shares it.
    pub(crate) fn remove(&mut self, body: &BodyId) -> Option<BRep> {
        self.bodies.remove(body).map(Arc::unwrap_or_clone)
    }

    /// Put `brep` in as `body`'s current BRep.
    pub(crate) fn insert(&mut self, body: BodyId, brep: BRep) {
        self.bodies.insert(body, Arc::new(brep));
    }

    /// The alive body ids (hash order: sort before relying on order).
    pub(crate) fn keys(&self) -> impl Iterator<Item = &BodyId> {
        self.bodies.keys()
    }

    /// The alive bodies (hash order: sort before relying on order).
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&BodyId, &BRep)> {
        self.bodies.iter().map(|(id, brep)| (id, brep.as_ref()))
    }

    /// The underlying shared map (read-only views such as
    /// [`crate::measure_query::MeasureContext`]).
    pub(crate) fn as_map(&self) -> &HashMap<BodyId, Arc<BRep>> {
        &self.bodies
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn untouched_bodies_stay_shared_and_a_taken_body_is_copied_once() {
        let cube = Arc::new(BRep::unit_cube(Tolerance::default()).expect("cube"));
        let other = Arc::new(BRep::unit_cube(Tolerance::default()).expect("cube"));
        let state = HashMap::from([(BodyId(0), cube.clone()), (BodyId(1), other.clone())]);
        let mut working = WorkingBodies::shared(state.clone());
        // Untouched: the working map holds the state's own allocation.
        let shared = working.get_shared(&BodyId(1)).expect("alive");
        assert!(Arc::ptr_eq(shared, &other));
        // Taken out while the state shares it: a copy, the state intact.
        let taken = working.remove(&BodyId(0)).expect("alive");
        assert_eq!(taken.faces.len(), cube.faces.len());
        assert!(Arc::strong_count(&cube) >= 2, "the state keeps its body");
        working.insert(BodyId(0), taken);
        assert!(!Arc::ptr_eq(
            working.get_shared(&BodyId(0)).expect("alive"),
            &cube
        ));
    }
}
