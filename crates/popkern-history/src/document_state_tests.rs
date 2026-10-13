// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use std::collections::{HashMap, HashSet};

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use slotmap::SlotMap;

use super::{DocumentState, Touched};
use crate::node::{BodyId, FeatureId};
use crate::working_bodies::WorkingBodies;

fn cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("unit cube")
}

/// A working map over fresh bodies.
fn working(bodies: Vec<(BodyId, BRep)>) -> WorkingBodies {
    WorkingBodies::shared(
        bodies
            .into_iter()
            .map(|(id, brep)| (id, std::sync::Arc::new(brep)))
            .collect::<HashMap<_, _>>(),
    )
}

fn features(n: usize) -> Vec<FeatureId> {
    let mut map: SlotMap<FeatureId, ()> = SlotMap::with_key();
    (0..n).map(|_| map.insert(())).collect()
}

#[test]
fn untouched_bodies_are_shared_and_keep_their_origin() {
    let f = features(2);
    let mut working = working(vec![(BodyId(0), cube()), (BodyId(1), cube())]);
    let first = DocumentState::empty().advanced(&working, Touched::All, f[0]);
    // A step changes body 1: it takes it out (a copy, since `first`
    // shares it) and puts the result back.
    let body = working.remove(&BodyId(1)).expect("alive");
    working.insert(BodyId(1), body);
    let touched = HashSet::from([BodyId(1)]);
    let second = first.advanced(&working, Touched::Bodies(&touched), f[1]);

    assert!(second.shares_body(&first, BodyId(0)), "body 0 untouched");
    assert!(!second.shares_body(&first, BodyId(1)), "body 1 copied");
    assert_eq!(second.changed_by(BodyId(0)), Some(f[0]));
    assert_eq!(second.changed_by(BodyId(1)), Some(f[1]));
    // Shared bodies count once across states.
    let mut seen = HashSet::new();
    let total = first.unique_bytes(&mut seen) + second.unique_bytes(&mut seen);
    let one = super::estimate_brep_bytes(&cube());
    assert_eq!(total, 3 * one, "two bodies + one new copy");
}

#[test]
fn consumed_bodies_leave_and_new_bodies_arrive() {
    let f = features(2);
    let both = working(vec![(BodyId(0), cube()), (BodyId(1), cube())]);
    let first = DocumentState::empty().advanced(&both, Touched::All, f[0]);
    // A boolean consumed body 1; nothing reported touched, body 2 is new.
    let mut after = both.clone();
    let _ = after.remove(&BodyId(1));
    after.insert(BodyId(2), cube());
    let none = HashSet::new();
    let second = first.advanced(&after, Touched::Bodies(&none), f[1]);
    assert_eq!(second.body_ids(), vec![BodyId(0), BodyId(2)]);
    assert!(second.shares_body(&first, BodyId(0)));
    assert_eq!(second.changed_by(BodyId(2)), Some(f[1]));
}

#[test]
fn active_body_prefers_body_zero_then_the_lowest_id() {
    let f = features(1);
    let low = working(vec![(BodyId(3), cube()), (BodyId(5), cube())]);
    let state = DocumentState::empty().advanced(&low, Touched::All, f[0]);
    assert!(state.active_body().is_some());
    assert_eq!(state.body_ids()[0], BodyId(3));
    assert!(DocumentState::empty().active_body().is_none());
}
