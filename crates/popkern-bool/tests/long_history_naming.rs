//! Naming cost over a long boolean history (PK-731).
//!
//! A disc cut by a row of slots: every boolean renames every result
//! entity, so a body's name pool carries its whole history. Three rules
//! keep a late boolean from paying for that history:
//!
//! - naming scheme 3 starts the result pool as the target's, so the
//!   target's names keep their ids (scheme 2 re-imported, and
//!   renumbered, every one of them per boolean);
//! - lineage is kept only where a reader can reach it — live names and
//!   their direct ancestors — and pruning it never changes what any
//!   live name descends from;
//! - the pool is shared copy-on-write, so the result's pool extends the
//!   target's by pointer.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, NameId, OperationId, RoleCategory};
use popkern_core::naming_lineage::prune_unreadable_lineage;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn slot(k: usize, n: usize, tol: Tolerance) -> BRep {
    let a = std::f64::consts::TAU * k as f64 / n as f64;
    let (ca, sa, w) = (a.cos(), a.sin(), 0.6);
    let pts = [
        Point3::new(22.0 * ca + w * sa, 22.0 * sa - w * ca, -1.0),
        Point3::new(27.0 * ca + w * sa, 27.0 * sa - w * ca, -1.0),
        Point3::new(27.0 * ca - w * sa, 27.0 * sa + w * ca, -1.0),
        Point3::new(22.0 * ca - w * sa, 22.0 * sa + w * ca, -1.0),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 4.0, tol).expect("slot prism")
}

/// Every named entity's name id.
fn live_names(brep: &BRep) -> Vec<NameId> {
    brep.vertices
        .values()
        .filter_map(|v| v.name)
        .chain(brep.edges.values().filter_map(|e| e.name))
        .chain(brep.faces.values().filter_map(|f| f.name))
        .collect()
}

/// A disc with `n` slots cut under naming scheme `scheme`, and the
/// body after every cut.
fn slotted(n: usize, scheme: u32) -> Vec<BRep> {
    let tol = Tolerance::default();
    let mut bodies = vec![BRep::cylinder(30.0, 1.5, tol).expect("disc")];
    for k in 0..n {
        let mut ctx = TrackingContext::new(OperationId(100 + k as u64)).with_naming_scheme(scheme);
        let target = bodies.last().cloned().expect("a body");
        let cut = boolean(
            BoolOperation::Subtraction,
            target,
            slot(k, n, tol),
            tol,
            Some(&mut ctx),
        )
        .expect("slot cut serves");
        bodies.push(cut);
    }
    bodies
}

#[test]
fn scheme_three_results_extend_the_target_pool() {
    let bodies = slotted(8, 3);
    for pair in bodies[1..].windows(2) {
        assert!(
            pair[1].name_pool.keeps_names_of(&pair[0].name_pool),
            "a boolean result keeps every target name at its id"
        );
    }
}

#[test]
fn scheme_two_keeps_the_import_rule() {
    // A scheme-2 document replays with the rules it was recorded under:
    // the result pool is built by import, which renumbers the target's
    // history, so it does not extend the target pool.
    let bodies = slotted(3, 2);
    assert!(!bodies[3].name_pool.keeps_names_of(&bodies[2].name_pool));
}

#[test]
fn lineage_lives_only_where_a_reader_reaches_it() {
    let bodies = slotted(12, 3);
    let last = bodies.last().expect("a body");
    let mut readable = std::collections::HashSet::new();
    for id in live_names(last) {
        readable.insert(id);
        if let Some(name) = last.name_pool.get(id) {
            readable.extend(name.ancestors.iter().copied());
        }
    }
    let stray: Vec<NameId> = last
        .name_pool
        .names()
        .map(|(id, _)| id)
        .filter(|id| !last.name_pool.lineage(*id).is_empty() && !readable.contains(id))
        .collect();
    assert!(stray.is_empty(), "unreadable lineage kept on {stray:?}");
}

#[test]
fn pruning_never_changes_what_a_live_name_descends_from() {
    let bodies = slotted(6, 3);
    let mut body = bodies.last().cloned().expect("a body");
    // Lineage on every name no reader reaches: what a pool held before
    // pruning existed.
    let marker = EntityName {
        operation: OperationId(9_999),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: Vec::new(),
    };
    let ids: Vec<NameId> = body.name_pool.names().map(|(id, _)| id).collect();
    for id in &ids {
        body.name_pool.add_lineage(*id, [marker.clone()]);
    }
    let live = live_names(&body);
    let before: Vec<Vec<EntityName>> = live.iter().map(|id| body.name_pool.ancestry(*id)).collect();
    prune_unreadable_lineage(&mut body);
    let after: Vec<Vec<EntityName>> = live.iter().map(|id| body.name_pool.ancestry(*id)).collect();
    assert_eq!(before, after, "every live name's ancestry is unchanged");
    let mut readable = std::collections::HashSet::new();
    for id in &live {
        readable.insert(*id);
        if let Some(name) = body.name_pool.get(*id) {
            readable.extend(name.ancestors.iter().copied());
        }
    }
    for id in &ids {
        let marked = body.name_pool.lineage(*id).contains(&marker);
        assert_eq!(
            marked,
            readable.contains(id),
            "{id:?}: kept iff a reader reaches it"
        );
    }
}
