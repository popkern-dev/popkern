use super::*;
use crate::naming::RoleCategory;

fn name(op: u64, category: RoleCategory, index: u32, ancestors: Vec<NameId>) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole { category, index },
        ancestors,
    }
}

/// The same lineage minted into two pools in a different order (an
/// upstream insert shifts every id): equal values.
#[test]
fn values_do_not_depend_on_pool_ids() {
    let mut a = NamePool::new();
    let cap = a.intern(name(4, RoleCategory::Cap, 0, vec![]), None);
    let rim_a = a.intern(name(5, RoleCategory::Modified, 14, vec![cap]), None);

    let mut b = NamePool::new();
    for k in 0..7 {
        b.intern(name(7, RoleCategory::Generated, k, vec![]), None);
    }
    let cap = b.intern(name(4, RoleCategory::Cap, 0, vec![]), None);
    let rim_b = b.intern(name(5, RoleCategory::Modified, 14, vec![cap]), None);
    assert_ne!(rim_a, rim_b, "the ids moved");
    assert_eq!(a.value_of(rim_a), b.value_of(rim_b));
}

#[test]
fn relaxations_compare_what_they_say() {
    let mut pool = NamePool::new();
    let cap0 = pool.intern(name(4, RoleCategory::Cap, 0, vec![]), None);
    let cap1 = pool.intern(name(4, RoleCategory::Cap, 1, vec![]), None);
    let recorded_id = pool.intern(name(5, RoleCategory::Modified, 14, vec![cap0]), None);
    let recorded = pool.value_of(recorded_id).expect("value");
    let renumbered_id = pool.intern(name(5, RoleCategory::Modified, 22, vec![cap0]), None);
    let renumbered = pool.value_of(renumbered_id).expect("value");
    let other_cap_id = pool.intern(name(5, RoleCategory::Modified, 22, vec![cap1]), None);
    let other_cap = pool.value_of(other_cap_id).expect("value");
    assert!(!same(&renumbered, &recorded, true, true));
    assert!(
        same(&renumbered, &recorded, false, true),
        "own index ignored"
    );
    assert!(
        !same(&other_cap, &recorded, false, true),
        "the ancestor differs"
    );
    assert!(
        same(&other_cap, &recorded, false, false),
        "only when every index is ignored"
    );
}

#[test]
fn depth_is_bounded() {
    let mut pool = NamePool::new();
    let mut id = pool.intern(name(1, RoleCategory::Lateral, 0, vec![]), None);
    for op in 2..10 {
        id = pool.intern(name(op, RoleCategory::Modified, 0, vec![id]), None);
    }
    let mut value = pool.value_of(id).expect("value");
    let mut depth = 0;
    while let Some(next) = value.ancestors.first().cloned() {
        value = next;
        depth += 1;
    }
    assert_eq!(depth, PROVENANCE_DEPTH);
}

/// A name imported from another pool keeps its ancestry by value: its
/// ancestor ids are the target pool's (scheme 2 booleans, PK-677).
#[test]
fn import_remaps_ancestors_into_the_target_pool() {
    let mut tool = NamePool::new();
    let circle = tool.intern(name(3, RoleCategory::Lateral, 0, vec![]), None);
    let cap_edge = tool.intern(name(4, RoleCategory::Cap, 0, vec![circle]), None);

    let mut result = NamePool::new();
    for k in 0..5 {
        result.intern(name(1, RoleCategory::Junction, k, vec![]), None);
    }
    let copied = result.import(&tool, cap_edge).expect("import");
    assert_eq!(result.value_of(copied), tool.value_of(cap_edge));
    let ancestor = result.get(copied).expect("copy").ancestors[0];
    assert_eq!(
        result.get(ancestor),
        tool.get(circle),
        "the ancestor id is the result pool's"
    );
    assert_eq!(
        result.import(&tool, cap_edge),
        Some(copied),
        "an identical name is the existing entry"
    );
}
