use super::*;

fn make_name(op: u64, cat: RoleCategory, idx: u32) -> EntityName {
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
fn intern_and_get_round_trip() {
    let mut pool = NamePool::new();
    let name = make_name(1, RoleCategory::Generated, 0);
    let id = pool.intern(name.clone(), None);
    let retrieved = pool.get(id).expect("name should be present");
    assert_eq!(*retrieved, name);
}

#[test]
fn dedup_same_name_same_id() {
    let mut pool = NamePool::new();
    let name = make_name(1, RoleCategory::Cap, 0);
    let id1 = pool.intern(name.clone(), None);
    let id2 = pool.intern(name.clone(), None);
    assert_eq!(id1, id2);
    assert_eq!(pool.len(), 1);
}

#[test]
fn different_names_different_ids() {
    let mut pool = NamePool::new();
    let a = make_name(1, RoleCategory::Generated, 0);
    let b = make_name(1, RoleCategory::Generated, 1);
    let id_a = pool.intern(a, None);
    let id_b = pool.intern(b, None);
    assert_ne!(id_a, id_b);
    assert_eq!(pool.len(), 2);
}

#[test]
fn find_existing_returns_some() {
    let mut pool = NamePool::new();
    let name = make_name(42, RoleCategory::Lateral, 3);
    let id = pool.intern(name.clone(), None);
    assert_eq!(pool.find(&name), Some(id));
}

#[test]
fn find_missing_returns_none() {
    let pool = NamePool::new();
    let name = make_name(99, RoleCategory::Modified, 0);
    assert_eq!(pool.find(&name), None);
}

#[test]
fn rebuild_index_restores_lookups() {
    let mut pool = NamePool::new();
    let name = make_name(7, RoleCategory::Intersection, 0);
    let id = pool.intern(name.clone(), None);

    // Simulate what happens after serde: index is gone.
    pool.index = Arc::new(HashMap::new());
    assert_eq!(pool.find(&name), None, "index cleared");

    pool.rebuild_index();
    assert_eq!(pool.find(&name), Some(id), "index restored");
}

#[test]
fn merge_imports_names_and_remaps_ancestors() {
    let mut pool_a = NamePool::new();
    let mut pool_b = NamePool::new();

    let base = make_name(1, RoleCategory::Cap, 0);
    let base_id_b = pool_b.intern(base.clone(), None);

    let derived = EntityName {
        operation: OperationId(2),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: vec![base_id_b],
    };
    pool_b.intern(derived.clone(), None);

    let remap = pool_a.merge(&pool_b);

    // Both names should now be in pool_a.
    assert_eq!(pool_a.len(), 2);

    // The base id in pool_b was remapped.
    let new_base_id = remap[&base_id_b];
    assert_eq!(pool_a.get(new_base_id).unwrap().operation, OperationId(1));

    // Find the derived name via its new ancestor id.
    let found_derived = pool_a
        .names()
        .map(|(_, n)| n)
        .find(|n| n.operation == OperationId(2));
    let found_derived = found_derived.expect("derived name should be present");
    assert_eq!(found_derived.ancestors[0], new_base_id);
}

#[test]
fn compact_removes_dead_names_and_remaps() {
    let mut pool = NamePool::new();
    let a = pool.intern(make_name(1, RoleCategory::Cap, 0), None);
    let b = pool.intern(make_name(2, RoleCategory::Cap, 0), None);
    let _c = pool.intern(make_name(3, RoleCategory::Cap, 0), None);

    let mut keep = HashSet::new();
    keep.insert(a);
    keep.insert(b);

    let remap = pool.compact(&keep);

    assert_eq!(pool.len(), 2, "dead name c removed");
    // a and b should be mapped to consecutive new ids.
    let new_a = remap[&a];
    let new_b = remap[&b];
    assert!(new_a.0 < 2);
    assert!(new_b.0 < 2);
    assert_ne!(new_a, new_b);
    // c is not in the remap.
    assert!(!remap.contains_key(&_c));
}

#[test]
fn compact_remaps_ancestor_references() {
    let mut pool = NamePool::new();
    // Insert a to-be-dropped name at index 0, then two survivors.
    let _dead = pool.intern(make_name(0, RoleCategory::Cap, 0), None);
    let parent = pool.intern(make_name(1, RoleCategory::Cap, 0), None);
    let child_name = EntityName {
        operation: OperationId(2),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: vec![parent],
    };
    let child = pool.intern(child_name, None);

    let mut keep = HashSet::new();
    keep.insert(parent);
    keep.insert(child);
    let remap = pool.compact(&keep);

    let new_parent = remap[&parent];
    let new_child = remap[&child];
    let child_entry = pool.get(new_child).expect("child present");
    assert_eq!(child_entry.ancestors[0], new_parent);
}

#[test]
fn is_empty_and_len() {
    let mut pool = NamePool::new();
    assert!(pool.is_empty());
    pool.intern(make_name(1, RoleCategory::Modified, 0), None);
    assert!(!pool.is_empty());
    assert_eq!(pool.len(), 1);
}
