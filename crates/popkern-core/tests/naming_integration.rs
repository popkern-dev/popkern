//! Integration tests for stable topological naming.
//!
//! These tests exercise the full lifecycle: naming, history tracking,
//! resolution, serialization, pool merge, and compaction.

use std::collections::HashSet;

use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, NameId, OperationId, RoleCategory};
use popkern_core::naming_criteria::{GeometricCriteria, capture_face, capture_vertex};
use popkern_core::naming_resolve::{
    Resolution, ResolvedEntity, resolve_name, resolve_with_criteria,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

// ---------------------------------------------------------------------------
// Test 1: Full naming lifecycle
// ---------------------------------------------------------------------------

#[test]
fn full_naming_lifecycle() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");

    let op = OperationId(1);
    let mut ctx = TrackingContext::new(op);

    // Name all 6 faces.
    let face_handles: Vec<_> = cube.faces.keys().collect();
    assert_eq!(face_handles.len(), 6, "unit cube has 6 faces");

    for (i, &fh) in face_handles.iter().enumerate() {
        let name = EntityName {
            operation: op,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: vec![],
        };
        let criteria = capture_face(&cube, fh).map(GeometricCriteria::Face);
        let name_id = cube.name_pool.intern(name.clone(), criteria);
        cube.faces[fh].name = Some(name_id);
        ctx.record_generated(name);
    }

    let history = ctx.finish();
    assert_eq!(history.generated.len(), 6);
    assert_eq!(history.operation, op);

    // Resolve each face.
    for &fh in &face_handles {
        let face = &cube.faces[fh];
        let name_id = face.name.expect("face should be named");
        let name = cube.name_pool.get(name_id).expect("name should exist");
        let resolution = resolve_name(name, &cube);
        match resolution {
            Resolution::Exact(ResolvedEntity::Face(h)) => {
                assert_eq!(h, fh, "resolved to wrong face");
            }
            other => panic!("expected Exact(Face), got {:?}", other),
        }
    }
}

// ---------------------------------------------------------------------------
// Test 2: Zero overhead — no naming
// ---------------------------------------------------------------------------

#[test]
fn zero_overhead_no_naming() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");

    // All name fields should be None.
    for (_, vertex) in &cube.vertices {
        assert!(vertex.name.is_none(), "vertex should have no name");
    }
    for (_, edge) in &cube.edges {
        assert!(edge.name.is_none(), "edge should have no name");
    }
    for (_, face) in &cube.faces {
        assert!(face.name.is_none(), "face should have no name");
    }

    // Pool should be empty.
    assert!(cube.name_pool.is_empty(), "name_pool should be empty");
}

// ---------------------------------------------------------------------------
// Test 3: Serialization round-trip
// ---------------------------------------------------------------------------

#[test]
fn serialization_round_trip() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");

    // Name all faces.
    let op = OperationId(10);
    let face_handles: Vec<_> = cube.faces.keys().collect();
    let mut names: Vec<EntityName> = Vec::new();

    for (i, &fh) in face_handles.iter().enumerate() {
        let name = EntityName {
            operation: op,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: vec![],
        };
        let criteria = capture_face(&cube, fh).map(GeometricCriteria::Face);
        let nid = cube.name_pool.intern(name.clone(), criteria);
        cube.faces[fh].name = Some(nid);
        names.push(name);
    }

    // Serialize → deserialize.
    let json = serde_json::to_string(&cube).expect("serialize");
    let mut cube2: BRep = serde_json::from_str(&json).expect("deserialize");
    cube2.name_pool.rebuild_index();

    // Verify all names still resolve.
    for name in &names {
        let resolution = resolve_name(name, &cube2);
        assert!(
            matches!(resolution, Resolution::Exact(ResolvedEntity::Face(_))),
            "name {:?} should resolve after round-trip, got {:?}",
            name.role.index,
            resolution
        );
    }
}

// ---------------------------------------------------------------------------
// Test 4: Merge pools
// ---------------------------------------------------------------------------

#[test]
fn merge_pools() {
    let mut cube_a = BRep::unit_cube(Tolerance::default()).expect("cube A");
    let mut cube_b = BRep::unit_cube(Tolerance::default()).expect("cube B");

    // Name faces in cube A with op 1.
    let op_a = OperationId(1);
    let faces_a: Vec<_> = cube_a.faces.keys().collect();
    let mut names_a = Vec::new();
    for (i, &fh) in faces_a.iter().enumerate() {
        let name = EntityName {
            operation: op_a,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: vec![],
        };
        let nid = cube_a.name_pool.intern(name.clone(), None);
        cube_a.faces[fh].name = Some(nid);
        names_a.push(name);
    }

    // Name faces in cube B with op 2.
    let op_b = OperationId(2);
    let faces_b: Vec<_> = cube_b.faces.keys().collect();
    let mut names_b = Vec::new();
    for (i, &fh) in faces_b.iter().enumerate() {
        let name = EntityName {
            operation: op_b,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: vec![],
        };
        let nid = cube_b.name_pool.intern(name.clone(), None);
        cube_b.faces[fh].name = Some(nid);
        names_b.push(name);
    }

    // Merge B's pool into A.
    let remap = cube_a.name_pool.merge(&cube_b.name_pool);
    assert_eq!(cube_a.name_pool.len(), 12, "6 from A + 6 from B");

    // Remap B's face names and move them into A's arena.
    for &fh in &faces_b {
        let old_nid = cube_b.faces[fh].name.expect("named");
        let new_nid = remap[&old_nid];
        // Verify the name is findable in A's pool.
        let stored = cube_a.name_pool.get(new_nid).expect("name in merged pool");
        assert_eq!(stored.operation, op_b);
    }

    // All original A names are still findable.
    for name in &names_a {
        assert!(
            cube_a.name_pool.find(name).is_some(),
            "name from A should survive merge"
        );
    }

    // All B names (after remapping) are findable.
    for name in &names_b {
        assert!(
            cube_a.name_pool.find(name).is_some(),
            "name from B should be in merged pool"
        );
    }
}

// ---------------------------------------------------------------------------
// Test 5: Compact
// ---------------------------------------------------------------------------

#[test]
fn compact_removes_dead_names() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");

    // Name all 6 faces.
    let op = OperationId(5);
    let face_handles: Vec<_> = cube.faces.keys().collect();
    let mut name_ids = Vec::new();

    for (i, &fh) in face_handles.iter().enumerate() {
        let name = EntityName {
            operation: op,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: vec![],
        };
        let nid = cube.name_pool.intern(name, None);
        cube.faces[fh].name = Some(nid);
        name_ids.push(nid);
    }
    assert_eq!(cube.name_pool.len(), 6);

    // "Delete" two faces by clearing their names.
    cube.faces[face_handles[0]].name = None;
    cube.faces[face_handles[1]].name = None;

    // Collect referenced ids.
    let referenced: HashSet<NameId> = cube.faces.values().filter_map(|f| f.name).collect();
    assert_eq!(referenced.len(), 4);

    // Compact.
    let remap = cube.name_pool.compact(&referenced);
    assert_eq!(cube.name_pool.len(), 4, "two dead names removed");

    // Remap surviving entity name fields.
    for (_, face) in &mut cube.faces {
        if let Some(ref mut nid) = face.name {
            if let Some(&new) = remap.get(nid) {
                *nid = new;
            }
        }
    }

    // Surviving faces still resolve.
    for &fh in &face_handles[2..] {
        let nid = cube.faces[fh].name.expect("surviving face should be named");
        let name = cube
            .name_pool
            .get(nid)
            .expect("name should exist after compact");
        let resolution = resolve_name(name, &cube);
        assert!(
            matches!(resolution, Resolution::Exact(ResolvedEntity::Face(_))),
            "surviving name should resolve, got {:?}",
            resolution
        );
    }
}

// ---------------------------------------------------------------------------
// Test 6: Geometric fallback
// ---------------------------------------------------------------------------

#[test]
fn geometric_fallback_finds_face_after_pool_clear() {
    let mut cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let fh = cube.faces.keys().next().expect("face");

    // Capture criteria before anything else.
    let fc = capture_face(&cube, fh).expect("face criteria");

    // Name the face.
    let op = OperationId(6);
    let name = EntityName {
        operation: op,
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: vec![],
    };
    let nid = cube
        .name_pool
        .intern(name, Some(GeometricCriteria::Face(fc.clone())));
    cube.faces[fh].name = Some(nid);

    // Clear the name from the entity AND rebuild a fresh pool (simulate rebuild).
    cube.faces[fh].name = None;

    // Use resolve_with_criteria directly — no pool lookup needed.
    let resolution = resolve_with_criteria(&GeometricCriteria::Face(fc), &cube, &cube.tolerance);
    assert!(
        matches!(resolution, Resolution::Exact(ResolvedEntity::Face(_))),
        "geometric fallback should find the face, got {:?}",
        resolution
    );
}

#[test]
fn geometric_fallback_finds_vertex_by_position() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let vh = cube.vertices.keys().next().expect("vertex");

    let vc = capture_vertex(&cube, vh).expect("vertex criteria");

    let resolution = resolve_with_criteria(&GeometricCriteria::Vertex(vc), &cube, &cube.tolerance);
    assert!(
        matches!(resolution, Resolution::Exact(ResolvedEntity::Vertex(_))),
        "geometric fallback should find vertex by position, got {:?}",
        resolution
    );
}
