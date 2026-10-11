//! Integration tests for extrude naming — verifies that the extrude Euler
//! operator correctly assigns stable names to lateral faces, caps, edges,
//! and vertices when a `TrackingContext` is provided.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use std::collections::HashSet;

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Find the bottom face after extrude: the face not in side_faces and not cap.
fn find_bottom_face(brep: &BRep, result: &popkern_core::euler::ExtrudeResult) -> FaceHandle {
    let mut accounted: HashSet<FaceHandle> = result.side_faces.iter().copied().collect();
    accounted.insert(result.cap_face);
    brep.faces
        .keys()
        .find(|fh| !accounted.contains(fh))
        .expect("bottom face should exist")
}

/// Build a square face and return (brep, face_handle).
fn build_square() -> (BRep, popkern_core::topology::FaceHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
    (brep, mef.face)
}

/// Build a triangle face and return (brep, face_handle).
fn build_triangle() -> (BRep, popkern_core::topology::FaceHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(0.5, 1.0, 0.0))
        .unwrap();
    let mef = brep.mef(mev2.half_edge_twin, mev1.half_edge).unwrap();
    (brep, mef.face)
}

// ---------------------------------------------------------------------------
// Test: Extrude named face produces named entities with ancestry
// ---------------------------------------------------------------------------

#[test]
fn extrude_named_face_produces_named_entities() {
    let (mut brep, face) = build_square();

    // Name the input face (simulating a prior operation).
    let prior_op = OperationId(1);
    let face_name = EntityName {
        operation: prior_op,
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: vec![],
    };
    let face_name_id = brep.name_pool.intern(face_name.clone(), None);
    brep.faces[face].name = Some(face_name_id);

    // Extrude with tracking.
    let extrude_op = OperationId(2);
    let mut ctx = TrackingContext::new(extrude_op);
    let result = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, Some(&mut ctx))
        .unwrap();
    let history = ctx.finish();

    // 4 lateral faces with Lateral role.
    assert_eq!(result.side_faces.len(), 4);
    for (i, &side_fh) in result.side_faces.iter().enumerate() {
        let nid = brep.faces[side_fh].name.expect("side face named");
        let name = brep.name_pool.get(nid).expect("name exists");
        assert_eq!(name.operation, extrude_op);
        assert_eq!(name.role.category, RoleCategory::Lateral);
        assert_eq!(name.role.index, i as u32);
    }

    // Cap face with Cap role and input face ancestor.
    let cap_nid = brep.faces[result.cap_face].name.expect("cap named");
    let cap_name = brep.name_pool.get(cap_nid).expect("cap name exists");
    assert_eq!(cap_name.operation, extrude_op);
    assert_eq!(cap_name.role.category, RoleCategory::Cap);
    assert_eq!(cap_name.ancestors.len(), 1);
    assert_eq!(cap_name.ancestors[0], face_name_id);

    // Bottom face with Modified role and input face ancestor.
    // Note: the input `face` handle becomes a side face after extrude.
    // The bottom face is the OTHER pre-existing face.
    let bottom_fh = find_bottom_face(&brep, &result);
    let bottom_nid = brep.faces[bottom_fh].name.expect("bottom named");
    let bottom_name = brep.name_pool.get(bottom_nid).expect("bottom exists");
    assert_eq!(bottom_name.operation, extrude_op);
    assert_eq!(bottom_name.role.category, RoleCategory::Modified);
    assert_eq!(bottom_name.ancestors.len(), 1);
    assert_eq!(bottom_name.ancestors[0], face_name_id);

    // New vertices named with Generated role.
    assert_eq!(result.new_vertices.len(), 4);
    for &vh in &result.new_vertices {
        let nid = brep.vertices[vh].name.expect("vertex named");
        let name = brep.name_pool.get(nid).unwrap();
        assert_eq!(name.operation, extrude_op);
        assert_eq!(name.role.category, RoleCategory::Generated);
    }

    // New edges named (8 = 4 vertical + 4 top).
    assert_eq!(result.new_edges.len(), 8);
    for &eh in &result.new_edges {
        assert!(brep.edges[eh].name.is_some(), "edge named");
    }

    // History has generated and modified entries.
    assert!(!history.generated.is_empty());
    assert!(!history.modified.is_empty());
}

// ---------------------------------------------------------------------------
// Test: Extrude without tracking produces no names
// ---------------------------------------------------------------------------

#[test]
fn extrude_without_tracking_produces_no_names() {
    let (mut brep, face) = build_square();

    let _result = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();

    for (_, f) in &brep.faces {
        assert!(f.name.is_none());
    }
    for (_, e) in &brep.edges {
        assert!(e.name.is_none());
    }
    for (_, v) in &brep.vertices {
        assert!(v.name.is_none());
    }
    assert!(brep.name_pool.is_empty());
}

// ---------------------------------------------------------------------------
// Test: Extrude triangle — correct counts
// ---------------------------------------------------------------------------

#[test]
fn extrude_triangle_naming_correct_counts() {
    let (mut brep, face) = build_triangle();

    let extrude_op = OperationId(10);
    let mut ctx = TrackingContext::new(extrude_op);
    let result = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, Some(&mut ctx))
        .unwrap();
    let history = ctx.finish();

    // 3 lateral faces.
    assert_eq!(result.side_faces.len(), 3);
    for &fh in &result.side_faces {
        assert!(brep.faces[fh].name.is_some());
    }

    // 1 cap.
    assert!(brep.faces[result.cap_face].name.is_some());

    // 3 new vertices + 6 new edges.
    assert_eq!(result.new_vertices.len(), 3);
    assert_eq!(result.new_edges.len(), 6);
    for &eh in &result.new_edges {
        assert!(brep.edges[eh].name.is_some());
    }

    // Bottom face named Modified with empty ancestors (was unnamed).
    let bottom_fh = find_bottom_face(&brep, &result);
    let bottom = brep
        .name_pool
        .get(brep.faces[bottom_fh].name.unwrap())
        .unwrap();
    assert_eq!(bottom.role.category, RoleCategory::Modified);
    assert!(bottom.ancestors.is_empty());

    // History: 3 laterals + 1 cap + 3 vertices + 6 edges + 1 bottom = 14.
    assert_eq!(history.generated.len(), 14);
}

// ---------------------------------------------------------------------------
// Test: Vertical edges have Lateral role, top edges have Cap role
// ---------------------------------------------------------------------------

#[test]
fn extrude_edge_roles_are_correct() {
    let (mut brep, face) = build_square();

    let op = OperationId(20);
    let mut ctx = TrackingContext::new(op);
    let result = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, Some(&mut ctx))
        .unwrap();
    let _history = ctx.finish();

    // First 4 edges are vertical (Lateral role).
    for i in 0..4 {
        let eh = result.new_edges[i];
        let nid = brep.edges[eh].name.unwrap();
        let name = brep.name_pool.get(nid).unwrap();
        assert_eq!(name.role.category, RoleCategory::Lateral, "edge {i}");
    }

    // Next 4 edges are top (Cap role).
    for i in 4..8 {
        let eh = result.new_edges[i];
        let nid = brep.edges[eh].name.unwrap();
        let name = brep.name_pool.get(nid).unwrap();
        assert_eq!(name.role.category, RoleCategory::Cap, "edge {i}");
    }
}
