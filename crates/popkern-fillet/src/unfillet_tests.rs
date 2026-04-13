use super::*;
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::collections::HashSet;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube() -> BRep {
    BRep::unit_cube(tol()).unwrap()
}

fn fillet_unit_cube_first_edge(radius: f64) -> (HashSet<FaceHandle>, BRep) {
    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let pre: HashSet<_> = brep.faces.keys().collect();
    let result = crate::fillet_edge(&brep, edge, radius, tol(), None).unwrap();
    (pre, result.brep)
}

fn find_blend_face(pre: &HashSet<FaceHandle>, post: &BRep) -> FaceHandle {
    post.faces
        .keys()
        .filter(|fh| !pre.contains(fh))
        .find(|fh| crate::identify_blend(post, *fh).is_some())
        .expect("post-fillet BRep should contain a recognised blend face")
}

// ---------------------------------------------------------------------------
// compute_unfillet_plan: positive cases
// ---------------------------------------------------------------------------

#[test]
fn compute_unfillet_plan_returns_plan_for_cube_fillet() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let blend = find_blend_face(&pre, &post);
    let plan = compute_unfillet_plan(&post, blend).expect("plan should compute");

    assert_eq!(plan.blend_face, blend);
    // Restored edge should have non-zero length.
    let dir = plan.restored_edge_direction;
    let dx = plan.restored_edge_end.x - plan.restored_edge_start.x;
    let dy = plan.restored_edge_end.y - plan.restored_edge_start.y;
    let dz = plan.restored_edge_end.z - plan.restored_edge_start.z;
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    assert!(
        len > 0.5 && len < 2.0,
        "restored edge length {len} should be ~1 for a unit cube"
    );
    // Direction vector is unit length.
    let dir_len = (dir.x * dir.x + dir.y * dir.y + dir.z * dir.z).sqrt();
    assert!((dir_len - 1.0).abs() < 1e-9);
}

#[test]
fn compute_unfillet_plan_identifies_two_supports() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let blend = find_blend_face(&pre, &post);
    let plan = compute_unfillet_plan(&post, blend).unwrap();
    assert_ne!(
        plan.support_a, plan.support_b,
        "supports must be two distinct faces"
    );
}

#[test]
fn compute_unfillet_plan_collects_edges_to_remove() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let blend = find_blend_face(&pre, &post);
    let plan = compute_unfillet_plan(&post, blend).unwrap();
    assert!(
        !plan.edges_to_remove.is_empty(),
        "blend face should have at least one boundary edge to remove"
    );
}

#[test]
fn compute_unfillet_plan_finds_long_boundary_to_at_least_one_support() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let blend = find_blend_face(&pre, &post);
    let plan = compute_unfillet_plan(&post, blend).unwrap();
    // The blend face must touch at least one of its detected
    // supports through a half-edge — otherwise the surgery would
    // have nothing to operate on.
    assert!(
        plan.long_he_to_a.is_some() || plan.long_he_to_b.is_some(),
        "expected at least one long boundary half-edge to a support"
    );
}

#[test]
fn compute_unfillet_plan_does_not_mutate_input() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let face_count_before = post.faces.len();
    let edge_count_before = post.edges.len();
    let vertex_count_before = post.vertices.len();
    let blend = find_blend_face(&pre, &post);
    let _ = compute_unfillet_plan(&post, blend).unwrap();
    assert_eq!(post.faces.len(), face_count_before);
    assert_eq!(post.edges.len(), edge_count_before);
    assert_eq!(post.vertices.len(), vertex_count_before);
}

// ---------------------------------------------------------------------------
// compute_unfillet_plan: negative cases
// ---------------------------------------------------------------------------

#[test]
fn compute_unfillet_plan_rejects_unknown_face() {
    let brep = cube();
    let bogus = FaceHandle::default();
    let err = compute_unfillet_plan(&brep, bogus).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

#[test]
fn compute_unfillet_plan_rejects_non_blend_face() {
    let brep = cube();
    let face = brep.faces.keys().next().unwrap();
    let err = compute_unfillet_plan(&brep, face).unwrap_err();
    let BlendError::InvalidParameter { reason, .. } = &err else {
        panic!("expected InvalidParameter, got {err:?}");
    };
    assert!(reason.contains("not recognised as a blend"));
}

// ---------------------------------------------------------------------------
// unfillet_face: stub returns the documented Round-5c hand-off
// ---------------------------------------------------------------------------

#[test]
fn unfillet_face_returns_documented_not_implemented_for_now() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let blend = find_blend_face(&pre, &post);
    let err = unfillet_face(&post, blend).unwrap_err();
    let BlendError::InvalidParameter { reason, .. } = &err else {
        panic!("expected InvalidParameter, got {err:?}");
    };
    // The remediation rationale must point users at the two
    // working alternatives so the hand-off is actionable.
    assert!(
        reason.contains("unfillet_via_record"),
        "rationale should mention the snapshot path: {reason}"
    );
    assert!(
        reason.contains("compute_unfillet_plan"),
        "rationale should mention the plan path: {reason}"
    );
}

#[test]
fn unfillet_face_propagates_input_validation_errors() {
    let brep = cube();
    let bogus = FaceHandle::default();
    let err = unfillet_face(&brep, bogus).unwrap_err();
    // Same shape as compute_unfillet_plan — invalid face is
    // surfaced the same way regardless of which entry point hit it.
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

// ---------------------------------------------------------------------------
// unfillet_via_record: the snapshot-based undo path
// ---------------------------------------------------------------------------

#[test]
fn unfillet_via_record_returns_input_brep() {
    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let pre_face_count = brep.faces.len();
    let pre_edge_count = brep.edges.len();
    let pre_vertex_count = brep.vertices.len();

    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let record = crate::fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();

    let restored = unfillet_via_record(&record);

    assert_eq!(restored.faces.len(), pre_face_count);
    assert_eq!(restored.edges.len(), pre_edge_count);
    assert_eq!(restored.vertices.len(), pre_vertex_count);
}

#[test]
fn unfillet_via_record_round_trips_through_re_execute() {
    // A FilletOperationRecord captures both ends of the round trip:
    // input (via undo_state / unfillet_via_record) and result
    // (via redo_state). re_execute should reproduce result from
    // input. This test confirms the two directions are consistent.
    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let record = crate::fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();

    let undo_via_helper = unfillet_via_record(&record);
    let undo_via_record_method = record.undo_state();

    // Both paths return the same BRep — identical face/edge/vertex
    // counts (and identical handles, since one is the same reference).
    assert_eq!(
        undo_via_helper.faces.len(),
        undo_via_record_method.faces.len()
    );
    assert_eq!(
        undo_via_helper.edges.len(),
        undo_via_record_method.edges.len()
    );
    assert_eq!(
        undo_via_helper.vertices.len(),
        undo_via_record_method.vertices.len()
    );

    // Re-executing the record from the unfilleted input should
    // reproduce the filleted result.
    let replay = record.re_execute().unwrap();
    assert_eq!(replay.faces.len(), record.result.faces.len());
}

#[test]
fn unfillet_via_record_works_for_multi_edge_fillets() {
    let brep = cube();
    let edges: Vec<_> = brep.edges.keys().take(3).collect();
    let specs: Vec<_> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();
    let record = crate::fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();

    let restored = unfillet_via_record(&record);
    assert_eq!(restored.faces.len(), brep.faces.len());
    assert_eq!(restored.edges.len(), brep.edges.len());
    assert_eq!(restored.vertices.len(), brep.vertices.len());
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

#[test]
fn compute_unfillet_plan_is_deterministic_across_calls() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let blend = find_blend_face(&pre, &post);
    let p1 = compute_unfillet_plan(&post, blend).unwrap();
    let p2 = compute_unfillet_plan(&post, blend).unwrap();

    let dx = (p1.restored_edge_start.x - p2.restored_edge_start.x).abs();
    let dy = (p1.restored_edge_start.y - p2.restored_edge_start.y).abs();
    let dz = (p1.restored_edge_start.z - p2.restored_edge_start.z).abs();
    assert!(dx + dy + dz < 1e-12);

    assert_eq!(p1.support_a, p2.support_a);
    assert_eq!(p1.support_b, p2.support_b);
}
