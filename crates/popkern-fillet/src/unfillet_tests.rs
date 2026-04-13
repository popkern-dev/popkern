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
// unfillet_face: real algorithmic surgery (Round 5c-1)
// ---------------------------------------------------------------------------

#[test]
fn unfillet_face_round_trip_restores_unit_cube_face_count() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let blend = find_blend_face(&pre, &post);
    let restored = unfillet_face(&post, blend).expect("unfillet_face should succeed");
    // The original cube has 6 faces, the filleted result has 7.
    // Unfillet should bring it back to 6.
    assert_eq!(
        restored.faces.len(),
        6,
        "unfillet should restore the cube to 6 faces, got {}",
        restored.faces.len()
    );
}

#[test]
fn unfillet_face_does_not_mutate_input_brep() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let blend = find_blend_face(&pre, &post);
    let pre_face_count = post.faces.len();
    let pre_edge_count = post.edges.len();
    let pre_vertex_count = post.vertices.len();
    let _ = unfillet_face(&post, blend).unwrap();
    // Input must be unchanged — surgery happens on a clone.
    assert_eq!(post.faces.len(), pre_face_count);
    assert_eq!(post.edges.len(), pre_edge_count);
    assert_eq!(post.vertices.len(), pre_vertex_count);
}

#[test]
fn unfillet_face_result_passes_brep_validation() {
    let (pre, post) = fillet_unit_cube_first_edge(0.1);
    let blend = find_blend_face(&pre, &post);
    let restored = unfillet_face(&post, blend).unwrap();
    // The restored BRep must validate as a closed solid.
    let validation = restored.validate();
    assert!(
        validation.is_ok(),
        "unfilleted result should validate, got {:?}",
        validation.err()
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

#[test]
fn unfillet_face_rejects_non_blend_face() {
    let brep = cube();
    let face = brep.faces.keys().next().unwrap();
    let err = unfillet_face(&brep, face).unwrap_err();
    let BlendError::InvalidParameter { reason, .. } = &err else {
        panic!("expected InvalidParameter, got {err:?}");
    };
    assert!(reason.contains("not recognised as a blend"));
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
// Round 5c spike: can boolean union do unfillet?
// ---------------------------------------------------------------------------

/// Spike test: try to unfillet a unit-cube fillet by union'ing the
/// filleted result with a fresh unit cube. If the boolean engine
/// handles "input ⊆ supercube" cleanly, the result should be the
/// fresh cube.
///
/// This is the **simplest possible** shape of the boolean-union
/// approach to unfillet: instead of building a precise filler
/// solid, just union with the original unrounded body. If it
/// works, Round 5c-1 has a clean foundation; if it fails the same
/// way as `fillet_edge_on_boolean_cube_result`, the boolean fix
/// is the prerequisite.
///
/// **Outcome (2026-04-13):** the boolean engine accepts the call
/// but returns 7 faces (same as the filleted input), failing to
/// resolve the coincident-face overlap. ALSO ran for ~10 minutes
/// per call, indicating a separate boolean-engine performance
/// issue. Marked `#[ignore]` so it doesn't slow CI.
#[test]
#[ignore = "boolean engine takes ~10min per call and returns wrong result; documented in commit history"]
fn spike_5c_unfillet_via_union_with_fresh_cube() {
    use popkern_bool::{BoolOperation, boolean};

    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let filleted = crate::fillet_edge(&brep, edge, 0.1, tol(), None)
        .unwrap()
        .brep;

    // Fresh unit cube — same shape, no fillet.
    let fresh = cube();

    eprintln!(
        "spike: filleted has {} faces, fresh has {} faces",
        filleted.faces.len(),
        fresh.faces.len()
    );

    let result = boolean(BoolOperation::Union, fresh, filleted, tol(), None);
    match result {
        Ok(merged) => {
            eprintln!("spike: union OK, result has {} faces", merged.faces.len());
            // If the union worked perfectly, the result should be
            // the original cube (6 faces). If it has more faces,
            // the boolean engine produced extra geometry but didn't
            // collapse the rounded notch.
            //
            // We don't assert face_count == 6 here because the
            // spike's purpose is to *observe* what happens, not
            // pass/fail. The eprintln output is the test result.
        }
        Err(e) => {
            eprintln!("spike: union FAILED: {e:?}");
        }
    }
}

/// Variant: union a small box positioned to fill the rounded notch.
/// More precise than the supercube approach — only the notch
/// region overlaps.
///
/// **Outcome (2026-04-13):** failed with `SplitFailed { reason:
/// "trim curve end point does not match any vertex on the face
/// boundary" }`, same precision-fragility family as the known
/// `fillet_edge_on_boolean_cube_result` bug. ~10 minutes per
/// call. Marked `#[ignore]`.
#[test]
#[ignore = "boolean engine fails with SplitFailed and runs for ~10min; same bug family as fillet_edge_on_boolean_cube_result"]
fn spike_5c_unfillet_via_union_with_filler_box() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;

    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let filleted = crate::fillet_edge(&brep, edge, 0.1, tol(), None)
        .unwrap()
        .brep;

    // Locate the blend face and find its vertices' bounding box —
    // that gives us the spatial extent of the notch.
    let pre: HashSet<_> = brep.faces.keys().collect();
    let blend = find_blend_face(&pre, &filleted);
    let plan = compute_unfillet_plan(&filleted, blend).unwrap();
    eprintln!(
        "spike: blend face is {:?}, restored edge from {:?} to {:?}",
        blend, plan.restored_edge_start, plan.restored_edge_end
    );

    // Build a filler box: 0.2 × 0.2 cross-section running along
    // the restored edge direction, oversized so it definitely
    // covers the notch.
    let unit = cube();
    let scaled = unit.transform(&Transform3::scale(1.2, 0.25, 0.25)).unwrap();
    let translated = scaled
        .transform(&Transform3::translation(-0.1, -0.1, -0.1))
        .unwrap();

    eprintln!(
        "spike: filler box has {} faces, {} vertices",
        translated.faces.len(),
        translated.vertices.len()
    );

    let result = boolean(BoolOperation::Union, filleted, translated, tol(), None);
    match result {
        Ok(merged) => {
            eprintln!(
                "spike: filler-box union OK, result has {} faces",
                merged.faces.len()
            );
        }
        Err(e) => {
            eprintln!("spike: filler-box union FAILED: {e:?}");
        }
    }
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
