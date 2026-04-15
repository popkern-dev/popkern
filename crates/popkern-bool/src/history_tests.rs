use super::*;
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_pair() -> (BRep, BRep) {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    (a, b)
}

#[test]
fn boolean_recorded_captures_inputs_and_result() {
    let (a, b) = cube_pair();
    let a_face_count = a.faces.len();
    let b_face_count = b.faces.len();

    let record = boolean_recorded(BoolOperation::Union, &a, &b, tol()).unwrap();

    assert_eq!(record.op, BoolOperation::Union);
    assert_eq!(record.input_a.faces.len(), a_face_count);
    assert_eq!(record.input_b.faces.len(), b_face_count);
    assert!(!record.result.faces.is_empty());
}

#[test]
fn undo_state_returns_original_inputs() {
    let (a, b) = cube_pair();
    let original_a_vertices = a.vertices.len();
    let original_b_vertices = b.vertices.len();

    let record = boolean_recorded(BoolOperation::Subtraction, &a, &b, tol()).unwrap();

    let (undo_a, undo_b) = record.undo_state();
    assert_eq!(undo_a.vertices.len(), original_a_vertices);
    assert_eq!(undo_b.vertices.len(), original_b_vertices);
}

#[test]
fn redo_state_returns_cached_result() {
    let (a, b) = cube_pair();
    let record = boolean_recorded(BoolOperation::Intersection, &a, &b, tol()).unwrap();

    let cached = record.redo_state();
    // The cached result reference must be equivalent to the
    // record's internal result — same face count, same vertex count.
    assert_eq!(cached.faces.len(), record.result.faces.len());
    assert_eq!(cached.vertices.len(), record.result.vertices.len());
}

#[test]
fn re_execute_reproduces_cached_result() {
    let (a, b) = cube_pair();
    let record = boolean_recorded(BoolOperation::Union, &a, &b, tol()).unwrap();

    // Re-execute should produce a topologically-equivalent BRep
    // (deterministic boolean given identical inputs and tolerance).
    let replay = record.re_execute().unwrap();
    assert_eq!(replay.faces.len(), record.result.faces.len());
    assert_eq!(replay.vertices.len(), record.result.vertices.len());
    assert_eq!(replay.edges.len(), record.result.edges.len());
}

#[test]
fn record_is_clone_for_history_list_storage() {
    // The record must be Clone so an application can deep-copy it
    // into a snapshot store. This test verifies the Clone impl
    // produces an independent copy that doesn't alias the original.
    let (a, b) = cube_pair();
    let record = boolean_recorded(BoolOperation::Union, &a, &b, tol()).unwrap();
    let cloned = record.clone();
    assert_eq!(cloned.input_a.faces.len(), record.input_a.faces.len());
    assert_eq!(cloned.input_b.faces.len(), record.input_b.faces.len());
    assert_eq!(cloned.result.faces.len(), record.result.faces.len());
    assert_eq!(cloned.op, record.op);
}

#[test]
fn history_list_navigation() {
    // Simulate a 3-step history: union → intersection → subtraction.
    // Each step records a boolean and the application can undo to
    // any earlier state.
    let (a, b) = cube_pair();
    let mut history: Vec<BoolOperationRecord> = Vec::new();

    let r1 = boolean_recorded(BoolOperation::Union, &a, &b, tol()).unwrap();
    history.push(r1);

    let r2 = boolean_recorded(BoolOperation::Intersection, &a, &b, tol()).unwrap();
    history.push(r2);

    let r3 = boolean_recorded(BoolOperation::Subtraction, &a, &b, tol()).unwrap();
    history.push(r3);

    assert_eq!(history.len(), 3);
    // Most recent record's redo_state is the current model.
    let current = history.last().unwrap().redo_state();
    assert_eq!(
        current.faces.len(),
        history.last().unwrap().result.faces.len()
    );

    // Undoing the most recent step: the previous record's
    // redo_state becomes the current model.
    history.pop(); // simulate undo
    assert_eq!(history.len(), 2);
    let after_undo = history.last().unwrap().redo_state();
    assert!(!after_undo.faces.is_empty());

    // Each record's undo_state restores the inputs of THAT step.
    for record in &history {
        let (undo_a, undo_b) = record.undo_state();
        assert_eq!(undo_a.vertices.len(), 8);
        assert_eq!(undo_b.vertices.len(), 8);
    }
}

#[test]
fn record_re_execute_is_deterministic_across_calls() {
    let (a, b) = cube_pair();
    let record = boolean_recorded(BoolOperation::Subtraction, &a, &b, tol()).unwrap();

    // Re-execute three times; every result should be
    // topologically identical.
    let r1 = record.re_execute().unwrap();
    let r2 = record.re_execute().unwrap();
    let r3 = record.re_execute().unwrap();

    assert_eq!(r1.faces.len(), r2.faces.len());
    assert_eq!(r2.faces.len(), r3.faces.len());
    assert_eq!(r1.vertices.len(), r2.vertices.len());
    assert_eq!(r2.vertices.len(), r3.vertices.len());
    assert_eq!(r1.edges.len(), r2.edges.len());
    assert_eq!(r2.edges.len(), r3.edges.len());
}

#[test]
fn record_preserves_tolerance_for_re_execution() {
    // A non-default tolerance must round-trip through the record
    // so re_execute uses the same tolerance the original call used.
    let custom_tol = Tolerance::default().with_linear(1e-4);
    let (a, b) = cube_pair();
    let record = boolean_recorded(BoolOperation::Union, &a, &b, custom_tol).unwrap();

    assert_eq!(record.tolerance.linear(), 1e-4);
    let replay = record.re_execute().unwrap();
    assert!(!replay.faces.is_empty());
}

#[test]
fn boolean_recorded_propagates_errors() {
    // Attempting a boolean on degenerate input should propagate
    // the BoolError through boolean_recorded.
    let empty = BRep::new(tol());
    let cube = BRep::unit_cube(tol()).unwrap();
    let result = boolean_recorded(BoolOperation::Union, &empty, &cube, tol());
    // Either DegenerateInput or another error category — the
    // important property is that the call doesn't panic and
    // returns Err.
    assert!(result.is_err());
}

#[test]
fn record_input_clones_are_independent_of_originals() {
    // After creating a record, mutating the original BReps must
    // NOT affect the snapshots inside the record.
    let (mut a, b) = cube_pair();
    let record = boolean_recorded(BoolOperation::Union, &a, &b, tol()).unwrap();

    // Add a stray vertex to the original `a`.
    a.vertices.insert(popkern_core::topology::Vertex {
        position: popkern_core::geometry::Point3::new(99.0, 99.0, 99.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });

    // The record's snapshot is unaffected.
    assert_eq!(record.input_a.vertices.len(), 8);
    assert_eq!(a.vertices.len(), 9);
}
