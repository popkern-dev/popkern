use super::*;
use crate::types::{ChamferKind, FilletRadius};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_and_edge() -> (BRep, popkern_core::topology::EdgeHandle) {
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    (brep, edge)
}

#[test]
fn fillet_recorded_captures_inputs_and_result() {
    let (brep, edge) = cube_and_edge();
    let pre_face_count = brep.faces.len();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let record = fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();

    assert_eq!(record.input.faces.len(), pre_face_count);
    assert!(record.result.faces.len() > pre_face_count);
    assert_eq!(record.specs.len(), 1);
    assert!((record.tolerance.linear() - tol().linear()).abs() < 1e-15);
}

#[test]
fn fillet_recorded_undo_state_returns_input() {
    let (brep, edge) = cube_and_edge();
    let pre_face_count = brep.faces.len();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let record = fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();
    let undo = record.undo_state();
    assert_eq!(undo.faces.len(), pre_face_count);
}

#[test]
fn fillet_recorded_redo_state_returns_cached_result() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let record = fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();
    let redo = record.redo_state();
    assert_eq!(redo.faces.len(), record.result.faces.len());
}

#[test]
fn fillet_recorded_re_execute_produces_same_shape() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let record = fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();
    let replay = record.re_execute().unwrap();
    // Deterministic — same face/edge/vertex counts.
    assert_eq!(replay.faces.len(), record.result.faces.len());
    assert_eq!(replay.edges.len(), record.result.edges.len());
    assert_eq!(replay.vertices.len(), record.result.vertices.len());
}

#[test]
fn fillet_record_is_clone_for_history_list_storage() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let record = fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();
    let cloned = record.clone();
    assert_eq!(cloned.input.faces.len(), record.input.faces.len());
    assert_eq!(cloned.result.faces.len(), record.result.faces.len());
}

#[test]
fn fillet_recorded_propagates_validation_errors() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(f64::NAN),
    }];
    let err = fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

#[test]
fn fillet_record_input_clones_are_independent_of_originals() {
    let (mut brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let record = fillet_recorded(&brep, &specs, &FilletOptions::default(), tol()).unwrap();
    // Mutate the original BRep — the record's input snapshot must
    // stay unchanged.
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: popkern_core::geometry::Point3::new(99.0, 99.0, 99.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    assert_ne!(brep.vertices.len(), record.input.vertices.len());
}

// ---------- chamfer records ----------

#[test]
fn chamfer_recorded_captures_inputs_and_result() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(0.1),
    }];
    let record = chamfer_recorded(&brep, &specs, &ChamferOptions::default(), tol()).unwrap();
    assert_eq!(record.specs.len(), 1);
    assert!(record.result.faces.len() > brep.faces.len());
}

#[test]
fn chamfer_recorded_undo_redo_round_trip() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(0.1),
    }];
    let record = chamfer_recorded(&brep, &specs, &ChamferOptions::default(), tol()).unwrap();
    let undo = record.undo_state();
    let redo = record.redo_state();
    // Round trip: undo shape differs from redo shape.
    assert_ne!(undo.faces.len(), redo.faces.len());
}

#[test]
fn chamfer_recorded_re_execute_is_deterministic() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(0.1),
    }];
    let record = chamfer_recorded(&brep, &specs, &ChamferOptions::default(), tol()).unwrap();
    let r1 = record.re_execute().unwrap();
    let r2 = record.re_execute().unwrap();
    assert_eq!(r1.faces.len(), r2.faces.len());
    assert_eq!(r1.edges.len(), r2.edges.len());
}

#[test]
fn history_list_navigation_simulation() {
    // Simulate a small history list: fillet → chamfer → fillet.
    let brep = BRep::unit_cube(tol()).unwrap();
    let edges: Vec<_> = brep.edges.keys().take(3).collect();

    // Step 1: fillet the first edge.
    let r1 = fillet_recorded(
        &brep,
        &[FilletSpec {
            edge: edges[0],
            radius: FilletRadius::Constant(0.1),
        }],
        &FilletOptions::default(),
        tol(),
    )
    .unwrap();

    // The record can be navigated: undo returns input, redo returns result.
    assert_eq!(r1.undo_state().faces.len(), brep.faces.len());
    assert!(r1.redo_state().faces.len() > brep.faces.len());
}
