//! Additional scene (VizServer) edge-case tests.

use super::*;
use popkern_core::tolerance::Tolerance;

#[test]
fn load_empty_brep_succeeds() {
    let mut server = VizServer::new();
    let brep = BRep::new(Tolerance::default());
    // Empty BRep has no faces -> tessellation produces an empty TriMesh.
    // This should succeed without error.
    let result = server.load_brep(brep);
    assert!(
        result.is_ok(),
        "loading empty brep should succeed: {result:?}"
    );
    assert!(server.brep.is_some());
    assert!(server.mesh.is_some());
    // Bounds should be None since mesh has no vertices.
    assert!(server.model_bounds().is_none());
}

#[test]
fn load_brep_resets_selection() {
    let mut server = VizServer::new();
    server.select(Selection::Face(FaceHandle::default()));
    assert!(matches!(server.selection, Selection::Face(_)));

    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    server.load_brep(brep).expect("load succeeds");
    // Selection should be reset to None after loading.
    assert!(matches!(server.selection, Selection::None));
}

#[test]
fn load_brep_fits_camera_to_model() {
    let mut server = VizServer::new();
    let initial_eye = server.camera.eye;
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    server.load_brep(brep).expect("load succeeds");
    // Camera should have been adjusted (eye moved).
    assert_ne!(
        server.camera.eye, initial_eye,
        "camera should be adjusted after load"
    );
}

#[test]
fn load_brep_twice_replaces_model() {
    let mut server = VizServer::new();
    let brep1 = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    server.load_brep(brep1).expect("first load succeeds");
    let first_face_count = server.brep.as_ref().map(|b| b.faces.len()).unwrap_or(0);

    let brep2 = BRep::new(Tolerance::default());
    server.load_brep(brep2).expect("second load succeeds");
    let second_face_count = server.brep.as_ref().map(|b| b.faces.len()).unwrap_or(0);

    assert_ne!(
        first_face_count, second_face_count,
        "second load should replace the model"
    );
    assert_eq!(second_face_count, 0);
}

#[test]
fn display_options_debug_and_clone() {
    let opts = DisplayOptions::default();
    let opts2 = opts.clone();
    assert_eq!(opts.show_wireframe, opts2.show_wireframe);
    assert_eq!(opts.show_normals, opts2.show_normals);
    assert_eq!(opts.color_mode, opts2.color_mode);
    // Verify Debug impl compiles.
    let s = format!("{:?}", opts2);
    assert!(s.contains("show_wireframe"));
}

#[test]
fn selection_debug_impl() {
    let sel = Selection::None;
    let s = format!("{:?}", sel);
    assert!(s.contains("None"));

    let sel2 = Selection::Face(FaceHandle::default());
    let s2 = format!("{:?}", sel2);
    assert!(s2.contains("Face"));
}
