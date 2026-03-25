//! Integration tests exercising the full viz pipeline.
//!
//! Tests that require a GPU adapter (wgpu headless rendering) are marked
//! `#[ignore]` because they will fail in CI environments without a GPU.

use std::path::Path;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_viz::commands::{Command, execute_command, run_command_loop};
use popkern_viz::loader::{load_model, save_brep_json};
use popkern_viz::scene::VizServer;
use slotmap::Key;

// ---------------------------------------------------------------------------
// 1. BRep JSON round-trip through viz
// ---------------------------------------------------------------------------

#[test]
fn brep_json_roundtrip_preserves_topology() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");

    let dir = std::env::temp_dir();
    let path = dir.join("popkern_viz_roundtrip_test.json");

    save_brep_json(&brep, &path).expect("save should succeed");
    let loaded = load_model(&path).expect("load should succeed");

    assert_eq!(loaded.faces.len(), brep.faces.len(), "face count mismatch");
    assert_eq!(loaded.edges.len(), brep.edges.len(), "edge count mismatch");
    assert_eq!(
        loaded.vertices.len(),
        brep.vertices.len(),
        "vertex count mismatch"
    );

    std::fs::remove_file(&path).ok();
}

#[test]
fn brep_json_roundtrip_through_viz_server() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let original_faces = brep.faces.len();
    let original_edges = brep.edges.len();
    let original_verts = brep.vertices.len();

    let dir = std::env::temp_dir();
    let path = dir.join("popkern_viz_server_roundtrip.json");

    // Save, load, and feed into VizServer
    save_brep_json(&brep, &path).expect("save should succeed");
    let loaded = load_model(&path).expect("load should succeed");

    let mut scene = VizServer::new();
    scene.load_brep(loaded).expect("load into scene succeeds");

    let scene_brep = scene.brep.as_ref().expect("model should be loaded");
    assert_eq!(scene_brep.faces.len(), original_faces);
    assert_eq!(scene_brep.edges.len(), original_edges);
    assert_eq!(scene_brep.vertices.len(), original_verts);

    std::fs::remove_file(&path).ok();
}

// ---------------------------------------------------------------------------
// 2. Command protocol round-trip: select + inspect
// ---------------------------------------------------------------------------

#[test]
fn command_protocol_select_and_inspect_via_api() {
    let mut scene = VizServer::new();
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    scene.load_brep(brep).expect("load succeeds");

    // Get a real face handle from the loaded BRep to use as the inspect target.
    // Slotmap keys are not sequential integers, so we extract the raw ID from
    // the first face key in the map.
    let first_face_key = scene
        .brep
        .as_ref()
        .expect("model loaded")
        .faces
        .keys()
        .next()
        .expect("cube has faces");
    let raw_id = first_face_key.data().as_ffi() as u32;

    // Execute select
    let select_result = execute_command(
        &mut scene,
        &Command::Select {
            entity_type: "face".into(),
            id: raw_id,
        },
    );
    assert_eq!(select_result["ok"], true);
    assert_eq!(select_result["data"]["selected"], "face");

    // Execute inspect — should return surface type for a real face
    let inspect_result = execute_command(
        &mut scene,
        &Command::Inspect {
            entity_type: "face".into(),
            id: raw_id,
        },
    );
    assert_eq!(inspect_result["ok"], true);
    assert!(
        inspect_result["data"]["surface"].is_string(),
        "inspect response should contain surface type"
    );
}

// ---------------------------------------------------------------------------
// 3. Command protocol: validate on a valid cube
// ---------------------------------------------------------------------------

#[test]
fn command_protocol_validate_valid_cube() {
    let mut scene = VizServer::new();
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    scene.load_brep(brep).expect("load succeeds");

    let input = "validate\nquit\n";
    let mut output = Vec::new();

    run_command_loop(&mut scene, input.as_bytes(), &mut output).expect("loop runs");

    let output_str = String::from_utf8(output).expect("valid utf8");
    let lines: Vec<&str> = output_str.trim().lines().collect();
    assert_eq!(lines.len(), 2);

    let val: serde_json::Value = serde_json::from_str(lines[0]).expect("valid json");
    assert_eq!(val["ok"], true);
    assert_eq!(val["data"]["valid"], true);

    // Errors array should be empty for a valid cube
    let errors = val["data"]["errors"].as_array().expect("errors is array");
    assert!(errors.is_empty(), "valid cube should have no errors");
}

// ---------------------------------------------------------------------------
// 4. Model loader fallback: nonexistent file
// ---------------------------------------------------------------------------

#[test]
fn load_nonexistent_file_returns_graceful_error() {
    let path = Path::new("/tmp/popkern_definitely_does_not_exist_xyz.json");
    let result = load_model(path);
    assert!(result.is_err(), "loading nonexistent file should error");

    let err_msg = result.unwrap_err().to_string();
    // Should be an IO error, not a panic
    assert!(
        err_msg.contains("io error") || err_msg.contains("No such file"),
        "error should indicate file not found, got: {err_msg}"
    );
}

#[test]
fn load_nonexistent_non_json_graceful_error() {
    let path = Path::new("/tmp/popkern_definitely_does_not_exist_xyz.stl");
    let result = load_model(path);
    assert!(result.is_err(), "loading unsupported format should error");
}

// ---------------------------------------------------------------------------
// 5. Camera fit-to-bounds
// ---------------------------------------------------------------------------

#[test]
fn camera_fit_to_bounds_reasonable_distance() {
    let mut scene = VizServer::new();
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    scene.load_brep(brep).expect("load succeeds");

    // After loading, camera should have been fitted to the model
    let eye = scene.camera.eye;
    let target = scene.camera.target;
    let distance = (eye - target).length();

    // The unit cube spans roughly [0,1]^3 so the bounding sphere radius
    // is about 0.87. The camera distance should be greater than this radius
    // and not absurdly far away.
    assert!(
        distance > 0.5,
        "camera should be at least 0.5 from target, got {distance}"
    );
    assert!(
        distance < 100.0,
        "camera should not be more than 100 from target, got {distance}"
    );

    // Target should be near the center of the unit cube
    let center_dist = (target - glam::Vec3::new(0.5, 0.5, 0.5)).length();
    assert!(
        center_dist < 1.0,
        "target should be near cube center, got distance {center_dist}"
    );
}

#[test]
fn fit_command_updates_camera() {
    let mut scene = VizServer::new();
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    scene.load_brep(brep).expect("load succeeds");

    // Move camera far away
    scene.camera.eye = glam::Vec3::new(1000.0, 1000.0, 1000.0);
    let far_dist = (scene.camera.eye - scene.camera.target).length();

    // Run fit command
    let input = "fit\nquit\n";
    let mut output = Vec::new();
    run_command_loop(&mut scene, input.as_bytes(), &mut output).expect("loop runs");

    let new_dist = (scene.camera.eye - scene.camera.target).length();
    assert!(
        new_dist < far_dist,
        "fit should bring camera closer: was {far_dist}, now {new_dist}"
    );
}

// ---------------------------------------------------------------------------
// Multi-command session
// ---------------------------------------------------------------------------

#[test]
fn full_agent_session() {
    let mut scene = VizServer::new();
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    scene.load_brep(brep).expect("load succeeds");

    let input = "\
        wireframe off\n\
        colors flat\n\
        camera front\n\
        highlight face 0,1 red\n\
        validate\n\
        fit\n\
        wireframe on\n\
        quit\n";
    let mut output = Vec::new();

    run_command_loop(&mut scene, input.as_bytes(), &mut output).expect("loop runs");

    let output_str = String::from_utf8(output).expect("valid utf8");
    let lines: Vec<&str> = output_str.trim().lines().collect();
    assert_eq!(lines.len(), 8, "expected 8 responses");

    // Every response should be valid JSON with an "ok" field
    for (i, line) in lines.iter().enumerate() {
        let val: serde_json::Value =
            serde_json::from_str(line).unwrap_or_else(|e| panic!("line {i} invalid json: {e}"));
        assert!(
            val.get("ok").is_some(),
            "line {i} missing 'ok' field: {val}"
        );
    }

    // All commands should succeed
    for (i, line) in lines.iter().enumerate() {
        let val: serde_json::Value = serde_json::from_str(line).expect("valid json");
        assert_eq!(val["ok"], true, "line {i} should be ok=true: {val}");
    }
}

// ---------------------------------------------------------------------------
// Headless render tests (require GPU, marked #[ignore])
// ---------------------------------------------------------------------------

/// Render a unit cube to PNG headlessly.
///
/// Requires a wgpu-compatible GPU adapter. Will fail in CI without a GPU.
#[cfg(feature = "gpu")]
#[test]
#[ignore = "requires GPU adapter for wgpu headless rendering"]
fn headless_render_unit_cube_to_png() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let camera = popkern_viz::camera::Camera::default_iso();
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_headless_test.png");

    let result = popkern_viz::headless::render_to_png_gpu(&brep, &camera, 320, 240, &path, true);

    assert!(result.is_ok(), "headless render should succeed: {result:?}");
    assert!(path.exists(), "PNG file should be created");
    let meta = std::fs::metadata(&path).expect("metadata");
    assert!(meta.len() > 0, "PNG should not be empty");

    std::fs::remove_file(&path).ok();
}
