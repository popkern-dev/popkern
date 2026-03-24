//! Integration tests for the command protocol.

use popkern_viz::commands::{CameraArgs, Command, parse_command, run_command_loop};
use popkern_viz::scene::VizServer;

#[test]
fn parse_all_command_variants() {
    let cases = vec![
        ("select face 3", true),
        ("camera front", true),
        ("camera 0 0 5 0 0 0 0 1 0", true),
        ("screenshot /tmp/view.png", true),
        ("inspect face 3", true),
        ("highlight face 1,2,3 red", true),
        ("validate", true),
        ("fit", true),
        ("wireframe on", true),
        ("wireframe off", true),
        ("colors surface", true),
        ("quit", true),
        ("unknown_cmd", false),
        ("", false),
    ];

    for (input, should_succeed) in cases {
        let result = parse_command(input);
        assert_eq!(
            result.is_ok(),
            should_succeed,
            "input={input:?}, result={result:?}"
        );
    }
}

#[test]
fn camera_preset_roundtrip() {
    let cmd = parse_command("camera top").unwrap();
    assert_eq!(
        cmd,
        Command::Camera {
            args: CameraArgs::Preset("top".into())
        }
    );
}

#[test]
fn camera_explicit_roundtrip() {
    let cmd = parse_command("camera 1.5 2.5 3.5 0 0 0 0 1 0").unwrap();
    match cmd {
        Command::Camera {
            args: CameraArgs::Explicit { eye, target, up },
        } => {
            assert!((eye[0] - 1.5).abs() < f32::EPSILON);
            assert!((eye[1] - 2.5).abs() < f32::EPSILON);
            assert!((eye[2] - 3.5).abs() < f32::EPSILON);
            assert!((target[0]).abs() < f32::EPSILON);
            assert!((up[1] - 1.0).abs() < f32::EPSILON);
        }
        other => panic!("expected Camera::Explicit, got {other:?}"),
    }
}

#[test]
fn highlight_parses_multiple_ids() {
    let cmd = parse_command("highlight edge 10,20,30 #ff0000").unwrap();
    assert_eq!(
        cmd,
        Command::Highlight {
            entity_type: "edge".into(),
            ids: vec![10, 20, 30],
            color: "#ff0000".into(),
        }
    );
}

#[test]
fn command_loop_full_session() {
    let input = "wireframe off\ncolors flat\nvalidate\nquit\n";
    let mut output = Vec::new();
    let mut scene = VizServer::new();

    run_command_loop(&mut scene, input.as_bytes(), &mut output).unwrap();

    let output_str = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = output_str.trim().lines().collect();
    assert_eq!(lines.len(), 4);

    // All should be valid JSON
    for line in &lines {
        let val: serde_json::Value = serde_json::from_str(line).unwrap();
        assert!(val.get("ok").is_some(), "response missing 'ok': {val}");
    }
}

#[test]
fn command_loop_with_loaded_model() {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let mut scene = VizServer::new();
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    scene.load_brep(brep).expect("load succeeds");

    let input = "validate\nfit\nquit\n";
    let mut output = Vec::new();

    run_command_loop(&mut scene, input.as_bytes(), &mut output).unwrap();

    let output_str = String::from_utf8(output).unwrap();
    let lines: Vec<&str> = output_str.trim().lines().collect();
    assert_eq!(lines.len(), 3);

    // validate should succeed
    let val: serde_json::Value = serde_json::from_str(lines[0]).unwrap();
    assert_eq!(val["ok"], true);

    // fit should succeed
    let fit: serde_json::Value = serde_json::from_str(lines[1]).unwrap();
    assert_eq!(fit["ok"], true);
}
