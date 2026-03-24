//! Command execution: dispatch [`Command`] variants against the [`VizServer`].

use serde_json::json;
use slotmap::KeyData;

use crate::camera::ViewPreset;
use crate::scene::{ColorMode, Selection, VizServer};

use super::{CameraArgs, Command};

/// Execute a command against the visualization server.
///
/// Returns a JSON value with `{ "ok": true/false, ... }` indicating
/// success or failure along with any result data.
pub fn execute_command(scene: &mut VizServer, cmd: &Command) -> serde_json::Value {
    match cmd {
        Command::Select { entity_type, id } => exec_select(scene, entity_type, *id),
        Command::Camera { args } => exec_camera(scene, args),
        Command::Screenshot { path } => exec_screenshot(path),
        Command::Inspect { entity_type, id } => exec_inspect(scene, entity_type, *id),
        Command::Highlight {
            entity_type,
            ids,
            color,
        } => exec_highlight(entity_type, ids, color),
        Command::Validate => exec_validate(scene),
        Command::Fit => exec_fit(scene),
        Command::Wireframe { on } => exec_wireframe(scene, *on),
        Command::Colors { mode } => exec_colors(scene, mode),
        Command::Quit => json!({"ok": true, "data": {"message": "quit"}}),
    }
}

/// Build an OK JSON response.
fn ok_response(data: serde_json::Value) -> serde_json::Value {
    json!({"ok": true, "data": data})
}

/// Build an error JSON response.
pub(super) fn err_response(msg: &str) -> serde_json::Value {
    json!({"ok": false, "error": msg})
}

/// Construct a slotmap key from a raw u32 index.
///
/// Uses version 1 (the first valid version for slotmap keys).
fn key_data(id: u32) -> KeyData {
    KeyData::from_ffi(id as u64 | (1u64 << 32))
}

/// Execute the `select` command.
fn exec_select(scene: &mut VizServer, entity_type: &str, id: u32) -> serde_json::Value {
    use popkern_core::topology::{EdgeHandle, FaceHandle, VertexHandle};

    match entity_type {
        "face" => {
            let handle = FaceHandle::from(key_data(id));
            scene.select(Selection::Face(handle));
            ok_response(json!({"selected": "face", "id": id}))
        }
        "edge" => {
            let handle = EdgeHandle::from(key_data(id));
            scene.select(Selection::Edge(handle));
            ok_response(json!({"selected": "edge", "id": id}))
        }
        "vertex" => {
            let handle = VertexHandle::from(key_data(id));
            scene.select(Selection::Vertex(handle));
            ok_response(json!({"selected": "vertex", "id": id}))
        }
        other => err_response(&format!("unknown entity type: {other}")),
    }
}

/// Execute the `camera` command.
fn exec_camera(scene: &mut VizServer, args: &CameraArgs) -> serde_json::Value {
    match args {
        CameraArgs::Preset(name) => {
            let preset = match name.as_str() {
                "front" => Some(ViewPreset::Front),
                "back" => Some(ViewPreset::Back),
                "top" => Some(ViewPreset::Top),
                "bottom" => Some(ViewPreset::Bottom),
                "left" => Some(ViewPreset::Left),
                "right" => Some(ViewPreset::Right),
                "iso" => Some(ViewPreset::Iso),
                _ => None,
            };
            match preset {
                Some(p) => {
                    scene.camera.set_view_preset(p);
                    ok_response(json!({"camera": name}))
                }
                None => err_response(&format!("unknown camera preset: {name}")),
            }
        }
        CameraArgs::Explicit { eye, target, up } => {
            scene.camera.eye = glam::Vec3::new(eye[0], eye[1], eye[2]);
            scene.camera.target = glam::Vec3::new(target[0], target[1], target[2]);
            scene.camera.up = glam::Vec3::new(up[0], up[1], up[2]);
            ok_response(json!({
                "eye": eye,
                "target": target,
                "up": up,
            }))
        }
    }
}

/// Execute the `screenshot` command (stub -- headless rendering wired in Task 11).
fn exec_screenshot(path: &str) -> serde_json::Value {
    ok_response(json!({"screenshot": path, "status": "queued"}))
}

/// Execute the `inspect` command.
fn exec_inspect(scene: &VizServer, entity_type: &str, id: u32) -> serde_json::Value {
    let Some(brep) = &scene.brep else {
        return err_response("no model loaded");
    };

    match entity_type {
        "face" => {
            let handle = popkern_core::topology::FaceHandle::from(key_data(id));
            match brep.faces.get(handle) {
                Some(face) => {
                    let surface_type = face
                        .surface
                        .as_ref()
                        .map(|s| format!("{:?}", s.geometry))
                        .unwrap_or_else(|| "none".into());
                    ok_response(json!({
                        "entity": "face",
                        "id": id,
                        "surface": surface_type,
                    }))
                }
                None => err_response(&format!("face {id} not found")),
            }
        }
        "edge" => {
            let handle = popkern_core::topology::EdgeHandle::from(key_data(id));
            match brep.edges.get(handle) {
                Some(_edge) => ok_response(json!({
                    "entity": "edge",
                    "id": id,
                })),
                None => err_response(&format!("edge {id} not found")),
            }
        }
        "vertex" => {
            let handle = popkern_core::topology::VertexHandle::from(key_data(id));
            match brep.vertices.get(handle) {
                Some(v) => ok_response(json!({
                    "entity": "vertex",
                    "id": id,
                    "position": [v.position.x, v.position.y, v.position.z],
                })),
                None => err_response(&format!("vertex {id} not found")),
            }
        }
        other => err_response(&format!("unknown entity type: {other}")),
    }
}

/// Execute the `highlight` command (visual effect only -- stub for now).
fn exec_highlight(entity_type: &str, ids: &[u32], color: &str) -> serde_json::Value {
    ok_response(json!({
        "highlighted": entity_type,
        "ids": ids,
        "color": color,
    }))
}

/// Execute the `validate` command.
fn exec_validate(scene: &VizServer) -> serde_json::Value {
    let Some(brep) = &scene.brep else {
        return err_response("no model loaded");
    };

    let errors = brep.validate_all();
    if errors.is_empty() {
        ok_response(json!({"valid": true, "errors": []}))
    } else {
        let msgs: Vec<String> = errors.iter().map(|e| e.to_string()).collect();
        ok_response(json!({"valid": false, "errors": msgs}))
    }
}

/// Execute the `fit` command.
fn exec_fit(scene: &mut VizServer) -> serde_json::Value {
    match scene.model_bounds() {
        Some((min, max)) => {
            scene.camera.fit_to_bounds(min, max);
            ok_response(json!({"fit": true}))
        }
        None => err_response("no model loaded to fit"),
    }
}

/// Execute the `wireframe` command.
fn exec_wireframe(scene: &mut VizServer, on: bool) -> serde_json::Value {
    scene.display.show_wireframe = on;
    ok_response(json!({"wireframe": on}))
}

/// Execute the `colors` command.
fn exec_colors(scene: &mut VizServer, mode: &str) -> serde_json::Value {
    match mode {
        "surface" => {
            scene.display.color_mode = ColorMode::Surface;
            ok_response(json!({"color_mode": "surface"}))
        }
        "flat" => {
            scene.display.color_mode = ColorMode::Flat;
            ok_response(json!({"color_mode": "flat"}))
        }
        "none" => {
            scene.display.color_mode = ColorMode::None;
            ok_response(json!({"color_mode": "none"}))
        }
        other => err_response(&format!("unknown color mode: {other}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exec_wireframe_on() {
        let mut scene = VizServer::new();
        let result = execute_command(&mut scene, &Command::Wireframe { on: true });
        assert_eq!(result["ok"], true);
        assert!(scene.display.show_wireframe);
    }

    #[test]
    fn exec_wireframe_off() {
        let mut scene = VizServer::new();
        let result = execute_command(&mut scene, &Command::Wireframe { on: false });
        assert_eq!(result["ok"], true);
        assert!(!scene.display.show_wireframe);
    }

    #[test]
    fn exec_colors_surface() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Colors {
                mode: "surface".into(),
            },
        );
        assert_eq!(result["ok"], true);
        assert_eq!(scene.display.color_mode, ColorMode::Surface);
    }

    #[test]
    fn exec_colors_flat() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Colors {
                mode: "flat".into(),
            },
        );
        assert_eq!(result["ok"], true);
        assert_eq!(scene.display.color_mode, ColorMode::Flat);
    }

    #[test]
    fn exec_colors_none() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Colors {
                mode: "none".into(),
            },
        );
        assert_eq!(result["ok"], true);
        assert_eq!(scene.display.color_mode, ColorMode::None);
    }

    #[test]
    fn exec_colors_unknown_mode() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Colors {
                mode: "rainbow".into(),
            },
        );
        assert_eq!(result["ok"], false);
    }

    #[test]
    fn exec_quit() {
        let mut scene = VizServer::new();
        let result = execute_command(&mut scene, &Command::Quit);
        assert_eq!(result["ok"], true);
    }

    #[test]
    fn exec_fit_no_model() {
        let mut scene = VizServer::new();
        let result = execute_command(&mut scene, &Command::Fit);
        assert_eq!(result["ok"], false);
    }

    #[test]
    fn exec_validate_no_model() {
        let mut scene = VizServer::new();
        let result = execute_command(&mut scene, &Command::Validate);
        assert_eq!(result["ok"], false);
    }

    #[test]
    fn exec_inspect_no_model() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Inspect {
                entity_type: "face".into(),
                id: 0,
            },
        );
        assert_eq!(result["ok"], false);
    }

    #[test]
    fn exec_camera_preset_front() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Camera {
                args: CameraArgs::Preset("front".into()),
            },
        );
        assert_eq!(result["ok"], true);
    }

    #[test]
    fn exec_camera_explicit() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Camera {
                args: CameraArgs::Explicit {
                    eye: [0.0, 0.0, 5.0],
                    target: [0.0, 0.0, 0.0],
                    up: [0.0, 1.0, 0.0],
                },
            },
        );
        assert_eq!(result["ok"], true);
        assert!((scene.camera.eye.z - 5.0).abs() < f32::EPSILON);
    }

    #[test]
    fn exec_camera_unknown_preset() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Camera {
                args: CameraArgs::Preset("diagonal".into()),
            },
        );
        assert_eq!(result["ok"], false);
    }

    #[test]
    fn exec_select_face() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Select {
                entity_type: "face".into(),
                id: 0,
            },
        );
        assert_eq!(result["ok"], true);
        assert!(matches!(scene.selection, Selection::Face(_)));
    }

    #[test]
    fn exec_select_unknown_type() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Select {
                entity_type: "body".into(),
                id: 0,
            },
        );
        assert_eq!(result["ok"], false);
    }

    #[test]
    fn exec_screenshot_queued() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Screenshot {
                path: "/tmp/test.png".into(),
            },
        );
        assert_eq!(result["ok"], true);
        assert_eq!(result["data"]["status"], "queued");
    }

    #[test]
    fn exec_highlight_returns_ok() {
        let mut scene = VizServer::new();
        let result = execute_command(
            &mut scene,
            &Command::Highlight {
                entity_type: "face".into(),
                ids: vec![1, 2],
                color: "red".into(),
            },
        );
        assert_eq!(result["ok"], true);
    }

    #[test]
    fn exec_validate_with_model() {
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let mut scene = VizServer::new();
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        scene.load_brep(brep).expect("load succeeds");

        let result = execute_command(&mut scene, &Command::Validate);
        assert_eq!(result["ok"], true);
        assert!(result["data"]["valid"].as_bool().is_some());
    }

    #[test]
    fn exec_fit_with_model() {
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let mut scene = VizServer::new();
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        scene.load_brep(brep).expect("load succeeds");

        let result = execute_command(&mut scene, &Command::Fit);
        assert_eq!(result["ok"], true);
    }
}
