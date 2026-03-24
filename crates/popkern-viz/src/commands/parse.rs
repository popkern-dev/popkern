//! Command parsing: text line to [`Command`] conversion.

use crate::error::VizError;

use super::{CameraArgs, Command};

/// Parse a single line of text into a [`Command`].
///
/// Returns a [`VizError`] if the line cannot be parsed.
pub fn parse_command(line: &str) -> Result<Command, VizError> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.is_empty() {
        return Err(VizError::Load("empty command".into()));
    }

    match parts[0] {
        "select" => parse_select(&parts),
        "camera" => parse_camera(&parts),
        "screenshot" => parse_screenshot(&parts),
        "inspect" => parse_inspect(&parts),
        "highlight" => parse_highlight(&parts),
        "validate" => Ok(Command::Validate),
        "fit" => Ok(Command::Fit),
        "wireframe" => parse_wireframe(&parts),
        "colors" => parse_colors(&parts),
        "quit" => Ok(Command::Quit),
        other => Err(VizError::Load(format!("unknown command: {other}"))),
    }
}

/// Parse `select <entity_type> <id>`.
fn parse_select(parts: &[&str]) -> Result<Command, VizError> {
    if parts.len() != 3 {
        return Err(VizError::Load(
            "select requires: select <entity_type> <id>".into(),
        ));
    }
    let id: u32 = parts[2]
        .parse()
        .map_err(|_| VizError::Load(format!("invalid id: {}", parts[2])))?;
    Ok(Command::Select {
        entity_type: parts[1].to_string(),
        id,
    })
}

/// Parse `camera <preset>` or `camera <9 floats>`.
fn parse_camera(parts: &[&str]) -> Result<Command, VizError> {
    if parts.len() == 2 {
        return Ok(Command::Camera {
            args: CameraArgs::Preset(parts[1].to_string()),
        });
    }
    if parts.len() == 10 {
        let floats: Vec<f32> = parts[1..]
            .iter()
            .map(|s| {
                s.parse::<f32>()
                    .map_err(|_| VizError::Load(format!("invalid float: {s}")))
            })
            .collect::<Result<Vec<_>, _>>()?;
        return Ok(Command::Camera {
            args: CameraArgs::Explicit {
                eye: [floats[0], floats[1], floats[2]],
                target: [floats[3], floats[4], floats[5]],
                up: [floats[6], floats[7], floats[8]],
            },
        });
    }
    Err(VizError::Load(
        "camera requires: camera <preset> or camera <eye x y z> <target x y z> <up x y z>".into(),
    ))
}

/// Parse `screenshot <path>`.
fn parse_screenshot(parts: &[&str]) -> Result<Command, VizError> {
    if parts.len() != 2 {
        return Err(VizError::Load(
            "screenshot requires: screenshot <path>".into(),
        ));
    }
    Ok(Command::Screenshot {
        path: parts[1].to_string(),
    })
}

/// Parse `inspect <entity_type> <id>`.
fn parse_inspect(parts: &[&str]) -> Result<Command, VizError> {
    if parts.len() != 3 {
        return Err(VizError::Load(
            "inspect requires: inspect <entity_type> <id>".into(),
        ));
    }
    let id: u32 = parts[2]
        .parse()
        .map_err(|_| VizError::Load(format!("invalid id: {}", parts[2])))?;
    Ok(Command::Inspect {
        entity_type: parts[1].to_string(),
        id,
    })
}

/// Parse `highlight <entity_type> <id1,id2,...> <color>`.
fn parse_highlight(parts: &[&str]) -> Result<Command, VizError> {
    if parts.len() != 4 {
        return Err(VizError::Load(
            "highlight requires: highlight <entity_type> <id1,id2,...> <color>".into(),
        ));
    }
    let ids: Vec<u32> = parts[2]
        .split(',')
        .map(|s| {
            s.parse::<u32>()
                .map_err(|_| VizError::Load(format!("invalid id: {s}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Command::Highlight {
        entity_type: parts[1].to_string(),
        ids,
        color: parts[3].to_string(),
    })
}

/// Parse `wireframe on|off`.
fn parse_wireframe(parts: &[&str]) -> Result<Command, VizError> {
    if parts.len() != 2 {
        return Err(VizError::Load(
            "wireframe requires: wireframe on|off".into(),
        ));
    }
    match parts[1] {
        "on" => Ok(Command::Wireframe { on: true }),
        "off" => Ok(Command::Wireframe { on: false }),
        other => Err(VizError::Load(format!(
            "wireframe: expected on|off, got {other}"
        ))),
    }
}

/// Parse `colors <mode>`.
fn parse_colors(parts: &[&str]) -> Result<Command, VizError> {
    if parts.len() != 2 {
        return Err(VizError::Load("colors requires: colors <mode>".into()));
    }
    Ok(Command::Colors {
        mode: parts[1].to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_select_face() {
        let cmd = parse_command("select face 3").unwrap();
        assert_eq!(
            cmd,
            Command::Select {
                entity_type: "face".into(),
                id: 3
            }
        );
    }

    #[test]
    fn parse_select_edge() {
        let cmd = parse_command("select edge 7").unwrap();
        assert_eq!(
            cmd,
            Command::Select {
                entity_type: "edge".into(),
                id: 7
            }
        );
    }

    #[test]
    fn parse_select_vertex() {
        let cmd = parse_command("select vertex 0").unwrap();
        assert_eq!(
            cmd,
            Command::Select {
                entity_type: "vertex".into(),
                id: 0
            }
        );
    }

    #[test]
    fn parse_camera_preset() {
        let cmd = parse_command("camera front").unwrap();
        assert_eq!(
            cmd,
            Command::Camera {
                args: CameraArgs::Preset("front".into())
            }
        );
    }

    #[test]
    fn parse_camera_explicit() {
        let cmd = parse_command("camera 0 0 5 0 0 0 0 1 0").unwrap();
        assert_eq!(
            cmd,
            Command::Camera {
                args: CameraArgs::Explicit {
                    eye: [0.0, 0.0, 5.0],
                    target: [0.0, 0.0, 0.0],
                    up: [0.0, 1.0, 0.0],
                }
            }
        );
    }

    #[test]
    fn parse_screenshot() {
        let cmd = parse_command("screenshot /tmp/view.png").unwrap();
        assert_eq!(
            cmd,
            Command::Screenshot {
                path: "/tmp/view.png".into()
            }
        );
    }

    #[test]
    fn parse_inspect() {
        let cmd = parse_command("inspect face 3").unwrap();
        assert_eq!(
            cmd,
            Command::Inspect {
                entity_type: "face".into(),
                id: 3
            }
        );
    }

    #[test]
    fn parse_highlight() {
        let cmd = parse_command("highlight face 1,2,3 red").unwrap();
        assert_eq!(
            cmd,
            Command::Highlight {
                entity_type: "face".into(),
                ids: vec![1, 2, 3],
                color: "red".into(),
            }
        );
    }

    #[test]
    fn parse_validate() {
        assert_eq!(parse_command("validate").unwrap(), Command::Validate);
    }

    #[test]
    fn parse_fit() {
        assert_eq!(parse_command("fit").unwrap(), Command::Fit);
    }

    #[test]
    fn parse_wireframe_on() {
        assert_eq!(
            parse_command("wireframe on").unwrap(),
            Command::Wireframe { on: true }
        );
    }

    #[test]
    fn parse_wireframe_off() {
        assert_eq!(
            parse_command("wireframe off").unwrap(),
            Command::Wireframe { on: false }
        );
    }

    #[test]
    fn parse_colors() {
        assert_eq!(
            parse_command("colors surface").unwrap(),
            Command::Colors {
                mode: "surface".into()
            }
        );
    }

    #[test]
    fn parse_quit() {
        assert_eq!(parse_command("quit").unwrap(), Command::Quit);
    }

    #[test]
    fn parse_unknown_command_is_error() {
        assert!(parse_command("foo").is_err());
    }

    #[test]
    fn parse_empty_is_error() {
        assert!(parse_command("").is_err());
        assert!(parse_command("   ").is_err());
    }

    #[test]
    fn parse_select_bad_id_is_error() {
        assert!(parse_command("select face abc").is_err());
    }

    #[test]
    fn parse_camera_wrong_arg_count_is_error() {
        assert!(parse_command("camera 0 0 5 0").is_err());
        assert!(parse_command("camera").is_err());
    }

    #[test]
    fn parse_wireframe_bad_arg_is_error() {
        assert!(parse_command("wireframe maybe").is_err());
    }
}
