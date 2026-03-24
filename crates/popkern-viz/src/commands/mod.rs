//! Stdin/stdout JSON command protocol for agent control.
//!
//! Parses line-based text commands from stdin, dispatches them against
//! a [`VizServer`], and returns JSON responses on stdout.

mod exec;
mod parse;

use crate::error::VizError;
use crate::scene::VizServer;

pub use exec::execute_command;
pub use parse::parse_command;

// ---------------------------------------------------------------------------
// Command types
// ---------------------------------------------------------------------------

/// Camera arguments: either a named preset or explicit eye/target/up vectors.
#[derive(Debug, Clone, PartialEq)]
pub enum CameraArgs {
    /// A named preset such as "front", "top", "iso", etc.
    Preset(String),
    /// Explicit camera parameters.
    Explicit {
        /// Camera position in world space.
        eye: [f32; 3],
        /// Point the camera looks at.
        target: [f32; 3],
        /// Up direction.
        up: [f32; 3],
    },
}

/// A command that can be dispatched to the visualization server.
#[derive(Debug, Clone, PartialEq)]
pub enum Command {
    /// Select a BRep entity by type and ID.
    Select {
        /// Entity type: "face", "edge", or "vertex".
        entity_type: String,
        /// Entity ID (slot-map index).
        id: u32,
    },
    /// Set the camera to a preset or explicit position.
    Camera {
        /// Camera arguments.
        args: CameraArgs,
    },
    /// Render a screenshot to a file.
    Screenshot {
        /// Output file path.
        path: String,
    },
    /// Inspect a BRep entity and return its metadata.
    Inspect {
        /// Entity type: "face", "edge", or "vertex".
        entity_type: String,
        /// Entity ID (slot-map index).
        id: u32,
    },
    /// Highlight a set of entities with a given color.
    Highlight {
        /// Entity type: "face", "edge", or "vertex".
        entity_type: String,
        /// Entity IDs to highlight.
        ids: Vec<u32>,
        /// CSS-style color string (e.g. "red", "#ff0000").
        color: String,
    },
    /// Run BRep validation and return results.
    Validate,
    /// Fit the camera to the model bounds.
    Fit,
    /// Toggle wireframe overlay on or off.
    Wireframe {
        /// Whether wireframe is enabled.
        on: bool,
    },
    /// Set the face coloring mode.
    Colors {
        /// Color mode: "surface", "flat", or "none".
        mode: String,
    },
    /// Quit the application.
    Quit,
}

/// Read commands from a [`std::io::BufRead`] source and execute them against the scene.
///
/// Writes JSON responses to the given writer, one per line.
/// Returns when a `quit` command is received or the input ends.
pub fn run_command_loop<R: std::io::BufRead, W: std::io::Write>(
    scene: &mut VizServer,
    reader: R,
    writer: &mut W,
) -> Result<(), VizError> {
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        match parse_command(&line) {
            Ok(cmd) => {
                let is_quit = cmd == Command::Quit;
                let response = execute_command(scene, &cmd);
                serde_json::to_writer(&mut *writer, &response)?;
                writeln!(writer)?;
                writer.flush()?;
                if is_quit {
                    break;
                }
            }
            Err(e) => {
                let response = exec::err_response(&e.to_string());
                serde_json::to_writer(&mut *writer, &response)?;
                writeln!(writer)?;
                writer.flush()?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn command_loop_processes_and_quits() {
        let input = "wireframe off\nquit\n";
        let mut output = Vec::new();
        let mut scene = VizServer::new();

        run_command_loop(&mut scene, input.as_bytes(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.trim().lines().collect();
        assert_eq!(lines.len(), 2);

        let r1: serde_json::Value = serde_json::from_str(lines[0]).unwrap();
        assert_eq!(r1["ok"], true);

        let r2: serde_json::Value = serde_json::from_str(lines[1]).unwrap();
        assert_eq!(r2["ok"], true);
    }

    #[test]
    fn command_loop_handles_errors() {
        let input = "badcommand\nquit\n";
        let mut output = Vec::new();
        let mut scene = VizServer::new();

        run_command_loop(&mut scene, input.as_bytes(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.trim().lines().collect();
        assert_eq!(lines.len(), 2);

        let r1: serde_json::Value = serde_json::from_str(lines[0]).unwrap();
        assert_eq!(r1["ok"], false);
    }

    #[test]
    fn command_loop_skips_empty_lines() {
        let input = "\n  \nquit\n";
        let mut output = Vec::new();
        let mut scene = VizServer::new();

        run_command_loop(&mut scene, input.as_bytes(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.trim().lines().collect();
        assert_eq!(lines.len(), 1); // only quit response
    }
}
