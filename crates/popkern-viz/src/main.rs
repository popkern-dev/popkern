//! CLI entry point for popkern-viz.

use popkern_viz::app::VizApp;

/// Print comprehensive CLI help text.
fn print_help() {
    eprintln!("popkern-viz — BRep debug viewer");
    eprintln!();
    eprintln!("USAGE:");
    eprintln!("    popkern-viz view <file>              Open interactive viewer");
    eprintln!("    popkern-viz screenshot <file>        Render to PNG (headless)");
    eprintln!("        --output <path>                  Output PNG path (default: output.png)");
    eprintln!("        --camera <preset|x,y,z>          Camera position (default: iso)");
    eprintln!("        --width <N>                      Image width (default: 1280)");
    eprintln!("        --height <N>                     Image height (default: 720)");
    eprintln!("        --no-wireframe                   Disable wireframe overlay");
    eprintln!();
    eprintln!("SUPPORTED FORMATS:");
    eprintln!("    .json       BRep JSON (full topology)");
    eprintln!("    .obj        Mesh only (no topology inspection)");
    eprintln!("    .stl        Mesh only (no topology inspection)");
    eprintln!("    .step/.stp  STEP (requires popkern-io, coming soon)");
    eprintln!();
    eprintln!("AGENT PROTOCOL:");
    eprintln!("    In interactive mode, send commands on stdin for programmatic control.");
    eprintln!("    Commands: select, camera, screenshot, inspect, highlight, validate,");
    eprintln!("             fit, wireframe, colors, quit");
    eprintln!("    Responses: one JSON object per line on stdout.");
}

fn main() {
    env_logger::init();

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 || args[1] == "--help" || args[1] == "-h" {
        print_help();
        std::process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    match args[1].as_str() {
        "view" => {
            let mut app = VizApp::new();

            // If a path is given, try to load it; otherwise show a unit cube.
            let brep = if let Some(file_arg) = args.get(2) {
                let path = std::path::Path::new(file_arg);
                match popkern_viz::loader::load_model(path) {
                    Ok(b) => b,
                    Err(e) => {
                        eprintln!("Failed to load {file_arg}: {e}");
                        eprintln!("Falling back to unit cube.");
                        popkern_core::topology::BRep::unit_cube(
                            popkern_core::tolerance::Tolerance::default(),
                        )
                        .expect("unit cube should build")
                    }
                }
            } else {
                popkern_core::topology::BRep::unit_cube(
                    popkern_core::tolerance::Tolerance::default(),
                )
                .expect("unit cube should build")
            };
            app.scene
                .load_brep(brep)
                .expect("tessellation should succeed");

            let native_options = eframe::NativeOptions {
                viewport: egui::ViewportBuilder::default()
                    .with_inner_size([1280.0, 720.0])
                    .with_title("popkern-viz"),
                depth_buffer: 32,
                ..Default::default()
            };

            eframe::run_native(
                "popkern-viz",
                native_options,
                Box::new(|_cc| Ok(Box::new(app))),
            )
            .unwrap();
        }
        "screenshot" => {
            let mut output_path = "output.png".to_string();
            let mut width = 1280u32;
            let mut height = 720u32;
            let mut camera_preset = "iso".to_string();
            let mut show_wireframe = true;

            // Parse optional CLI flags: --output, --width, --height, --camera, --no-wireframe
            let mut i = 2;
            while i < args.len() {
                match args[i].as_str() {
                    "--output" | "-o" => {
                        i += 1;
                        if let Some(val) = args.get(i) {
                            output_path = val.clone();
                        }
                    }
                    "--width" => {
                        i += 1;
                        if let Some(val) = args.get(i) {
                            width = val.parse().expect("--width must be a number");
                        }
                    }
                    "--height" => {
                        i += 1;
                        if let Some(val) = args.get(i) {
                            height = val.parse().expect("--height must be a number");
                        }
                    }
                    "--camera" => {
                        i += 1;
                        if let Some(val) = args.get(i) {
                            camera_preset = val.clone();
                        }
                    }
                    "--no-wireframe" => {
                        show_wireframe = false;
                    }
                    other => {
                        // Treat bare argument as output path
                        output_path = other.to_string();
                    }
                }
                i += 1;
            }

            let brep = popkern_core::topology::BRep::unit_cube(
                popkern_core::tolerance::Tolerance::default(),
            )
            .expect("unit cube should build");

            let mut camera = popkern_viz::camera::Camera::default_iso();
            // Apply camera preset
            let preset = match camera_preset.as_str() {
                "front" => Some(popkern_viz::camera::ViewPreset::Front),
                "back" => Some(popkern_viz::camera::ViewPreset::Back),
                "top" => Some(popkern_viz::camera::ViewPreset::Top),
                "bottom" => Some(popkern_viz::camera::ViewPreset::Bottom),
                "left" => Some(popkern_viz::camera::ViewPreset::Left),
                "right" => Some(popkern_viz::camera::ViewPreset::Right),
                "iso" => Some(popkern_viz::camera::ViewPreset::Iso),
                _ => {
                    eprintln!("Unknown camera preset: {camera_preset}, using iso");
                    Some(popkern_viz::camera::ViewPreset::Iso)
                }
            };
            if let Some(p) = preset {
                camera.set_view_preset(p);
            }

            match popkern_viz::headless::render_to_png(
                &brep,
                &camera,
                width,
                height,
                std::path::Path::new(&output_path),
                show_wireframe,
            ) {
                Ok(()) => eprintln!("Screenshot saved to {output_path}"),
                Err(e) => {
                    eprintln!("Headless render failed: {e}");
                    std::process::exit(1);
                }
            }
        }
        _ => {
            eprintln!("Unknown command: {}", args[1]);
            std::process::exit(1);
        }
    }
}
