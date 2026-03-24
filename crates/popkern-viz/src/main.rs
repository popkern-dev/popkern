//! CLI entry point for popkern-viz.

use popkern_viz::app::VizApp;

fn main() {
    env_logger::init();

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: popkern-viz <command> [options]");
        eprintln!("Commands:");
        eprintln!("  view [file]        Open interactive viewer");
        eprintln!("  screenshot <file>  Render to PNG (headless)");
        std::process::exit(1);
    }

    match args[1].as_str() {
        "view" => {
            let mut app = VizApp::new();

            // If a path is given, load a model; otherwise show a unit cube.
            if args.get(2).is_some() {
                // File loading will be implemented in Task 12.
                // For now, always load a unit cube as demo geometry.
                let brep = popkern_core::topology::BRep::unit_cube(
                    popkern_core::tolerance::Tolerance::default(),
                )
                .expect("unit cube should build");
                app.scene
                    .load_brep(brep)
                    .expect("tessellation should succeed");
            } else {
                // No file — show a unit cube by default.
                let brep = popkern_core::topology::BRep::unit_cube(
                    popkern_core::tolerance::Tolerance::default(),
                )
                .expect("unit cube should build");
                app.scene
                    .load_brep(brep)
                    .expect("tessellation should succeed");
            }

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
            eprintln!("Headless rendering not yet implemented");
        }
        _ => {
            eprintln!("Unknown command: {}", args[1]);
            std::process::exit(1);
        }
    }
}
