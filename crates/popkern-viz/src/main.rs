//! CLI entry point for popkern-viz.

fn main() {
    env_logger::init();

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: popkern-viz <command> [options]");
        eprintln!("Commands:");
        eprintln!("  view <file>        Open interactive viewer");
        eprintln!("  screenshot <file>  Render to PNG (headless)");
        std::process::exit(1);
    }

    match args[1].as_str() {
        "view" => {
            eprintln!("Interactive viewer not yet implemented");
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
