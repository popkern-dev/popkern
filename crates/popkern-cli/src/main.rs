//! `popkern` binary entry point.

fn main() -> std::process::ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    popkern_cli::main_with(&args)
}
