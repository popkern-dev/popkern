//! `popkern` binary entry point.

/// The static musl release builds' allocator (see Cargo.toml).
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> std::process::ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    popkern_cli::main_with(&args)
}
