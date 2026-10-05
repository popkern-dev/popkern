//! `popkern run model.rhai` — the minimal live loop (design §10 phase
//! 3e): run a script under budgets, print the replay report, export
//! STEP / PNG / `.pkn`, and optionally rebuild on save. No server, no
//! protocol, no viewer (those are the daemon's, phase 6).
//!
//! Exit codes distinguish what went wrong: `0` served; `2` the script
//! itself (usage, syntax, argument errors, budget exceeded); `3` a
//! kernel decline (a feature failed); `4` a gate or export failure
//! (not watertight, validation, I/O); `5` a malformed project folder;
//! `6` a saved document whose applied history diverged, under
//! `--strict`.
//!
//! `popkern build <dir>` builds a project folder of migration scripts
//! (live-projects-and-history.md); `popkern new <dir>` scaffolds one.

#![forbid(unsafe_code)]

pub mod args;
pub mod project;
pub mod render;
pub mod report;
pub mod run;
pub mod watch;

pub use args::{Cli, Command, parse_args};
pub use report::{Outcome, RunSummary};
pub use run::{CliError, execute, exit_code};
pub use watch::Watcher;

use std::process::ExitCode;

/// The whole command: parse, execute (or watch), print, exit.
pub fn main_with(args: &[String]) -> ExitCode {
    let cli = match parse_args(args) {
        Ok(cli) => cli,
        Err(error) => {
            eprintln!("popkern: {error}\n\n{}", args::USAGE);
            return ExitCode::from(2);
        }
    };
    match cli.command {
        Command::Help => {
            println!("{}", args::USAGE);
            ExitCode::SUCCESS
        }
        Command::Build => ExitCode::from(project::build_command(&cli)),
        Command::New => ExitCode::from(project::new_command(&cli)),
        Command::Run => {
            if cli.watch {
                return watch::watch_loop(&cli);
            }
            let summary = execute(&cli);
            report::print(&cli, &summary);
            ExitCode::from(exit_code(&summary))
        }
    }
}
