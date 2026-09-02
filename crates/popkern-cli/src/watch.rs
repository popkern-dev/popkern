//! Rebuild on save: poll the script's modification time, debounce,
//! rerun, keep the last good outputs. Deliberately bounded — no
//! server, no protocol, no viewer (design §10 phase 3e).

use std::process::ExitCode;
use std::time::{Duration, SystemTime};

use popkern_core::time::Stopwatch;

use crate::args::Cli;
use crate::report::{self, Outcome, RunSummary};
use crate::run::execute;

/// Poll interval.
pub const POLL: Duration = Duration::from_millis(200);
/// A change must be quiet this long before it triggers a rebuild.
pub const DEBOUNCE: Duration = Duration::from_millis(150);

/// A file watcher a host polls; testable without threads.
pub struct Watcher {
    cli: Cli,
    last_modified: Option<SystemTime>,
    pending_since: Option<Stopwatch>,
    /// Summaries of every run so far, newest last.
    pub history: Vec<RunSummary>,
}

impl Watcher {
    /// Start watching `cli.script`; the first `check` runs the script.
    pub fn new(cli: Cli) -> Watcher {
        Watcher {
            cli,
            last_modified: None,
            pending_since: None,
            history: Vec::new(),
        }
    }

    fn modified(&self) -> Option<SystemTime> {
        std::fs::metadata(&self.cli.script)
            .and_then(|m| m.modified())
            .ok()
    }

    /// Whether the newest run served.
    pub fn last_served(&self) -> bool {
        matches!(
            self.history.last().map(|s| &s.outcome),
            Some(Outcome::Served)
        )
    }

    /// Poll once: run the script when it changed and the change has
    /// settled. Returns the new summary when a run happened.
    pub fn check(&mut self) -> Option<&RunSummary> {
        let modified = self.modified();
        let changed = modified != self.last_modified;
        if changed {
            self.last_modified = modified;
            self.pending_since = Some(Stopwatch::start());
            return None;
        }
        // Without a monotonic clock (never the case natively) a change
        // is due at the next poll.
        let due = self
            .pending_since
            .as_ref()
            .is_some_and(|since| since.elapsed().is_none_or(|e| e >= DEBOUNCE));
        if !due {
            return None;
        }
        self.pending_since = None;
        let summary = execute(&self.cli);
        self.history.push(summary);
        self.history.last()
    }
}

/// The blocking watch loop (Ctrl-C ends the process).
pub fn watch_loop(cli: &Cli) -> ExitCode {
    let mut watcher = Watcher::new(cli.clone());
    eprintln!(
        "popkern: watching {} (Ctrl-C to stop)",
        cli.script.display()
    );
    loop {
        if let Some(summary) = watcher.check() {
            report::print(cli, summary);
        }
        std::thread::sleep(POLL);
    }
}
