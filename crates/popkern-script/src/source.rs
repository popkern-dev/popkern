//! Where in a script each feature was recorded (live-projects-and-history.md
//! §5.6): a host that tracks positions sets the script being run and,
//! as evaluation steps, the current call position; the recorder notes
//! that position against every feature it records. Session-local
//! provenance for the history scrubber and the editor — never saved.

use std::cell::RefCell;

use serde::Serialize;

/// A position in a named script.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SourceSpan {
    /// The script's name (`003_rounds.rhai`).
    pub script: String,
    /// 1-based line.
    pub line: u32,
    /// 1-based column.
    pub column: u32,
}

thread_local! {
    static SCRIPT: RefCell<Option<String>> = const { RefCell::new(None) };
    static POSITION: RefCell<Option<(u32, u32)>> = const { RefCell::new(None) };
}

/// Begin tracking `script` (`None` stops tracking).
pub fn set_script(script: Option<&str>) {
    SCRIPT.with(|s| *s.borrow_mut() = script.map(str::to_string));
    POSITION.with(|p| *p.borrow_mut() = None);
}

/// The evaluator reached `line`:`column` of the tracked script.
pub fn set_position(line: u32, column: u32) {
    POSITION.with(|p| *p.borrow_mut() = Some((line, column)));
}

/// The current position, when a script is tracked.
pub fn current() -> Option<SourceSpan> {
    let script = SCRIPT.with(|s| s.borrow().clone())?;
    let (line, column) = POSITION.with(|p| *p.borrow())?;
    Some(SourceSpan {
        script,
        line,
        column,
    })
}
