//! Script positions for recorded features (PK-661): with
//! [`crate::HostOptions::track_sources`] the engine's debugger steps
//! every expression and publishes its position to
//! [`popkern_script::source`], where the recorder notes it against
//! each feature it records. The innermost expression evaluated before
//! a recording call is that call's own argument or receiver, so the
//! noted line is the line of the call.

use rhai::Engine;
use rhai::debugger::DebuggerCommand;

/// Step every expression and publish its position.
#[allow(deprecated)] // rhai marks the debugger API volatile, not deprecated.
pub(crate) fn enable(engine: &mut Engine) {
    engine.register_debugger(
        // The default debugger starts at `Start`; answering every event
        // with `StepInto` then stops at every expression.
        |_, debugger| debugger,
        |_, _, _, _, position| {
            if let (Some(line), Some(column)) = (position.line(), position.position()) {
                popkern_script::source::set_position(line as u32, column as u32);
            }
            Ok(DebuggerCommand::StepInto)
        },
    );
}
