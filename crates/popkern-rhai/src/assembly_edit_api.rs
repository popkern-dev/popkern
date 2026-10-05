//! `a.edit_part(name, |m| …)` (PK-678): edit a part definition of an
//! assembly — typically in a project migration that continues it
//! (`let a = previous();`) — through the same model verbs a part
//! script uses (`m.insert_after(…)`, `m.feature(…)`, fillets, …).
//!
//! ```rhai
//! let a = previous();
//! a.edit_part("pin", |m| {
//!     m.insert_after(m.feature("extrude-1"), |m, body| body.chamfer(edges().at_height(1.0), 0.02));
//! });
//! a
//! ```
//!
//! Every instance of the part follows; mates are solved again on the
//! edited geometry. Transactional: a closure that raises, a part left
//! without a body, or mates that no longer hold leave the assembly as
//! it was.

use rhai::{Dynamic, Engine, EvalAltResult, FnPtr, NativeCallContext};

use crate::assembly_api::{IntoAssemblyResult, ScriptAssembly, in_part};
use crate::error::argument_error;
use crate::model_api::ScriptModel;

fn edit_part(
    ctx: NativeCallContext,
    a: &mut ScriptAssembly,
    name: &str,
    edit: FnPtr,
) -> Result<(), Box<EvalAltResult>> {
    let (id, model) = {
        let mut state = a.borrow_mut();
        let id = state
            .assembly
            .definition_named(name)
            .ok_or_else(|| argument_error("part", format!("no part '{name}' in this assembly")))?;
        let model = state.assembly.part_model_mut(id).assembly()?.clone();
        (id, model)
    };
    let handle = ScriptModel::new(model);
    let outcome: Result<Dynamic, _> = edit.call_within_context(&ctx, (handle.clone(),));
    let _: Dynamic = outcome.map_err(|error| in_part(name, error))?;
    a.borrow_mut()
        .assembly
        .replace_part_model(id, handle.into_inner())
        .assembly()
}

/// Register `edit_part`.
pub fn register(engine: &mut Engine) {
    engine.register_fn("edit_part", edit_part);
}
