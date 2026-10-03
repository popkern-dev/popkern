//! Assembly parameters in the Rhai skin (PK-638):
//!
//! ```rhai
//! a.param("width", 0.4, "m");
//! a.drive("width", [sun, planet, ring]);        // each part's own "width"
//! a.drive("bore", carrier, "output_bore_radius"); // a differently named part parameter
//! ```
//!
//! Setting `width` (a slider, `--param width=…`, `a.set_param`) sets the
//! driven parameter in every listed part and rebuilds each once.

use popkern_history::unit::Unit;
use rhai::{Array, Dynamic, Engine, EvalAltResult};

use crate::args;
use crate::assembly_api::{IntoAssemblyResult, ScriptAssembly, ScriptDefinition};
use crate::error::argument_error;

fn param(
    a: &mut ScriptAssembly,
    name: &str,
    value: Dynamic,
    unit: Dynamic,
) -> Result<(), Box<EvalAltResult>> {
    let unit = if unit.is_unit() {
        Unit::Meters
    } else {
        args::unit("unit", &unit)?
    };
    let value = args::float("value", &value)?;
    a.borrow_mut()
        .assembly
        .add_param(name, value, unit)
        .assembly()
}

fn drive_each(a: &mut ScriptAssembly, name: &str, parts: Array) -> Result<(), Box<EvalAltResult>> {
    for part in parts {
        let def = part
            .try_cast::<ScriptDefinition>()
            .ok_or_else(|| argument_error("parts", "expected an array of part definitions"))?;
        a.borrow_mut()
            .assembly
            .drive(name, def.id, name)
            .assembly()?;
    }
    Ok(())
}

/// Register the parameter verbs.
pub fn register(engine: &mut Engine) {
    engine
        .register_fn(
            "param",
            |a: &mut ScriptAssembly, name: &str, value: Dynamic| {
                param(a, name, value, Dynamic::UNIT)
            },
        )
        .register_fn("param", param)
        .register_fn("drive", drive_each)
        .register_fn(
            "drive",
            |a: &mut ScriptAssembly, name: &str, def: ScriptDefinition, part_param: &str| {
                a.borrow_mut()
                    .assembly
                    .drive(name, def.id, part_param)
                    .assembly()
            },
        )
        .register_fn(
            "set_param",
            |a: &mut ScriptAssembly, name: &str, value: Dynamic| {
                let value = args::float("value", &value)?;
                a.borrow_mut()
                    .assembly
                    .set_assembly_param(name, value)
                    .assembly()
                    .map(|_| ())
            },
        );
}
