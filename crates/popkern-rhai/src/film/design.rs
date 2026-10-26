//! `s.evolve`: a generative design converging (PK-959). The design's
//! layer (`"design <label>"`) shows the design's recorded fields as
//! display meshes, never bodies; `show`, `fade` and `highlight` reach it
//! by that selector like any body.

use rhai::{Dynamic, EvalAltResult, NativeCallContext};

use super::actions::options_or_unit;
use super::api::ScriptStep;
use super::model::{Channel, Ease};
use crate::args::{self, Options};
use crate::error::argument_error;

const EVOLVE_OPTIONS: &[&str] = &["from", "to", "at", "duration", "ease"];

/// An update count: a number in `0..=updates`, or `"end"` (the final
/// design).
fn update(argument: &str, value: &Dynamic, updates: usize) -> Result<f64, Box<EvalAltResult>> {
    if let Some(text) = value.clone().try_cast::<rhai::ImmutableString>() {
        return match text.as_str() {
            "end" => Ok(updates as f64),
            other => Err(argument_error(
                argument,
                format!("an update number or \"end\", got \"{other}\""),
            )),
        };
    }
    let n = args::float(argument, value)?;
    if !(n.is_finite() && (0.0..=updates as f64).contains(&n)) {
        return Err(argument_error(
            argument,
            format!(
                "an update between 0 and {updates} (the design's solve ran {updates}), got {n}"
            ),
        ));
    }
    Ok(n)
}

/// `s.evolve(design, #{ from, to, at, duration, ease })`: the design
/// `design` (its feature label, or `"design <label>"`) runs from update
/// `from` (default 0, the first recorded) to `to` (default `"end"`).
/// The time spreads evenly over the recorded designs, which are denser
/// early, where the design changes most; the ease defaults to linear.
pub(super) fn evolve(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    design: Dynamic,
    options: Option<Dynamic>,
) -> Result<(), Box<EvalAltResult>> {
    let context = s.context();
    let name = args::text("design", &design)?;
    let selector = if name.starts_with("design ") {
        name.clone()
    } else {
        format!("design {name}")
    };
    let Some(found) = context.design(&selector) else {
        let known: Vec<&str> = context.designs.iter().map(|d| d.label.as_str()).collect();
        let reason = if known.is_empty() {
            "the document has no generative design to evolve".to_string()
        } else {
            format!(
                "the document has no generative design '{name}'; it has {}",
                known.join(", ")
            )
        };
        return Err(argument_error("design", reason));
    };
    let updates = found.updates;
    let options = Options::parse("evolve", &options_or_unit(options), EVOLVE_OPTIONS)?;
    let from = options
        .get("from")
        .map(|v| update("from", v, updates))
        .transpose()?
        .unwrap_or(0.0);
    let to = options
        .get("to")
        .map(|v| update("to", v, updates))
        .transpose()?
        .unwrap_or(updates as f64);
    let (start, duration, ease) = s.timing(&options)?;
    let ease = if options.get("ease").is_some() {
        ease
    } else {
        Ease::Linear
    };
    s.push(
        &ctx,
        (start, duration, ease),
        Channel::Evolve {
            target: selector,
            from,
            to,
        },
    );
    Ok(())
}
