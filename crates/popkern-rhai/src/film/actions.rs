//! The film actions a step's closure calls (docs/architecture/
//! presentation-films.md §4), split from `api.rs` under the 500-line
//! budget. Each validates its names against the document snapshot and
//! records one track at its script line.

use rhai::{Dynamic, EvalAltResult, NativeCallContext};

use super::api::{ScriptStep, ease, length3, seconds, target, unit3};
use super::context::FilmContext;
use super::model::{CameraTo, Channel, Ease, SectionKey};
use crate::args::{self, Options};
use crate::error::argument_error;
use crate::profile_api::ScriptAxis;

const TIMING: &[&str] = &["at", "duration", "ease"];
const MOVE_OPTIONS: &[&str] = &["from", "to", "at", "duration", "ease"];
const CAMERA_ACTION: &[&str] = &[
    "azimuth",
    "elevation",
    "zoom",
    "target",
    "at",
    "duration",
    "ease",
];
const SECTION_OPTIONS: &[&str] = &[
    "from", "to", "from_at", "to_at", "combine", "at", "duration", "ease",
];

pub(super) fn options_or_unit(value: Option<Dynamic>) -> Dynamic {
    value.unwrap_or(Dynamic::UNIT)
}

pub(super) fn orbit(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    degrees: Dynamic,
    options: Option<Dynamic>,
) -> Result<(), Box<EvalAltResult>> {
    let degrees = args::float("degrees", &degrees)?;
    let options = Options::parse("orbit", &options_or_unit(options), TIMING)?;
    let timing = s.timing(&options)?;
    s.push(
        &ctx,
        (
            timing.0,
            timing.1,
            options.with("ease", ease)?.unwrap_or(Ease::Linear),
        ),
        Channel::Orbit { degrees },
    );
    Ok(())
}

pub(super) fn camera(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    options: Dynamic,
) -> Result<(), Box<EvalAltResult>> {
    let options = Options::parse("camera", &options, CAMERA_ACTION)?;
    let context = s.context();
    let to = CameraTo {
        azimuth: options.with("azimuth", args::float)?,
        elevation: options.with("elevation", args::float)?,
        zoom: options.with("zoom", args::float)?,
        target: options
            .get("target")
            .map(|v| target(&context, "target", v))
            .transpose()?,
    };
    if to.zoom.is_some_and(|z| !(z.is_finite() && z > 0.0)) {
        return Err(argument_error("zoom", "a positive factor"));
    }
    if to.elevation.is_some_and(|e| !(-89.0..=89.0).contains(&e)) {
        return Err(argument_error("elevation", "degrees in -89..89"));
    }
    let timing = s.timing(&options)?;
    s.push(&ctx, timing, Channel::Camera { to });
    Ok(())
}

pub(super) fn movement(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    path: Dynamic,
    options: Dynamic,
) -> Result<(), Box<EvalAltResult>> {
    let context = s.context();
    let target = target(&context, "target", &path)?;
    let options = Options::parse("move", &options, MOVE_OPTIONS)?;
    let from = options.with("from", |a, v| length3(&context, a, v))?;
    let to = options.with("to", |a, v| length3(&context, a, v))?;
    let (from, to) = match (from, to) {
        (Some(f), None) => (f, [0.0; 3]),
        (None, Some(t)) => ([0.0; 3], t),
        (Some(f), Some(t)) => (f, t),
        (None, None) => {
            return Err(argument_error(
                "options",
                "move needs from: [x, y, z] (fly in) or to: [x, y, z] (fly out)",
            ));
        }
    };
    let timing = s.timing(&options)?;
    s.push(&ctx, timing, Channel::Move { target, from, to });
    Ok(())
}

pub(super) fn turn(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    path: Dynamic,
    axis: ScriptAxis,
    from: Dynamic,
    to: Dynamic,
    options: Option<Dynamic>,
) -> Result<(), Box<EvalAltResult>> {
    let context = s.context();
    let target = target(&context, "target", &path)?;
    let ScriptAxis::Literal { origin, direction } = axis else {
        return Err(argument_error(
            "axis",
            "a film turns about a literal axis_at(point, direction)",
        ));
    };
    let n =
        (direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2])
            .sqrt();
    if !(n.is_finite() && n > 1e-12) {
        return Err(argument_error("axis", "a non-zero direction"));
    }
    let direction = direction.map(|c| c / n);
    let options = Options::parse("turn", &options_or_unit(options), TIMING)?;
    let timing = s.timing(&options)?;
    s.push(
        &ctx,
        timing,
        Channel::Turn {
            target,
            origin: origin.map(|c| c * context.unit_scale),
            direction,
            from: args::float("from", &from)?,
            to: args::float("to", &to)?,
        },
    );
    Ok(())
}

pub(super) fn joint(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    param: &str,
    from: Dynamic,
    to: Dynamic,
    options: Option<Dynamic>,
) -> Result<(), Box<EvalAltResult>> {
    let context = s.context();
    if !context.joints.iter().any(|j| j == param) {
        let reason = if context.kind == "part" {
            "a joint is an assembly parameter that only moves occurrences; a part document has none"
                .to_string()
        } else if context.params.iter().any(|p| p == param) {
            format!(
                "'{param}' drives a part's geometry, so sweeping it would rebuild every frame; a joint only moves occurrences"
            )
        } else if context.joints.is_empty() {
            "the assembly declares no joint parameters".to_string()
        } else {
            format!(
                "the assembly has no joint '{param}'; its joints are {}",
                context.joints.join(", ")
            )
        };
        return Err(argument_error("param", reason));
    }
    let options = Options::parse("joint", &options_or_unit(options), TIMING)?;
    let timing = s.timing(&options)?;
    s.push(
        &ctx,
        timing,
        Channel::Joint {
            param: param.to_string(),
            from: args::float("from", &from)?,
            to: args::float("to", &to)?,
        },
    );
    Ok(())
}

pub(super) fn explode(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    from: Dynamic,
    to: Dynamic,
    options: Option<Dynamic>,
) -> Result<(), Box<EvalAltResult>> {
    let (from, to) = (args::float("from", &from)?, args::float("to", &to)?);
    if !(from.is_finite() && to.is_finite() && from >= 0.0 && to >= 0.0) {
        return Err(argument_error(
            "from",
            "explode factors are 0 (assembled) or more",
        ));
    }
    if s.context().kind != "assembly" {
        return Err(argument_error(
            "explode",
            "explode pulls an assembly's parts apart; a part document has none",
        ));
    }
    let options = Options::parse("explode", &options_or_unit(options), TIMING)?;
    let timing = s.timing(&options)?;
    s.push(&ctx, timing, Channel::Explode { from, to });
    Ok(())
}

pub(super) fn section_key(
    context: &FilmContext,
    argument: &str,
    value: &Dynamic,
    at: Option<f64>,
) -> Result<SectionKey, Box<EvalAltResult>> {
    if value.is_unit() {
        return Ok(SectionKey::Off);
    }
    if let Some(text) = value.clone().try_cast::<rhai::ImmutableString>() {
        return match text.as_str() {
            "x" | "y" | "z" => {
                let at = at.unwrap_or(0.5);
                if !(0.0..=1.0).contains(&at) {
                    return Err(argument_error(
                        argument,
                        format!("a position 0..1 across the bounds, got {at}"),
                    ));
                }
                Ok(SectionKey::Axis {
                    axis: text.to_string(),
                    at,
                })
            }
            "off" => Ok(SectionKey::Off),
            other => Err(argument_error(
                argument,
                format!("\"x\", \"y\", \"z\", \"off\" or #{{ normal, offset }}, got \"{other}\""),
            )),
        };
    }
    let options = Options::parse("section plane", value, &["normal", "offset"])?;
    let normal = options
        .with("normal", unit3)?
        .ok_or_else(|| argument_error(argument, "a plane is #{ normal: [x, y, z], offset: d }"))?;
    let offset = options.with("offset", args::float)?.unwrap_or(0.0) * context.unit_scale;
    Ok(SectionKey::Plane { normal, offset })
}

pub(super) fn section(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    options: Dynamic,
) -> Result<(), Box<EvalAltResult>> {
    let context = s.context();
    let options = Options::parse("section", &options, SECTION_OPTIONS)?;
    let from_at = options.with("from_at", args::float)?;
    let to_at = options.with("to_at", args::float)?;
    let from = match options.get("from") {
        Some(v) => section_key(&context, "from", v, from_at)?,
        None => SectionKey::Off,
    };
    // No `to`: the same plane, moved to `to_at` when that is given.
    let to = match (options.get("to"), &from, to_at) {
        (Some(v), _, _) => section_key(&context, "to", v, to_at)?,
        (None, SectionKey::Axis { axis, .. }, Some(at)) => {
            section_key(&context, "to_at", &Dynamic::from(axis.clone()), Some(at))?
        }
        (None, _, _) => from.clone(),
    };
    let combine = match options.with("combine", args::text)? {
        Some(c) if c == "any" || c == "all" => c,
        Some(c) => {
            return Err(argument_error(
                "combine",
                format!("\"any\" or \"all\", got \"{c}\""),
            ));
        }
        None => "any".into(),
    };
    let timing = s.timing(&options)?;
    s.push(&ctx, timing, Channel::Section { from, to, combine });
    Ok(())
}

pub(super) fn fade(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    path: Dynamic,
    from: Dynamic,
    to: Dynamic,
    options: Option<Dynamic>,
) -> Result<(), Box<EvalAltResult>> {
    let context = s.context();
    let target = target(&context, "target", &path)?;
    let (from, to) = (args::float("from", &from)?, args::float("to", &to)?);
    if !((0.0..=1.0).contains(&from) && (0.0..=1.0).contains(&to)) {
        return Err(argument_error("from", "opacities are between 0 and 1"));
    }
    let options = Options::parse("fade", &options_or_unit(options), TIMING)?;
    let timing = s.timing(&options)?;
    s.push(&ctx, timing, Channel::Fade { target, from, to });
    Ok(())
}

pub(super) fn show(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    path: Dynamic,
    visible: bool,
    options: Option<Dynamic>,
) -> Result<(), Box<EvalAltResult>> {
    let context = s.context();
    let target = target(&context, "target", &path)?;
    let options = Options::parse("show", &options_or_unit(options), &["at"])?;
    let at = options.with("at", seconds)?.unwrap_or(0.0);
    s.push(
        &ctx,
        (s.start + at, 0.0, Ease::Linear),
        Channel::Show { target, visible },
    );
    Ok(())
}

pub(super) fn highlight(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    path: Dynamic,
    options: Option<Dynamic>,
) -> Result<(), Box<EvalAltResult>> {
    let context = s.context();
    let target = target(&context, "target", &path)?;
    let options = Options::parse("highlight", &options_or_unit(options), TIMING)?;
    let timing = s.timing(&options)?;
    s.push(&ctx, timing, Channel::Highlight { target });
    Ok(())
}
