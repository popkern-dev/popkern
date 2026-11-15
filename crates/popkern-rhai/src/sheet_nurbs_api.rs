//! `m.sheet_nurbs(points[, options])`: the sheet of a B-spline surface
//! from its control net (`points[i][j]`, rows along u), with the options
//! `degree` (one integer for both directions, or `[du, dv]`; default 3),
//! `weights` (the net's shape; default all one), `knots` (`[ku, kv]`,
//! each clamped; default clamped uniform), `label` and the placement
//! `at`, `axis`, `ref`.

use popkern_script::Placement;
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::IntoScriptResult;
use crate::model_api::ScriptModel;

const OPTIONS: &[&str] = &["label", "at", "axis", "ref", "degree", "weights", "knots"];

fn rows<T>(
    argument: &str,
    value: &Dynamic,
    item: impl Fn(&Dynamic) -> Result<T, Box<EvalAltResult>>,
) -> Result<Vec<Vec<T>>, Box<EvalAltResult>> {
    let outer = value
        .clone()
        .into_array()
        .map_err(|_| crate::error::argument_error(argument, "a list of rows"))?;
    outer
        .iter()
        .map(|row| {
            row.clone()
                .into_array()
                .map_err(|_| crate::error::argument_error(argument, "each row a list"))?
                .iter()
                .map(&item)
                .collect()
        })
        .collect()
}

fn floats(argument: &str, value: &Dynamic) -> Result<Vec<f64>, Box<EvalAltResult>> {
    value
        .clone()
        .into_array()
        .map_err(|_| crate::error::argument_error(argument, "a list of numbers"))?
        .iter()
        .map(|v| args::float(argument, v))
        .collect()
}

fn degree(value: &Dynamic) -> Result<[u32; 2], Box<EvalAltResult>> {
    let one = |v: &Dynamic| {
        v.as_int()
            .ok()
            .and_then(|n| u32::try_from(n).ok())
            .filter(|n| *n >= 1)
            .ok_or_else(|| crate::error::argument_error("degree", "an integer of at least 1"))
    };
    if let Ok(pair) = value.clone().into_array() {
        if pair.len() != 2 {
            return Err(crate::error::argument_error("degree", "[du, dv]"));
        }
        return Ok([one(&pair[0])?, one(&pair[1])?]);
    }
    let d = one(value)?;
    Ok([d, d])
}

fn sheet_nurbs(
    m: &mut ScriptModel,
    p: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("sheet_nurbs", o, OPTIONS)?;
    let points = rows("points", p, |v| args::vec3("points", v))?;
    let weights: Vec<Vec<f64>> = options
        .with("weights", |a, v| rows(a, v, |w| args::float(a, w)))?
        .unwrap_or_default();
    let degree = options.with("degree", |_, v| degree(v))?.unwrap_or([3, 3]);
    let knots: [Vec<f64>; 2] = options
        .with("knots", |a, v| {
            let pair = v
                .clone()
                .into_array()
                .map_err(|_| crate::error::argument_error(a, "[ku, kv]"))?;
            if pair.len() != 2 {
                return Err(crate::error::argument_error(a, "[ku, kv]"));
            }
            Ok([floats(a, &pair[0])?, floats(a, &pair[1])?])
        })?
        .unwrap_or_default();
    let mut placement = Placement::default();
    if let Some(at) = options.with("at", args::vec3)? {
        placement.at = at;
    }
    if let Some(axis) = options.with("axis", args::vec3)? {
        placement.axis = axis;
    }
    if let Some(reference) = options.with("ref", args::vec3)? {
        placement.reference = reference;
    }
    let body = m.borrow_mut().sheet_nurbs(
        &points,
        &weights,
        degree,
        knots,
        placement,
        &options.label()?,
    );
    Ok(ScriptBody::new(m.clone(), body.script()?))
}

/// Register `sheet_nurbs`.
pub fn register(engine: &mut Engine) {
    let none = Dynamic::UNIT;
    engine.register_fn("sheet_nurbs", move |m: &mut ScriptModel, p: Dynamic| {
        sheet_nurbs(m, &p, &none)
    });
    engine.register_fn(
        "sheet_nurbs",
        |m: &mut ScriptModel, p: Dynamic, o: Dynamic| sheet_nurbs(m, &p, &o),
    );
}
