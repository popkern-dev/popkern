//! Wire verbs (PK-1083): `m.wire_line(a, b)`, `m.wire_polyline(points)`,
//! `m.wire_arc(start, via, end)`, `m.wire_circle(r)`,
//! `m.wire_helix(r, pitch, turns)` and `m.wire_spline(points)`, each a
//! wire body with an options map (`label`, the placement `at`, `axis`,
//! `ref`; `closed` for a polyline, `left_hand` for a helix, `degree`
//! and `control` for a spline). A wire is a path: `b.sweep(wire)` and `wire.pipe(r)`
//! take its one edge without an edge query.

use popkern_script::{Placement, Scalar};
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::IntoScriptResult;
use crate::model_api::ScriptModel;

const PLACED: &[&str] = &["label", "at", "axis", "ref"];
const POLYLINE_OPTIONS: &[&str] = &["label", "at", "axis", "ref", "closed"];
const HELIX_OPTIONS: &[&str] = &["label", "at", "axis", "ref", "left_hand"];
const SPLINE_OPTIONS: &[&str] = &["label", "at", "axis", "ref", "degree", "control"];

fn placement(options: &Options) -> Result<Placement, Box<EvalAltResult>> {
    let mut p = Placement::default();
    if let Some(at) = options.with("at", args::vec3)? {
        p.at = at;
    }
    if let Some(axis) = options.with("axis", args::vec3)? {
        p.axis = axis;
    }
    if let Some(reference) = options.with("ref", args::vec3)? {
        p.reference = reference;
    }
    Ok(p)
}

fn points(argument: &str, value: &Dynamic) -> Result<Vec<[f64; 3]>, Box<EvalAltResult>> {
    let list = value
        .clone()
        .into_array()
        .map_err(|_| crate::error::argument_error(argument, "a list of [x, y, z] points"))?;
    list.iter().map(|p| args::vec3(argument, p)).collect()
}

fn served(
    m: &ScriptModel,
    body: Result<popkern_script::Body, popkern_script::ScriptError>,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    Ok(ScriptBody::new(m.clone(), body.script()?))
}

fn line(
    m: &mut ScriptModel,
    a: &Dynamic,
    b: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("wire_line", o, PLACED)?;
    let (a, b) = (args::vec3("from", a)?, args::vec3("to", b)?);
    let body = m
        .borrow_mut()
        .wire_line(a, b, placement(&options)?, &options.label()?);
    served(m, body)
}

fn polyline(
    m: &mut ScriptModel,
    p: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("wire_polyline", o, POLYLINE_OPTIONS)?;
    let pts = points("points", p)?;
    let closed = options.with("closed", args::boolean)?.unwrap_or(false);
    let body = m
        .borrow_mut()
        .wire_polyline(&pts, closed, placement(&options)?, &options.label()?);
    served(m, body)
}

fn arc(
    m: &mut ScriptModel,
    [a, v, b]: [&Dynamic; 3],
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("wire_arc", o, PLACED)?;
    let p = [
        args::vec3("start", a)?,
        args::vec3("via", v)?,
        args::vec3("end", b)?,
    ];
    let body = m
        .borrow_mut()
        .wire_arc(p, placement(&options)?, &options.label()?);
    served(m, body)
}

fn circle(m: &mut ScriptModel, r: &Dynamic, o: &Dynamic) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("wire_circle", o, PLACED)?;
    let r = args::scalar("radius", r)?;
    let body = m
        .borrow_mut()
        .wire_circle(r, placement(&options)?, &options.label()?);
    served(m, body)
}

fn helix(
    m: &mut ScriptModel,
    [r, p, t]: [&Dynamic; 3],
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("wire_helix", o, HELIX_OPTIONS)?;
    let (r, p, t): (Scalar, Scalar, Scalar) = (
        args::scalar("radius", r)?,
        args::scalar("pitch", p)?,
        args::scalar("turns", t)?,
    );
    let left = options.with("left_hand", args::boolean)?.unwrap_or(false);
    let body = m
        .borrow_mut()
        .wire_helix(r, p, t, left, placement(&options)?, &options.label()?);
    served(m, body)
}

fn spline(m: &mut ScriptModel, p: &Dynamic, o: &Dynamic) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("wire_spline", o, SPLINE_OPTIONS)?;
    let pts = points("points", p)?;
    let degree = options
        .with("degree", |a, v| {
            v.as_int()
                .ok()
                .and_then(|d| u32::try_from(d).ok())
                .filter(|d| *d >= 1)
                .ok_or_else(|| crate::error::argument_error(a, "a positive integer"))
        })?
        .unwrap_or(3);
    let control = options.with("control", args::boolean)?.unwrap_or(false);
    let body = m.borrow_mut().wire_spline(
        &pts,
        degree,
        control,
        placement(&options)?,
        &options.label()?,
    );
    served(m, body)
}

/// Register the wire verbs and the edge-query-free sweep and pipe.
pub fn register(engine: &mut Engine) {
    let none = Dynamic::UNIT;
    let n = none.clone();
    engine.register_fn(
        "wire_line",
        move |m: &mut ScriptModel, a: Dynamic, b: Dynamic| line(m, &a, &b, &n),
    );
    engine.register_fn(
        "wire_line",
        |m: &mut ScriptModel, a: Dynamic, b: Dynamic, o: Dynamic| line(m, &a, &b, &o),
    );
    let n = none.clone();
    engine.register_fn("wire_polyline", move |m: &mut ScriptModel, p: Dynamic| {
        polyline(m, &p, &n)
    });
    engine.register_fn(
        "wire_polyline",
        |m: &mut ScriptModel, p: Dynamic, o: Dynamic| polyline(m, &p, &o),
    );
    let n = none.clone();
    engine.register_fn(
        "wire_arc",
        move |m: &mut ScriptModel, a: Dynamic, v: Dynamic, b: Dynamic| arc(m, [&a, &v, &b], &n),
    );
    engine.register_fn(
        "wire_arc",
        |m: &mut ScriptModel, a: Dynamic, v: Dynamic, b: Dynamic, o: Dynamic| {
            arc(m, [&a, &v, &b], &o)
        },
    );
    let n = none.clone();
    engine.register_fn("wire_circle", move |m: &mut ScriptModel, r: Dynamic| {
        circle(m, &r, &n)
    });
    engine.register_fn(
        "wire_circle",
        |m: &mut ScriptModel, r: Dynamic, o: Dynamic| circle(m, &r, &o),
    );
    let n = none.clone();
    engine.register_fn(
        "wire_helix",
        move |m: &mut ScriptModel, r: Dynamic, p: Dynamic, t: Dynamic| helix(m, [&r, &p, &t], &n),
    );
    engine.register_fn(
        "wire_helix",
        |m: &mut ScriptModel, r: Dynamic, p: Dynamic, t: Dynamic, o: Dynamic| {
            helix(m, [&r, &p, &t], &o)
        },
    );
    let n = none;
    engine.register_fn("wire_spline", move |m: &mut ScriptModel, p: Dynamic| {
        spline(m, &p, &n)
    });
    engine.register_fn(
        "wire_spline",
        |m: &mut ScriptModel, p: Dynamic, o: Dynamic| spline(m, &p, &o),
    );
    // A wire's one edge is the path without a query.
    engine.register_fn("sweep", |b: &mut ScriptBody, path: ScriptBody| {
        sweep_all(b, path, rhai::Map::new())
    });
    engine.register_fn(
        "sweep",
        |b: &mut ScriptBody, path: ScriptBody, o: rhai::Map| sweep_all(b, path, o),
    );
    engine.register_fn("pipe", |b: &mut ScriptBody, r: Dynamic| {
        pipe_all(b, r, rhai::Map::new())
    });
    engine.register_fn("pipe", |b: &mut ScriptBody, r: Dynamic, o: rhai::Map| {
        pipe_all(b, r, o)
    });
}

/// A sweep's `corner` option: `"mitre"` (the default) or `"round"`.
pub(crate) fn sweep_corner(
    argument: &str,
    value: &Dynamic,
) -> Result<popkern_history::feature_sweep_path::SweepCorner, Box<EvalAltResult>> {
    use popkern_history::feature_sweep_path::SweepCorner;
    match args::text(argument, value)?.as_str() {
        "mitre" | "miter" => Ok(SweepCorner::Mitre),
        "round" => Ok(SweepCorner::Round),
        other => Err(crate::error::argument_error(
            argument,
            format!("\"mitre\" or \"round\", not \"{other}\""),
        )),
    }
}

fn sweep_all(
    b: &mut ScriptBody,
    path: ScriptBody,
    options: rhai::Map,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    b.sweep_along(
        path,
        popkern_script::selection::edges(),
        &Dynamic::from_map(options),
    )
}

fn pipe_all(
    b: &mut ScriptBody,
    radius: Dynamic,
    options: rhai::Map,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    b.pipe_along(
        popkern_script::selection::edges(),
        radius,
        &Dynamic::from_map(options),
    )
}
