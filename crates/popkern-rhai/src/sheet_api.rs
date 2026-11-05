//! Sheet verbs (PK-1083): `m.sheet_rect(x, y)`, `m.sheet_disc(r)`,
//! `m.sheet_polygon(points)`, each with `label` and the placement `at`,
//! `axis`, `ref` (`centred` for a rectangle), and
//! `m.sheet_fill(wire[, edges])`, the face a closed planar wire bounds.

use popkern_script::{Placement, Scalar};
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::IntoScriptResult;
use crate::model_api::ScriptModel;

const PLACED: &[&str] = &["label", "at", "axis", "ref"];
const RECT_OPTIONS: &[&str] = &["label", "at", "axis", "ref", "centred"];
const LABEL_ONLY: &[&str] = &["label"];

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

fn served(
    m: &ScriptModel,
    body: Result<popkern_script::Body, popkern_script::ScriptError>,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    Ok(ScriptBody::new(m.clone(), body.script()?))
}

fn rect(
    m: &mut ScriptModel,
    x: &Dynamic,
    y: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("sheet_rect", o, RECT_OPTIONS)?;
    let size: [Scalar; 2] = [args::scalar("x", x)?, args::scalar("y", y)?];
    let centred = options.with("centred", args::boolean)?.unwrap_or(false);
    let body = m
        .borrow_mut()
        .sheet_rect(size, centred, placement(&options)?, &options.label()?);
    served(m, body)
}

fn disc(m: &mut ScriptModel, r: &Dynamic, o: &Dynamic) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("sheet_disc", o, PLACED)?;
    let r = args::scalar("radius", r)?;
    let body = m
        .borrow_mut()
        .sheet_disc(r, placement(&options)?, &options.label()?);
    served(m, body)
}

fn polygon(
    m: &mut ScriptModel,
    p: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("sheet_polygon", o, PLACED)?;
    let list = p
        .clone()
        .into_array()
        .map_err(|_| crate::error::argument_error("points", "a list of [x, y] points"))?;
    let points = list
        .iter()
        .map(|v| {
            let q = args::point("points", v, 2)?;
            Ok([q[0], q[1]])
        })
        .collect::<Result<Vec<[f64; 2]>, Box<EvalAltResult>>>()?;
    let body = m
        .borrow_mut()
        .sheet_polygon(&points, placement(&options)?, &options.label()?);
    served(m, body)
}

fn fill(
    m: &mut ScriptModel,
    wire: ScriptBody,
    edges: popkern_script::selection::EdgeSelector,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("sheet_fill", o, LABEL_ONLY)?;
    let body = m
        .borrow_mut()
        .sheet_fill(wire.body, &edges, &options.label()?);
    served(m, body)
}

/// Register the sheet verbs.
pub fn register(engine: &mut Engine) {
    let none = Dynamic::UNIT;
    let n = none.clone();
    engine.register_fn(
        "sheet_rect",
        move |m: &mut ScriptModel, x: Dynamic, y: Dynamic| rect(m, &x, &y, &n),
    );
    engine.register_fn(
        "sheet_rect",
        |m: &mut ScriptModel, x: Dynamic, y: Dynamic, o: Dynamic| rect(m, &x, &y, &o),
    );
    let n = none.clone();
    engine.register_fn("sheet_disc", move |m: &mut ScriptModel, r: Dynamic| {
        disc(m, &r, &n)
    });
    engine.register_fn(
        "sheet_disc",
        |m: &mut ScriptModel, r: Dynamic, o: Dynamic| disc(m, &r, &o),
    );
    let n = none.clone();
    engine.register_fn("sheet_polygon", move |m: &mut ScriptModel, p: Dynamic| {
        polygon(m, &p, &n)
    });
    engine.register_fn(
        "sheet_polygon",
        |m: &mut ScriptModel, p: Dynamic, o: Dynamic| polygon(m, &p, &o),
    );
    let n = none.clone();
    engine.register_fn("sheet_fill", move |m: &mut ScriptModel, w: ScriptBody| {
        fill(m, w, popkern_script::selection::edges(), &n)
    });
    let n = none;
    engine.register_fn(
        "sheet_fill",
        move |m: &mut ScriptModel, w: ScriptBody, e: popkern_script::selection::EdgeSelector| {
            fill(m, w, e, &n)
        },
    );
    engine.register_fn(
        "sheet_fill",
        |m: &mut ScriptModel,
         w: ScriptBody,
         e: popkern_script::selection::EdgeSelector,
         o: Dynamic| fill(m, w, e, &o),
    );
}
