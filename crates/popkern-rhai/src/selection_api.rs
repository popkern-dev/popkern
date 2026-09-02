//! `edges()` / `faces()` query builders (spec §6.1) with every
//! predicate spelled as in the grammar; enum arguments are strings.

use popkern_script::{EdgeSelector, FaceSelector, Scalar};
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args;

fn edge_between(
    sel: &mut EdgeSelector,
    a: Dynamic,
    b: Dynamic,
) -> Result<EdgeSelector, Box<EvalAltResult>> {
    Ok(sel
        .clone()
        .between(args::surface_kind("a", &a)?, args::surface_kind("b", &b)?))
}

fn edge_at_height(sel: &mut EdgeSelector, z: Dynamic) -> Result<EdgeSelector, Box<EvalAltResult>> {
    Ok(sel.clone().at_height(args::float("z", &z)?))
}

fn edge_direction(sel: &mut EdgeSelector, d: Dynamic) -> Result<EdgeSelector, Box<EvalAltResult>> {
    Ok(sel.clone().direction(args::dir("direction", &d)?))
}

fn edge_chord(
    sel: &mut EdgeSelector,
    min: Dynamic,
    max: Dynamic,
) -> Result<EdgeSelector, Box<EvalAltResult>> {
    Ok(sel
        .clone()
        .chord_length_between(args::float("min", &min)?, args::float("max", &max)?))
}

fn face_kind(sel: &mut FaceSelector, k: Dynamic) -> Result<FaceSelector, Box<EvalAltResult>> {
    Ok(sel.clone().kind(args::surface_kind("kind", &k)?))
}

fn face_facing(sel: &mut FaceSelector, d: Dynamic) -> Result<FaceSelector, Box<EvalAltResult>> {
    Ok(sel.clone().facing(args::dir("direction", &d)?))
}

fn face_facing_neg(sel: &mut FaceSelector, d: Dynamic) -> Result<FaceSelector, Box<EvalAltResult>> {
    Ok(sel.clone().facing_neg(args::dir("direction", &d)?))
}

fn face_at_height(sel: &mut FaceSelector, z: Dynamic) -> Result<FaceSelector, Box<EvalAltResult>> {
    Ok(sel.clone().at_height(args::float("z", &z)?))
}

/// Register the selection grammar and the scalar-binding helpers.
pub fn register(engine: &mut Engine) {
    engine
        .register_type_with_name::<EdgeSelector>("EdgeSelector")
        .register_fn("edges", popkern_script::selection::edges)
        .register_fn("circular", |s: &mut EdgeSelector| s.clone().circular())
        .register_fn("straight", |s: &mut EdgeSelector| s.clone().straight())
        .register_fn("conic", |s: &mut EdgeSelector| s.clone().conic())
        .register_fn("nurbs", |s: &mut EdgeSelector| s.clone().nurbs())
        .register_fn("between", edge_between)
        .register_fn("at_height", edge_at_height)
        .register_fn("direction", edge_direction)
        .register_fn("chord_length_between", edge_chord)
        .register_fn("not_seam", |s: &mut EdgeSelector| s.clone().not_seam())
        .register_fn("labeled", |s: &mut EdgeSelector, l: &str| {
            s.clone().labeled(l)
        })
        .register_fn("named", |s: &mut EdgeSelector, t: &str| {
            s.clone().spelled(t)
        })
        .register_fn("or", |s: &mut EdgeSelector, other: EdgeSelector| {
            s.clone().or(other)
        })
        .register_fn("allow_empty", |s: &mut EdgeSelector| {
            s.clone().allow_empty()
        })
        .register_fn("describe", |s: &mut EdgeSelector| s.describe())
        .register_fn("to_string", |s: &mut EdgeSelector| s.describe());
    engine
        .register_type_with_name::<FaceSelector>("FaceSelector")
        .register_fn("faces", popkern_script::selection::faces)
        .register_fn("kind", face_kind)
        .register_fn("facing", face_facing)
        .register_fn("facing_neg", face_facing_neg)
        .register_fn("at_height", face_at_height)
        .register_fn("labeled", |s: &mut FaceSelector, l: &str| {
            s.clone().labeled(l)
        })
        .register_fn("named", |s: &mut FaceSelector, t: &str| {
            s.clone().spelled(t)
        })
        .register_fn("allow_empty", |s: &mut FaceSelector| {
            s.clone().allow_empty()
        })
        .register_fn("describe", |s: &mut FaceSelector| s.describe())
        .register_fn("to_string", |s: &mut FaceSelector| s.describe());
    engine
        .register_type_with_name::<Scalar>("Scalar")
        .register_fn("p", |name: &str| Scalar::Param(name.to_string()))
        .register_fn("expr", |text: &str| Scalar::Expr(text.to_string()))
        .register_fn("to_string", |s: &mut Scalar| match s {
            Scalar::Literal(v) => v.to_string(),
            Scalar::Param(name) => format!("p({name:?})"),
            Scalar::Expr(text) => format!("expr({text:?})"),
        });
}
