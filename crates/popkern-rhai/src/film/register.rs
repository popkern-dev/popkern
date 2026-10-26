//! Registering the film surface on an engine (split from `api.rs`
//! under the 500-line budget): `document()` and its getters, the
//! storyboard and step builders, and every action.

use std::rc::Rc;

use rhai::{Dynamic, Engine, EvalAltResult, FnPtr, Map, NativeCallContext};

use super::actions::{
    camera, explode, fade, highlight, joint, movement, orbit, section, show, turn,
};
use super::api::{ScriptFilm, ScriptFilmDocument, ScriptStep, add_step, storyboard};
use super::context::FilmContext;
use crate::error::argument_error;
use crate::profile_api::ScriptAxis;

fn document_paths(d: &mut ScriptFilmDocument) -> rhai::Array {
    d.0.paths.iter().map(|p| Dynamic::from(p.clone())).collect()
}

fn document_joints(d: &mut ScriptFilmDocument) -> rhai::Array {
    d.0.joints
        .iter()
        .map(|p| Dynamic::from(p.clone()))
        .collect()
}

/// Register the film surface on `engine` for a film named `name`
/// reading `context`.
pub fn register(engine: &mut Engine, context: Rc<FilmContext>, name: &str) {
    engine
        .register_type_with_name::<ScriptFilm>("Film")
        .register_type_with_name::<ScriptStep>("FilmStep")
        .register_type_with_name::<ScriptFilmDocument>("FilmDocument");
    let doc = Rc::clone(&context);
    engine.register_fn("document", move || ScriptFilmDocument(Rc::clone(&doc)));
    engine.register_get("kind", |d: &mut ScriptFilmDocument| d.0.kind.clone());
    engine.register_get("name", |d: &mut ScriptFilmDocument| d.0.name.clone());
    engine.register_get("bodies", |d: &mut ScriptFilmDocument| d.0.bodies as i64);
    engine.register_fn("instances", document_paths);
    engine.register_fn(
        "value",
        |d: &mut ScriptFilmDocument, name: &str| -> Result<f64, Box<EvalAltResult>> {
            d.0.values.get(name).copied().ok_or_else(|| {
                let names: Vec<&str> = d.0.values.keys().map(String::as_str).take(24).collect();
                argument_error(
                    "name",
                    format!(
                        "the document has no parameter '{name}'; it has {}",
                        names.join(", ")
                    ),
                )
            })
        },
    );
    engine.register_fn("joints", document_joints);
    let (c1, c2) = (Rc::clone(&context), Rc::clone(&context));
    let (n1, n2) = (name.to_string(), name.to_string());
    engine.register_fn("storyboard", move || {
        storyboard(Rc::clone(&c1), &n1, &Dynamic::UNIT)
    });
    engine.register_fn("storyboard", move |options: Map| {
        storyboard(Rc::clone(&c2), &n2, &Dynamic::from_map(options))
    });
    engine.register_fn(
        "step",
        |ctx: NativeCallContext,
         f: &mut ScriptFilm,
         title: &str,
         subtitle: &str,
         seconds: Dynamic,
         body: FnPtr| { add_step(&ctx, f, "step", title, subtitle, &seconds, Some(body)) },
    );
    engine.register_fn(
        "step",
        |ctx: NativeCallContext, f: &mut ScriptFilm, title: &str, seconds: Dynamic, body: FnPtr| {
            add_step(&ctx, f, "step", title, "", &seconds, Some(body))
        },
    );
    engine.register_fn(
        "title",
        |ctx: NativeCallContext,
         f: &mut ScriptFilm,
         title: &str,
         subtitle: &str,
         seconds: Dynamic,
         body: FnPtr| { add_step(&ctx, f, "title", title, subtitle, &seconds, Some(body)) },
    );
    engine.register_fn(
        "title",
        |ctx: NativeCallContext,
         f: &mut ScriptFilm,
         title: &str,
         subtitle: &str,
         seconds: Dynamic| { add_step(&ctx, f, "title", title, subtitle, &seconds, None) },
    );
    engine.register_fn(
        "hold",
        |ctx: NativeCallContext, f: &mut ScriptFilm, seconds: Dynamic| {
            add_step(&ctx, f, "hold", "", "", &seconds, None)
        },
    );
    engine.register_fn(
        "hold",
        |ctx: NativeCallContext, f: &mut ScriptFilm, seconds: Dynamic, body: FnPtr| {
            add_step(&ctx, f, "hold", "", "", &seconds, Some(body))
        },
    );
    engine.register_fn(
        "orbit",
        |ctx: NativeCallContext, s: &mut ScriptStep, d: Dynamic| orbit(ctx, s, d, None),
    );
    engine.register_fn(
        "orbit",
        |ctx: NativeCallContext, s: &mut ScriptStep, d: Dynamic, o: Dynamic| {
            orbit(ctx, s, d, Some(o))
        },
    );
    engine.register_fn("camera", camera);
    engine.register_fn("move", movement);
    engine.register_fn("insert", super::fastener::insert);
    engine.register_fn("screw_in", super::fastener::screw_in);
    engine.register_fn(
        "turn",
        |ctx: NativeCallContext,
         s: &mut ScriptStep,
         p: Dynamic,
         a: ScriptAxis,
         f: Dynamic,
         t: Dynamic| { turn(ctx, s, p, a, f, t, None) },
    );
    engine.register_fn(
        "turn",
        |ctx: NativeCallContext,
         s: &mut ScriptStep,
         p: Dynamic,
         a: ScriptAxis,
         f: Dynamic,
         t: Dynamic,
         o: Dynamic| turn(ctx, s, p, a, f, t, Some(o)),
    );
    engine.register_fn(
        "joint",
        |ctx: NativeCallContext, s: &mut ScriptStep, p: &str, f: Dynamic, t: Dynamic| {
            joint(ctx, s, p, f, t, None)
        },
    );
    engine.register_fn(
        "joint",
        |ctx: NativeCallContext,
         s: &mut ScriptStep,
         p: &str,
         f: Dynamic,
         t: Dynamic,
         o: Dynamic| { joint(ctx, s, p, f, t, Some(o)) },
    );
    engine.register_fn(
        "explode",
        |ctx: NativeCallContext, s: &mut ScriptStep, f: Dynamic, t: Dynamic| {
            explode(ctx, s, f, t, None)
        },
    );
    engine.register_fn(
        "explode",
        |ctx: NativeCallContext, s: &mut ScriptStep, f: Dynamic, t: Dynamic, o: Dynamic| {
            explode(ctx, s, f, t, Some(o))
        },
    );
    engine.register_fn("section", section);
    engine.register_fn(
        "fade",
        |ctx: NativeCallContext, s: &mut ScriptStep, p: Dynamic, f: Dynamic, t: Dynamic| {
            fade(ctx, s, p, f, t, None)
        },
    );
    engine.register_fn(
        "fade",
        |ctx: NativeCallContext,
         s: &mut ScriptStep,
         p: Dynamic,
         f: Dynamic,
         t: Dynamic,
         o: Dynamic| { fade(ctx, s, p, f, t, Some(o)) },
    );
    engine.register_fn(
        "show",
        |ctx: NativeCallContext, s: &mut ScriptStep, p: Dynamic, v: bool| show(ctx, s, p, v, None),
    );
    engine.register_fn(
        "show",
        |ctx: NativeCallContext, s: &mut ScriptStep, p: Dynamic, v: bool, o: Dynamic| {
            show(ctx, s, p, v, Some(o))
        },
    );
    engine.register_fn(
        "highlight",
        |ctx: NativeCallContext, s: &mut ScriptStep, p: Dynamic| highlight(ctx, s, p, None),
    );
    engine.register_fn(
        "evolve",
        |ctx: NativeCallContext, s: &mut ScriptStep, d: Dynamic| {
            super::design::evolve(ctx, s, d, None)
        },
    );
    engine.register_fn(
        "evolve",
        |ctx: NativeCallContext, s: &mut ScriptStep, d: Dynamic, o: Dynamic| {
            super::design::evolve(ctx, s, d, Some(o))
        },
    );
    engine.register_fn(
        "highlight",
        |ctx: NativeCallContext, s: &mut ScriptStep, p: Dynamic, o: Dynamic| {
            highlight(ctx, s, p, Some(o))
        },
    );
}
