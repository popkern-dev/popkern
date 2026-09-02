//! Timeline / migration verbs and configurations (spec §§4.3, 7) on
//! the `Model` value: `insert_before`/`insert_after` (a closure that
//! records exactly one feature at the anchor), `reorder`,
//! `rollback_to`/`rollback_to_end`, `undo`/`redo`, and the
//! configuration surface. Every verb is transactional in the
//! recorder; the skin only spells it.

use popkern_script::{Body, Model, ScriptError};
use rhai::{Dynamic, Engine, EvalAltResult, FnPtr, NativeCallContext};

use crate::args;
use crate::body_api::ScriptBody;
use crate::error::{IntoScriptResult, argument_error};
use crate::feature_api::ScriptFeature;
use crate::model_api::ScriptModel;
use crate::values::report_map;

/// Run a script closure `|m, body| ...` against the recorder's `&mut
/// Model`: the session is moved into a temporary handle for the
/// duration of the call and moved back afterwards, so the closure's
/// method calls do not re-borrow the shared cell.
fn call_with_model(
    ctx: &NativeCallContext,
    build: &FnPtr,
    model: &mut Model,
    body: Body,
) -> Result<Body, ScriptError> {
    let taken = std::mem::replace(model, Model::new());
    let handle = ScriptModel::new(taken);
    let outcome: Result<ScriptBody, Box<EvalAltResult>> =
        build.call_within_context(ctx, (handle.clone(), ScriptBody::new(handle.clone(), body)));
    *model = handle.into_inner();
    match outcome {
        Ok(served) => Ok(served.body),
        Err(error) => Err(ScriptError::Argument {
            argument: "insert".into(),
            reason: error.to_string(),
        }),
    }
}

fn insert(
    ctx: NativeCallContext,
    m: &mut ScriptModel,
    anchor: ScriptFeature,
    build: FnPtr,
    before: bool,
) -> Result<ScriptFeature, Box<EvalAltResult>> {
    let handle = m.clone();
    let mut model = handle.borrow_mut();
    let closure = |model: &mut Model, body: Body| call_with_model(&ctx, &build, model, body);
    let recorded = if before {
        model.insert_before(&anchor.handle, closure)
    } else {
        model.insert_after(&anchor.handle, closure)
    }
    .script()?;
    drop(model);
    Ok(ScriptFeature::new(m.clone(), recorded))
}

fn config_name(value: &Dynamic) -> Result<Option<String>, Box<EvalAltResult>> {
    if value.is_unit() {
        Ok(None)
    } else {
        args::text("config", value).map(Some)
    }
}

/// Register the migration verbs and configuration surface.
pub fn register(engine: &mut Engine) {
    engine
        .register_fn(
            "insert_before",
            |ctx: NativeCallContext, m: &mut ScriptModel, anchor: ScriptFeature, build: FnPtr| {
                insert(ctx, m, anchor, build, true)
            },
        )
        .register_fn(
            "insert_after",
            |ctx: NativeCallContext, m: &mut ScriptModel, anchor: ScriptFeature, build: FnPtr| {
                insert(ctx, m, anchor, build, false)
            },
        )
        .register_fn(
            "reorder",
            |m: &mut ScriptModel, feature: ScriptFeature, index: i64| {
                if index < 0 {
                    return Err(argument_error("index", "timeline index must be >= 0"));
                }
                m.borrow_mut()
                    .reorder(&feature.handle, index as usize)
                    .script()
            },
        )
        .register_fn(
            "rollback_to",
            |m: &mut ScriptModel, feature: ScriptFeature| {
                m.borrow_mut()
                    .rollback_to(&feature.handle)
                    .script()
                    .map(|r| report_map(&r))
            },
        )
        .register_fn("rollback_to_end", |m: &mut ScriptModel| {
            m.borrow_mut()
                .rollback_to_end()
                .script()
                .map(|r| report_map(&r))
        })
        .register_fn("undo", |m: &mut ScriptModel| m.borrow_mut().undo().script())
        .register_fn("redo", |m: &mut ScriptModel| m.borrow_mut().redo().script())
        .register_fn(
            "add_config",
            |m: &mut ScriptModel, name: &str, description: &str| {
                m.borrow_mut().add_config(name, description).script()
            },
        )
        .register_fn("set_active_config", |m: &mut ScriptModel, name: Dynamic| {
            let name = config_name(&name)?;
            m.borrow_mut().set_active_config(name.as_deref()).script()
        })
        .register_fn("active_config", |m: &mut ScriptModel| {
            m.borrow_mut()
                .active_config()
                .map(Dynamic::from)
                .unwrap_or(Dynamic::UNIT)
        })
        .register_fn(
            "set_override",
            |m: &mut ScriptModel, config: &str, param: &str, value: Dynamic| {
                let value = args::float("value", &value)?;
                m.borrow_mut().set_override(config, param, value).script()
            },
        )
        .register_fn(
            "clear_override",
            |m: &mut ScriptModel, config: &str, param: &str| {
                m.borrow_mut().clear_override(config, param).script()
            },
        )
        .register_fn(
            "set_suppression",
            |m: &mut ScriptModel, config: &str, feature: ScriptFeature, suppressed: bool| {
                m.borrow_mut()
                    .set_suppression(config, &feature.handle, suppressed)
                    .script()
            },
        )
        .register_fn("configs", |m: &mut ScriptModel| {
            m.borrow_mut()
                .configs()
                .into_iter()
                .map(Dynamic::from)
                .collect::<rhai::Array>()
        })
        .register_fn("to_rhai", |m: &mut ScriptModel| {
            let mut model = m.borrow_mut();
            crate::emit::to_rhai(&mut model).map_err(|e| {
                crate::error::raise(serde_json::json!({
                    "kind": "NotApplicable",
                    "message": e.to_string(),
                    "detail": e.to_string(),
                }))
            })
        });
}
