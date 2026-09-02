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
        .register_fn("set_expr", |m: &mut ScriptModel, name: &str, text: &str| {
            m.borrow_mut().set_expr(name, text).script()
        })
        .register_fn("clear_expr", |m: &mut ScriptModel, name: &str| {
            m.borrow_mut().clear_expr(name).script()
        })
        .register_fn(
            "set_bounds",
            |m: &mut ScriptModel, name: &str, min: Dynamic, max: Dynamic| {
                let (min, max) = (args::float("min", &min)?, args::float("max", &max)?);
                m.borrow_mut().set_bounds(name, min, max).script()
            },
        )
        .register_fn("clear_bounds", |m: &mut ScriptModel, name: &str| {
            m.borrow_mut().clear_bounds(name).script()
        })
        .register_fn(
            "rename_param",
            |m: &mut ScriptModel, old: &str, new: &str| {
                m.borrow_mut().rename_param(old, new).script()
            },
        )
        .register_fn(
            "add_derived_config",
            |m: &mut ScriptModel, name: &str, description: &str, parent: &str| {
                m.borrow_mut()
                    .add_derived_config(name, description, parent)
                    .script()
            },
        )
        .register_fn(
            "duplicate_config",
            |m: &mut ScriptModel, source: &str, new: &str| {
                m.borrow_mut().duplicate_config(source, new).script()
            },
        )
        .register_fn("delete_config", |m: &mut ScriptModel, name: &str| {
            m.borrow_mut().delete_config(name).script().map(|names| {
                names
                    .into_iter()
                    .map(Dynamic::from)
                    .collect::<rhai::Array>()
            })
        })
        .register_fn(
            "rename_config",
            |m: &mut ScriptModel, old: &str, new: &str| {
                m.borrow_mut().rename_config(old, new).script()
            },
        )
        .register_fn(
            "clear_suppression",
            |m: &mut ScriptModel, config: &str, feature: ScriptFeature| {
                m.borrow_mut()
                    .clear_suppression(config, &feature.handle)
                    .script()
            },
        )
        .register_fn("config", |m: &mut ScriptModel, name: &str| {
            let info = m.borrow_mut().config(name).script()?;
            let mut map = rhai::Map::new();
            map.insert("name".into(), info.name.into());
            map.insert("description".into(), info.description.into());
            map.insert(
                "parent".into(),
                info.parent.map(Dynamic::from).unwrap_or(Dynamic::UNIT),
            );
            map.insert(
                "overrides".into(),
                Dynamic::from_array(info.overrides.into_iter().map(Dynamic::from).collect()),
            );
            map.insert(
                "suppressions".into(),
                Dynamic::from_array(
                    info.suppressions
                        .into_iter()
                        .map(|(label, on)| {
                            let mut row = rhai::Map::new();
                            row.insert("feature".into(), label.into());
                            row.insert("suppressed".into(), on.into());
                            Dynamic::from_map(row)
                        })
                        .collect(),
                ),
            );
            Ok::<_, Box<EvalAltResult>>(Dynamic::from_map(map))
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
