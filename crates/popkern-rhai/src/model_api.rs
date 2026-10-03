//! The `Model` value: session, document, parameters, labels, bodies,
//! and report access (spec §§2, 4.1, 6.3, 8, 9).

use std::cell::{RefCell, RefMut};
use std::rc::Rc;

use popkern_core::tolerance::Tolerance;
use popkern_history::deferred_certification::ContinuityCertification;
use popkern_history::unit::Unit;
use popkern_script::{EdgeSelector, FaceSelector, Model, Plane};
use rhai::{Dynamic, Engine, EvalAltResult, FnPtr, NativeCallContext};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::{IntoScriptResult, argument_error};
use crate::feature_api::ScriptFeature;
use crate::profile_api::ScriptProfile;
use crate::values::{label_map, param_map, report_map};

/// A shared handle to the recorder session.
#[derive(Clone)]
pub struct ScriptModel(Rc<RefCell<Model>>);

impl std::fmt::Debug for ScriptModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Model")
    }
}

impl ScriptModel {
    /// Wrap a recorder session.
    pub fn new(model: Model) -> ScriptModel {
        ScriptModel(Rc::new(RefCell::new(model)))
    }

    /// Exclusive access to the session.
    pub fn borrow_mut(&self) -> RefMut<'_, Model> {
        self.0.borrow_mut()
    }

    /// Take the session out of the handle (the handle keeps an empty
    /// document).
    pub fn into_inner(self) -> Model {
        std::mem::replace(&mut *self.0.borrow_mut(), Model::new())
    }

    fn feature(&self, handle: popkern_script::FeatureHandle) -> ScriptFeature {
        ScriptFeature::new(self.clone(), handle)
    }
}

const MODEL_OPTIONS: &[&str] = &["unit", "tolerance"];

fn new_model(
    options: &Dynamic,
    certification: ContinuityCertification,
) -> Result<ScriptModel, Box<EvalAltResult>> {
    let options = Options::parse("model", options, MODEL_OPTIONS)?;
    let mut model = match options.with("tolerance", args::float)? {
        Some(t) if t > 0.0 => Model::with_tolerance(Tolerance::default().with_linear(t)),
        Some(t) => {
            return Err(argument_error(
                "tolerance",
                format!("must be positive, got {t}"),
            ));
        }
        None => Model::new(),
    };
    if let Some(unit) = options.with("unit", args::unit)? {
        model = model.unit(unit);
    }
    model.set_continuity_certification(certification);
    Ok(ScriptModel::new(model))
}

fn open(
    path: &str,
    certification: ContinuityCertification,
) -> Result<ScriptModel, Box<EvalAltResult>> {
    let mut model = Model::open(path).script()?;
    model.set_continuity_certification(certification);
    Ok(ScriptModel::new(model))
}

/// (Re-)register `model()`, `model(options)` and `open(path)` so every
/// model a script creates or opens certifies its blend seams under
/// `certification`. Registering again replaces the previous
/// constructors (same names and arities).
pub fn register_constructors(engine: &mut Engine, certification: ContinuityCertification) {
    engine
        .register_fn("model", move || new_model(&Dynamic::UNIT, certification))
        .register_fn("model", move |options: Dynamic| {
            new_model(&options, certification)
        })
        .register_fn("open", move |path: &str| open(path, certification));
}

fn tool(
    ctx: NativeCallContext,
    m: &mut ScriptModel,
    build: FnPtr,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let previous = m.borrow_mut().active().ok();
    let model = m.clone();
    // An error inside the closure is the failing call's own (class,
    // payload, line), not a generic closure-call failure (PK-634).
    let body: ScriptBody = build
        .call_within_context(&ctx, (model,))
        .map_err(crate::error::innermost)?;
    if let Some(previous) = previous {
        m.borrow_mut().set_active(previous).script()?;
    }
    Ok(body)
}

fn param(
    m: &mut ScriptModel,
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
    m.borrow_mut().param(name, value, unit).script()
}

fn param_expr(
    m: &mut ScriptModel,
    name: &str,
    text: &str,
    unit: Dynamic,
) -> Result<(), Box<EvalAltResult>> {
    let unit = if unit.is_unit() {
        Unit::Meters
    } else {
        args::unit("unit", &unit)?
    };
    m.borrow_mut().param_expr(name, text, unit).script()
}

fn loft(
    m: &mut ScriptModel,
    sections: rhai::Array,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("loft", options, &["label", "mode", "closed"])?;
    let bodies = sections
        .iter()
        .map(|d| {
            d.clone()
                .try_cast::<ScriptBody>()
                .map(|sb| sb.body)
                .ok_or_else(|| {
                    argument_error("sections", "each loft section must be a Body (a sketch)")
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let loft_options = popkern_script::LoftOptions {
        closed: options.with("closed", args::boolean)?.unwrap_or(false),
        mode: options
            .with("mode", args::mode)?
            .unwrap_or(popkern_history::feature::OperationMode::NewBody),
    };
    let body = m
        .borrow_mut()
        .loft_with(&bodies, loft_options, &options.label()?)
        .script()?;
    Ok(ScriptBody::new(m.clone(), body))
}

/// Register the model value and its global constructors (models
/// certify inline; a host re-registers the constructors with its own
/// policy through [`register_constructors`]).
pub fn register(engine: &mut Engine) {
    register_constructors(engine, ContinuityCertification::Inline);
    engine
        .register_type_with_name::<ScriptModel>("Model")
        .register_fn("save", |m: &mut ScriptModel, path: &str| {
            m.borrow_mut().save(path).script()
        })
        .register_fn("replay", |m: &mut ScriptModel| {
            m.borrow_mut().replay().script().map(|r| report_map(&r))
        })
        .register_fn("report", |m: &mut ScriptModel| {
            report_map(m.borrow_mut().report())
        })
        .register_fn("to_step", |m: &mut ScriptModel, path: &str| {
            m.borrow_mut().to_step(path).script()
        })
        .register_fn("assert_watertight", |m: &mut ScriptModel| {
            m.borrow_mut().assert_watertight().script()
        })
        .register_fn(
            "assert_step_roundtrip",
            |m: &mut ScriptModel, path: &str| m.borrow_mut().assert_step_roundtrip(path).script(),
        )
        .register_get("tolerance", |m: &mut ScriptModel| {
            m.borrow_mut().tolerance().linear()
        })
        .register_fn("loft", |m: &mut ScriptModel, sections: rhai::Array| {
            loft(m, sections, &Dynamic::UNIT)
        })
        .register_fn(
            "loft",
            |m: &mut ScriptModel, sections: rhai::Array, o: Dynamic| loft(m, sections, &o),
        )
        .register_fn("profile", |m: &mut ScriptModel, plane: Plane| {
            ScriptProfile::new(m.clone(), plane)
        })
        .register_fn("tool", tool)
        .register_fn("active", |m: &mut ScriptModel| {
            let body = m.borrow_mut().active().script()?;
            Ok::<_, Box<EvalAltResult>>(ScriptBody::new(m.clone(), body))
        })
        .register_fn("set_active", |m: &mut ScriptModel, body: ScriptBody| {
            m.borrow_mut().set_active(body.body).script()
        })
        .register_fn("bodies", |m: &mut ScriptModel| {
            let bodies = m.borrow_mut().alive_bodies();
            bodies
                .into_iter()
                .map(|b| Dynamic::from(ScriptBody::new(m.clone(), b)))
                .collect::<rhai::Array>()
        })
        .register_fn("body", |m: &mut ScriptModel, label: &str| {
            let body = m.borrow_mut().body(label).script()?;
            Ok::<_, Box<EvalAltResult>>(ScriptBody::new(m.clone(), body))
        })
        .register_fn(
            "label_body",
            |m: &mut ScriptModel, name: &str, body: ScriptBody| {
                m.borrow_mut().label_body(name, body.body).script()
            },
        )
        .register_fn(
            "param",
            |m: &mut ScriptModel, name: &str, value: Dynamic| param(m, name, value, Dynamic::UNIT),
        )
        .register_fn("param", param)
        .register_fn(
            "param_expr",
            |m: &mut ScriptModel, name: &str, text: &str| param_expr(m, name, text, Dynamic::UNIT),
        )
        .register_fn("param_expr", param_expr)
        .register_fn(
            "set_param",
            |m: &mut ScriptModel, name: &str, value: Dynamic| {
                let value = args::float("value", &value)?;
                m.borrow_mut().set_param(name, value).script()
            },
        )
        .register_fn("params", |m: &mut ScriptModel| {
            m.borrow_mut()
                .params()
                .iter()
                .map(param_map)
                .collect::<rhai::Array>()
        })
        .register_fn("feature", |m: &mut ScriptModel, label: &str| {
            let handle = m.borrow_mut().feature(label).script()?;
            Ok::<_, Box<EvalAltResult>>(m.feature(handle))
        })
        .register_fn("features", |m: &mut ScriptModel| {
            let handles = m.borrow_mut().features();
            handles
                .into_iter()
                .map(|h| Dynamic::from(m.feature(h)))
                .collect::<rhai::Array>()
        })
        .register_fn(
            "rename",
            |m: &mut ScriptModel, feature: ScriptFeature, label: &str| {
                m.borrow_mut().rename(&feature.handle, label).script()
            },
        )
        .register_fn(
            "label_edge",
            |m: &mut ScriptModel, name: &str, sel: EdgeSelector| {
                m.borrow_mut().label_edge(name, &sel).script().map(|_| ())
            },
        )
        .register_fn(
            "label_face",
            |m: &mut ScriptModel, name: &str, sel: FaceSelector| {
                m.borrow_mut().label_face(name, &sel).script().map(|_| ())
            },
        )
        .register_fn("unlabel", |m: &mut ScriptModel, name: &str| {
            m.borrow_mut().unlabel(name).script().map(|_| ())
        })
        .register_fn(
            "rename_label",
            |m: &mut ScriptModel, old: &str, new: &str| {
                m.borrow_mut().rename_label(old, new).script()
            },
        )
        .register_fn("labels", |m: &mut ScriptModel| {
            m.borrow_mut()
                .labels()
                .iter()
                .map(|(l, n)| label_map(l, n))
                .collect::<rhai::Array>()
        })
        .register_fn(
            "label_entity",
            |m: &mut ScriptModel, name: &str, spelled: &str| {
                m.borrow_mut()
                    .label_spelled(name, spelled)
                    .script()
                    .map(|_| ())
            },
        )
        .register_fn("entity_spelling", |m: &mut ScriptModel, label: &str| {
            let model = m.borrow_mut();
            let entity = model.entity_by_label(label).script()?;
            model.spell_entity(&entity).script()
        })
        .register_fn("to_string", |_: &mut ScriptModel| "Model".to_string());
    crate::migration_api::register(engine);
}
