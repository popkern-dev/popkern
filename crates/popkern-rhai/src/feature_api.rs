//! The `Feature` value: a recorded feature's handle with the stable
//! interrogation and edit verbs of spec §§7–8 (`describe`,
//! `evidence`, `modify`, `suppress`/`unsuppress`, `delete`,
//! `dependents`).

use popkern_script::{EdgeSelector, Edit, FaceSelector, FeatureHandle};
use rhai::{Dynamic, Engine, EvalAltResult, Map};

use crate::args;
use crate::error::{IntoScriptResult, argument_error};
use crate::model_api::ScriptModel;
use crate::values::to_value;

/// A recorded feature, bound to its model.
#[derive(Debug, Clone)]
pub struct ScriptFeature {
    pub(crate) model: ScriptModel,
    pub(crate) handle: FeatureHandle,
}

/// Keys a `modify` map accepts: scalar roles plus the selection and
/// continuity edits.
const MODIFY_KEYS: &[&str] = &["continuity", "edges", "faces"];

impl ScriptFeature {
    pub(crate) fn new(model: ScriptModel, handle: FeatureHandle) -> ScriptFeature {
        ScriptFeature { model, handle }
    }

    fn edit_from(map: &Map) -> Result<Edit, Box<EvalAltResult>> {
        let mut edit = Edit::new();
        for (key, value) in map {
            let key = key.as_str();
            edit =
                match key {
                    "continuity" => edit.continuity(args::continuity(key, value)?),
                    "edges" => {
                        edit.edges(value.clone().try_cast::<EdgeSelector>().ok_or_else(|| {
                            argument_error("edges", "expected an edges() selection")
                        })?)
                    }
                    "faces" => {
                        edit.faces(value.clone().try_cast::<FaceSelector>().ok_or_else(|| {
                            argument_error("faces", "expected a faces() selection")
                        })?)
                    }
                    role => edit.scalar(role, args::scalar(role, value)?),
                };
        }
        let _ = MODIFY_KEYS;
        Ok(edit)
    }

    fn modify(&self, edits: Dynamic) -> Result<(), Box<EvalAltResult>> {
        let map = edits
            .try_cast::<Map>()
            .ok_or_else(|| argument_error("edits", "modify expects a map #{ role: value, .. }"))?;
        let edit = Self::edit_from(&map)?;
        self.model.borrow_mut().modify(&self.handle, edit).script()
    }

    fn param(&self, role: &str) -> Result<Dynamic, Box<EvalAltResult>> {
        let name = format!(
            "{}.{role}",
            popkern_history::labels::qualified_prefix(&self.handle.label)
        );
        let model = self.model.borrow_mut();
        model
            .params()
            .into_iter()
            .find(|p| p.name == name)
            .map(|p| Dynamic::from(p.value))
            .ok_or_else(|| {
                argument_error(
                    "role",
                    format!("feature '{}' has no parameter '{role}'", self.handle.label),
                )
            })
    }
}

/// Register the feature value.
pub fn register(engine: &mut Engine) {
    engine
        .register_type_with_name::<ScriptFeature>("Feature")
        .register_get("label", |f: &mut ScriptFeature| f.handle.label.clone())
        .register_get("operation", |f: &mut ScriptFeature| {
            f.handle.operation.clone()
        })
        .register_fn("describe", |f: &mut ScriptFeature| {
            f.model.borrow_mut().describe(&f.handle).script()
        })
        .register_fn("evidence", |f: &mut ScriptFeature| {
            f.model
                .borrow_mut()
                .evidence(&f.handle)
                .map(to_value)
                .unwrap_or(Dynamic::UNIT)
        })
        .register_fn("param", |f: &mut ScriptFeature, role: &str| f.param(role))
        .register_fn("modify", |f: &mut ScriptFeature, edits: Dynamic| {
            f.modify(edits)
        })
        .register_fn("suppress", |f: &mut ScriptFeature| {
            f.model.borrow_mut().suppress(&f.handle).script()
        })
        .register_fn("unsuppress", |f: &mut ScriptFeature| {
            f.model.borrow_mut().unsuppress(&f.handle).script()
        })
        .register_fn("delete", |f: &mut ScriptFeature| {
            f.model
                .borrow_mut()
                .delete(&f.handle)
                .script()
                .map(|deleted| deleted.len() as i64)
        })
        .register_fn("dependents", |f: &mut ScriptFeature| {
            let handles = f.model.borrow_mut().dependents(&f.handle);
            handles
                .into_iter()
                .map(|h| Dynamic::from(ScriptFeature::new(f.model.clone(), h)))
                .collect::<rhai::Array>()
        })
        .register_fn("to_string", |f: &mut ScriptFeature| {
            format!("{} ({})", f.handle.label, f.handle.operation)
        });
}
