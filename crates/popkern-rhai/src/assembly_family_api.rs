//! Part families (design §6.1, PK-630): `a.part(name, defaults, |m, p| …)`
//! declares a parametric definition and builds its default variant;
//! `definition.variant(#{ … })` returns the variant for a parameter set,
//! running the closure only for a set not built before.
//!
//! ```rhai
//! let bolt = a.part("bolt", #{ size: "M6", length: 14.0 }, |m, p| {
//!     m.bolt(p.size, p.length, #{ label: "bolt" });
//! });
//! a.place(bolt.variant(#{ size: "M4", length: 10.0 }), "small");
//! ```

use std::collections::BTreeMap;

use popkern_assembly::{FamilyVariant, ParamValue};
use rhai::{Dynamic, EvalAltResult, FnPtr, Map, NativeCallContext};

use crate::assembly_api::{
    Family, IntoAssemblyResult, ScriptAssembly, ScriptDefinition, record_part,
};
use crate::error::argument_error;

fn param_value(key: &str, value: &Dynamic) -> Result<ParamValue, Box<EvalAltResult>> {
    if let Ok(b) = value.as_bool() {
        return Ok(ParamValue::Bool(b));
    }
    if let Ok(i) = value.as_int() {
        return Ok(ParamValue::Number(i as f64));
    }
    if let Ok(f) = value.as_float() {
        return Ok(ParamValue::Number(f));
    }
    if let Ok(s) = value.clone().into_string() {
        return Ok(ParamValue::Text(s));
    }
    Err(argument_error(
        key,
        format!(
            "a family parameter must be a number, string or bool, got {}",
            value.type_name()
        ),
    ))
}

fn to_map(params: &BTreeMap<String, ParamValue>) -> Map {
    params
        .iter()
        .map(|(k, v)| {
            let value = match v {
                ParamValue::Bool(b) => Dynamic::from_bool(*b),
                ParamValue::Number(n) => Dynamic::from_float(*n),
                ParamValue::Text(s) => Dynamic::from(s.clone()),
            };
            (k.as_str().into(), value)
        })
        .collect()
}

/// `a.part(name, defaults, |m, p| …)`.
pub(crate) fn family(
    ctx: NativeCallContext,
    a: &mut ScriptAssembly,
    name: &str,
    defaults: Map,
    build: FnPtr,
) -> Result<ScriptDefinition, Box<EvalAltResult>> {
    let mut params = BTreeMap::new();
    for (key, value) in &defaults {
        params.insert(key.to_string(), param_value(key, value)?);
    }
    let model = record_part(&ctx, a, name, &build, Some(to_map(&params)))?;
    let variant = FamilyVariant {
        family: name.to_string(),
        params: params.clone(),
    };
    let id = a
        .borrow_mut()
        .assembly
        .add_family_variant(name, model, variant)
        .assembly()?;
    let index = {
        let mut state = a.borrow_mut();
        state.families.push(Family {
            name: name.to_string(),
            defaults: params,
            builder: build,
        });
        state.families.len() - 1
    };
    Ok(ScriptDefinition {
        owner: a.clone(),
        id,
        family: Some(index),
    })
}

/// `definition.variant(#{ … })` (`with` is a Rhai keyword).
pub(crate) fn variant(
    ctx: NativeCallContext,
    d: &mut ScriptDefinition,
    overrides: Map,
) -> Result<ScriptDefinition, Box<EvalAltResult>> {
    let Some(index) = d.family else {
        return Err(argument_error(
            "variant",
            "this part was not declared with parameters: declare it as a.part(name, #{ … }, |m, p| …)",
        ));
    };
    let (family_name, mut params, builder) = {
        let state = d.owner.borrow_mut();
        let family = &state.families[index];
        (
            family.name.clone(),
            family.defaults.clone(),
            family.builder.clone(),
        )
    };
    let declared: Vec<String> = params.keys().cloned().collect();
    for (key, value) in &overrides {
        let value = param_value(key, value)?;
        let Some(default) = params.get(key.as_str()) else {
            return Err(argument_error(
                key,
                format!(
                    "'{family_name}' has no parameter '{key}' (parameters: {})",
                    declared.join(", ")
                ),
            ));
        };
        if std::mem::discriminant(default) != std::mem::discriminant(&value) {
            return Err(argument_error(
                key,
                format!(
                    "'{family_name}' parameter '{key}' is {}, got {}",
                    default.type_name(),
                    value.type_name()
                ),
            ));
        }
        params.insert(key.to_string(), value);
    }
    let existing = {
        let state = d.owner.borrow_mut();
        state.assembly.definitions().iter().position(|def| {
            def.family
                .as_ref()
                .is_some_and(|v| v.matches(&family_name, &params))
        })
    };
    if let Some(i) = existing {
        return Ok(ScriptDefinition {
            owner: d.owner.clone(),
            id: popkern_assembly::DefinitionId(i as u32),
            family: Some(index),
        });
    }
    let defaults = d.owner.borrow_mut().families[index].defaults.clone();
    let variant = FamilyVariant {
        family: family_name,
        params: params.clone(),
    };
    let name = variant.display_name(&defaults);
    let model = record_part(&ctx, &d.owner, &name, &builder, Some(to_map(&params)))?;
    let id = d
        .owner
        .borrow_mut()
        .assembly
        .add_family_variant(&name, model, variant)
        .assembly()?;
    Ok(ScriptDefinition {
        owner: d.owner.clone(),
        id,
        family: Some(index),
    })
}
