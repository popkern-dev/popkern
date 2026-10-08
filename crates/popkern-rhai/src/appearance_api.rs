//! Part appearance in the Rhai skin (six-axis-arm §4.5):
//!
//! ```rhai
//! m.appearance(#{ color: "#2b5cab", metal: 0.6, rough: 0.35 });
//! a.appearance(housing, #{ color: "aluminium" });
//! ```
//!
//! A colour is `"#rrggbb"`, a name from [`NAMED`], or `[r, g, b]` in
//! `[0, 1]`; `metal` (default 0) and `rough` (default 0.5) lie in
//! `[0, 1]`. Appearance is display and exchange data: STEP carries the
//! colour, popkern live the finish too.

use popkern_history::appearance::Appearance;
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::assembly_api::{IntoAssemblyResult, ScriptDefinition};
use crate::error::argument_error;
use crate::model_api::ScriptModel;

/// Named colours (sRGB, `[0, 1]`): the materials a mechanism shows.
pub const NAMED: &[(&str, [f64; 3])] = &[
    ("aluminium", [0.77, 0.78, 0.80]),
    ("steel", [0.55, 0.56, 0.58]),
    ("black", [0.06, 0.06, 0.07]),
    ("white", [0.93, 0.93, 0.92]),
    ("copper", [0.72, 0.45, 0.20]),
    ("brass", [0.78, 0.64, 0.30]),
    ("blue", [0.16, 0.36, 0.67]),
    ("red", [0.75, 0.15, 0.13]),
    ("orange", [0.95, 0.47, 0.10]),
    ("yellow", [0.95, 0.80, 0.15]),
    ("green", [0.20, 0.55, 0.30]),
    ("gray", [0.45, 0.46, 0.48]),
];

fn colour(value: &Dynamic) -> Result<[f64; 3], Box<EvalAltResult>> {
    let accepted = || {
        let names: Vec<&str> = NAMED.iter().map(|(n, _)| *n).collect();
        format!(
            "expected \"#rrggbb\", [r, g, b] in [0, 1], or one of: {}",
            names.join(", ")
        )
    };
    if let Ok(text) = value.clone().into_immutable_string() {
        let text = text.trim().to_ascii_lowercase();
        if let Some(hex) = text.strip_prefix('#') {
            let byte = |i: usize| {
                hex.get(i..i + 2)
                    .and_then(|h| u8::from_str_radix(h, 16).ok())
                    .map(|b| f64::from(b) / 255.0)
            };
            return match (hex.len(), byte(0), byte(2), byte(4)) {
                (6, Some(r), Some(g), Some(b)) => Ok([r, g, b]),
                _ => Err(argument_error("color", accepted())),
            };
        }
        return NAMED
            .iter()
            .find(|(n, _)| *n == text)
            .map(|(_, c)| *c)
            .ok_or_else(|| argument_error("color", accepted()));
    }
    let a = value
        .clone()
        .try_cast::<rhai::Array>()
        .filter(|a| a.len() == 3)
        .ok_or_else(|| argument_error("color", accepted()))?;
    Ok([
        args::float("color", &a[0])?,
        args::float("color", &a[1])?,
        args::float("color", &a[2])?,
    ])
}

/// Parse an appearance option map.
fn parse(map: &Dynamic) -> Result<Appearance, Box<EvalAltResult>> {
    let o = Options::parse("appearance", map, &["color", "metal", "rough"])?;
    let color = o
        .with("color", |_, v| colour(v))?
        .ok_or_else(|| argument_error("color", "an appearance needs a color"))?;
    let metal = o.with("metal", args::float)?.unwrap_or(0.0);
    let rough = o.with("rough", args::float)?.unwrap_or(0.5);
    Appearance::new(color, metal, rough).map_err(|e| argument_error(e.field, e.to_string()))
}

/// Register `appearance` on models and assembly definitions.
pub fn register(engine: &mut Engine) {
    engine
        .register_fn("appearance", |m: &mut ScriptModel, map: Dynamic| {
            let appearance = parse(&map)?;
            m.borrow_mut().set_appearance(Some(appearance));
            Ok::<_, Box<EvalAltResult>>(())
        })
        .register_fn(
            "appearance",
            |a: &mut crate::assembly_api::ScriptAssembly, d: ScriptDefinition, map: Dynamic| {
                let appearance = parse(&map)?;
                let mut state = a.borrow_mut();
                state
                    .assembly
                    .part_model_mut(d.id)
                    .assembly()?
                    .set_appearance(Some(appearance));
                Ok::<_, Box<EvalAltResult>>(())
            },
        );
}
