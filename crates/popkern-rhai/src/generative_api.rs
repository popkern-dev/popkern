//! Generative design (PK-956): `body.generative(#{ ... })` replaces the
//! body (the design space) with the stiffest design at a volume
//! fraction under the map's supports and loads.
//!
//! ```text
//! let bracket = space.generative(#{
//!     supports: [#{ faces: faces().labeled("wall"), kind: "fixed" }],
//!     loads: [#{ faces: faces().labeled("eye"), force: 200.0, direction: [0.0, 0.0, -1.0] }],
//!     keep: [faces().labeled("wall"), faces().labeled("eye")],
//!     volume_fraction: 0.25,
//!     filter: 3.0,
//!     build_direction: "+y",
//! });
//! ```

use popkern_history::feature_generative::{BuildDirectionDesc, GenerativeSupportKindDesc};
use popkern_script::selection::FaceSelector;
use popkern_script::{GenerativeLoad, GenerativeSpec};
use rhai::{Array, Dynamic, Engine, EvalAltResult, Map};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::{IntoScriptResult, argument_error};

/// Keys a generative map may carry.
const GENERATIVE_OPTIONS: &[&str] = &[
    "material",
    "supports",
    "loads",
    "keep",
    "volume_fraction",
    "filter",
    "keep_depth",
    "cell",
    "build_direction",
    "max_iterations",
    "robust",
    "max_cells",
    "label",
];
const SUPPORT_KEYS: &[&str] = &["faces", "kind"];
const LOAD_KEYS: &[&str] = &["faces", "force", "direction"];

/// The material a design defaults to (popkern live's printed parts).
const DEFAULT_MATERIAL: &str = "PLA (FDM)";
/// Updates a design runs at most by default (the bored-block fixture
/// settles in 212 at 0.5 mm cells).
const DEFAULT_MAX_ITERATIONS: u32 = 300;

fn selector(argument: &str, value: &Dynamic) -> Result<FaceSelector, Box<EvalAltResult>> {
    value.clone().try_cast::<FaceSelector>().ok_or_else(|| {
        argument_error(
            argument,
            format!(
                "expected a face selection (faces()...), got {}",
                value.type_name()
            ),
        )
    })
}

fn list(argument: &str, value: &Dynamic) -> Result<Array, Box<EvalAltResult>> {
    value.clone().into_array().map_err(|_| {
        argument_error(
            argument,
            format!("expected an array, got {}", value.type_name()),
        )
    })
}

fn entry(argument: &str, value: &Dynamic, keys: &[&str]) -> Result<Map, Box<EvalAltResult>> {
    let map = value.clone().try_cast::<Map>().ok_or_else(|| {
        argument_error(
            argument,
            format!("expected a map #{{ {} }}", keys.join(", ")),
        )
    })?;
    if let Some(k) = map.keys().find(|k| !keys.contains(&k.as_str())) {
        return Err(argument_error(
            argument,
            format!("unknown key '{k}'; expected {}", keys.join(", ")),
        ));
    }
    Ok(map)
}

fn required<'a>(
    map: &'a Map,
    argument: &str,
    key: &str,
) -> Result<&'a Dynamic, Box<EvalAltResult>> {
    map.get(key)
        .ok_or_else(|| argument_error(argument, format!("missing '{key}'")))
}

fn support(
    value: &Dynamic,
) -> Result<(FaceSelector, GenerativeSupportKindDesc), Box<EvalAltResult>> {
    let map = entry("supports", value, SUPPORT_KEYS)?;
    let faces = selector("supports.faces", required(&map, "supports", "faces")?)?;
    let kind = match map.get("kind") {
        Some(k) => args::spelling(
            "supports.kind",
            k,
            &[
                ("fixed", GenerativeSupportKindDesc::Fixed),
                ("surface_roller", GenerativeSupportKindDesc::SurfaceRoller),
            ],
        )?,
        None => GenerativeSupportKindDesc::Fixed,
    };
    Ok((faces, kind))
}

fn load(value: &Dynamic) -> Result<GenerativeLoad, Box<EvalAltResult>> {
    let map = entry("loads", value, LOAD_KEYS)?;
    let direction = args::vec3("loads.direction", required(&map, "loads", "direction")?)?;
    if direction.iter().all(|c| *c == 0.0) {
        return Err(argument_error("loads.direction", "the direction is zero"));
    }
    Ok(GenerativeLoad {
        faces: selector("loads.faces", required(&map, "loads", "faces")?)?,
        force: args::scalar("loads.force", required(&map, "loads", "force")?)?,
        direction,
    })
}

/// `"+x" | "-x" | "+y" | "-y" | "+z" | "-z"`.
fn build_direction(
    argument: &str,
    value: &Dynamic,
) -> Result<BuildDirectionDesc, Box<EvalAltResult>> {
    let d = |axis: u8, negative: bool| BuildDirectionDesc { axis, negative };
    args::spelling(
        argument,
        value,
        &[
            ("+x", d(0, false)),
            ("-x", d(0, true)),
            ("+y", d(1, false)),
            ("-y", d(1, true)),
            ("+z", d(2, false)),
            ("-z", d(2, true)),
        ],
    )
}

fn count(argument: &str, value: &Dynamic) -> Result<u32, Box<EvalAltResult>> {
    let v = args::float(argument, value)?;
    if !(v.is_finite() && v >= 1.0 && v.fract() == 0.0 && v <= f64::from(u32::MAX)) {
        return Err(argument_error(
            argument,
            format!("expected a whole number of at least 1, got {v}"),
        ));
    }
    Ok(v as u32)
}

fn generative(body: &mut ScriptBody, options: &Dynamic) -> Result<ScriptBody, Box<EvalAltResult>> {
    let o = Options::parse("generative", options, GENERATIVE_OPTIONS)?;
    let need = |key: &str| {
        o.get(key)
            .ok_or_else(|| argument_error(key, format!("generative needs '{key}'")))
    };
    let supports = list("supports", need("supports")?)?
        .iter()
        .map(support)
        .collect::<Result<Vec<_>, _>>()?;
    let loads = list("loads", need("loads")?)?
        .iter()
        .map(load)
        .collect::<Result<Vec<_>, _>>()?;
    let keep = match o.get("keep") {
        Some(v) => list("keep", v)?
            .iter()
            .map(|s| selector("keep", s))
            .collect::<Result<Vec<_>, _>>()?,
        None => Vec::new(),
    };
    let spec = GenerativeSpec {
        material: o
            .with("material", args::text)?
            .unwrap_or_else(|| DEFAULT_MATERIAL.to_string()),
        supports,
        loads,
        keep,
        volume_fraction: args::scalar("volume_fraction", need("volume_fraction")?)?,
        filter: args::scalar("filter", need("filter")?)?,
        keep_depth: o.with("keep_depth", args::scalar)?,
        cell: o.with("cell", args::scalar)?,
        build_direction: o.with("build_direction", build_direction)?,
        max_iterations: o
            .with("max_iterations", count)?
            .unwrap_or(DEFAULT_MAX_ITERATIONS),
        robust: o.with("robust", args::boolean)?.unwrap_or(false),
        max_cells: o.with("max_cells", count)?,
    };
    let served = body
        .body
        .generative(&mut body.model.borrow_mut(), &spec, &o.label()?)
        .script()?;
    Ok(ScriptBody::new(body.model.clone(), served))
}

/// Register `generative` on bodies.
pub fn register(engine: &mut Engine) {
    engine.register_fn("generative", |b: &mut ScriptBody, options: Dynamic| {
        generative(b, &options)
    });
}
