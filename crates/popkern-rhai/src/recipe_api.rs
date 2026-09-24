//! Fillet recipes, chamfer kinds, and shell options as script values
//! (spec §5 `fillet`/`chamfer`/`shell` optional arguments): parsed
//! from a trailing options map, with the parametrised forms
//! (`conic(rho)`, `linear(end)`, `by_stations([...])`,
//! `two_distance(face, d2)`, `distance_angle(face, angle)`) as
//! constructor functions.

use popkern_history::feature_blend::{BlendProfileKind, FilletRecipeDesc};
use popkern_script::{ChamferSpec, FaceSelector, RadiusLaw, Scalar, ShellOptions};
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::error::argument_error;

/// A blend profile value (`conic(rho)`, or the strings `circular` /
/// `curvature_continuous`).
#[derive(Debug, Clone, Copy)]
pub struct ScriptProfileKind(pub BlendProfileKind);

/// A radius-law value (`constant()`, `linear(end)`, `by_stations([[t, r], …])`).
#[derive(Debug, Clone)]
pub struct ScriptRadiusLaw(pub RadiusLaw);

/// A chamfer-kind value (`equal()`, `two_distance(face, d2)`,
/// `distance_angle(face, angle)`).
#[derive(Debug, Clone)]
pub struct ScriptChamferKind(pub ChamferSpec);

/// Options a `fillet` call accepts.
pub const FILLET_OPTIONS: &[&str] = &[
    "label",
    "continuity",
    "cap",
    "profile",
    "interference",
    "law",
    "corner",
];

/// Options a `chamfer` call accepts.
pub const CHAMFER_OPTIONS: &[&str] = &["label", "kind"];

/// Options a `shell` call accepts.
pub const SHELL_OPTIONS: &[&str] = &["label", "direction", "join", "trim"];

fn profile_kind(argument: &str, value: &Dynamic) -> Result<BlendProfileKind, Box<EvalAltResult>> {
    if let Some(kind) = value.clone().try_cast::<ScriptProfileKind>() {
        return Ok(kind.0);
    }
    let word = args::text(argument, value)?;
    match word.as_str() {
        "circular" => Ok(BlendProfileKind::Circular),
        "curvature_continuous" => Ok(BlendProfileKind::CurvatureContinuous),
        "conic" => Err(argument_error(
            argument,
            "conic needs a shoulder fullness: profile: conic(rho)",
        )),
        other => Err(argument_error(
            argument,
            format!(
                "unknown profile '{other}'; expected circular | conic(rho) | curvature_continuous"
            ),
        )),
    }
}

fn radius_law(argument: &str, value: &Dynamic) -> Result<RadiusLaw, Box<EvalAltResult>> {
    if let Some(law) = value.clone().try_cast::<ScriptRadiusLaw>() {
        return Ok(law.0);
    }
    let word = args::text(argument, value)?;
    match word.as_str() {
        "constant" => Ok(RadiusLaw::Constant),
        other => Err(argument_error(
            argument,
            format!(
                "unknown law '{other}'; expected constant | linear(end) | by_stations([[t, r], ..])"
            ),
        )),
    }
}

fn chamfer_kind(argument: &str, value: &Dynamic) -> Result<ChamferSpec, Box<EvalAltResult>> {
    if let Some(kind) = value.clone().try_cast::<ScriptChamferKind>() {
        return Ok(kind.0);
    }
    let word = args::text(argument, value)?;
    match word.as_str() {
        "equal" => Ok(ChamferSpec::Equal),
        other => Err(argument_error(
            argument,
            format!(
                "unknown chamfer kind '{other}'; expected equal | two_distance(face, d2) | distance_angle(face, angle)"
            ),
        )),
    }
}

/// The fillet recipe, radius law, and label from an options map.
pub fn fillet_options(
    value: &Dynamic,
) -> Result<(FilletRecipeDesc, RadiusLaw, String), Box<EvalAltResult>> {
    let options = Options::parse("fillet", value, FILLET_OPTIONS)?;
    let mut recipe = FilletRecipeDesc::default();
    if let Some(c) = options.with("continuity", args::continuity)? {
        recipe.continuity = c;
    }
    if let Some(c) = options.with("cap", args::cap)? {
        recipe.cap = c;
    }
    if let Some(p) = options.with("profile", profile_kind)? {
        recipe.profile = p;
    }
    if let Some(i) = options.with("interference", args::interference)? {
        recipe.interference = i;
    }
    if let Some(c) = options.with("corner", args::corner)? {
        recipe.corner = c;
    }
    let law = options
        .with("law", radius_law)?
        .unwrap_or(RadiusLaw::Constant);
    Ok((recipe, law, options.label()?))
}

/// The chamfer kind and label from an options map.
pub fn chamfer_options(value: &Dynamic) -> Result<(ChamferSpec, String), Box<EvalAltResult>> {
    let options = Options::parse("chamfer", value, CHAMFER_OPTIONS)?;
    let kind = options
        .with("kind", chamfer_kind)?
        .unwrap_or(ChamferSpec::Equal);
    Ok((kind, options.label()?))
}

/// The shell options and label from an options map.
pub fn shell_options(value: &Dynamic) -> Result<(ShellOptions, String), Box<EvalAltResult>> {
    let options = Options::parse("shell", value, SHELL_OPTIONS)?;
    let mut shell = ShellOptions::default();
    if let Some(d) = options.with("direction", args::shell_direction)? {
        shell.direction = d;
    }
    if let Some(j) = options.with("join", args::join_style)? {
        shell.join = j;
    }
    if let Some(t) = options.with("trim", args::trim_policy)? {
        shell.trim = t;
    }
    Ok((shell, options.label()?))
}

fn stations(value: Dynamic) -> Result<Vec<(f64, Scalar)>, Box<EvalAltResult>> {
    let items = value
        .into_array()
        .map_err(|_| argument_error("law", "by_stations expects [[t, radius], ..]"))?;
    items
        .iter()
        .map(|item| {
            let pair = item
                .clone()
                .into_array()
                .map_err(|_| argument_error("law", "each station is a [t, radius] pair"))?;
            if pair.len() != 2 {
                return Err(argument_error("law", "each station is a [t, radius] pair"));
            }
            Ok((
                args::float("t", &pair[0])?,
                args::scalar("radius", &pair[1])?,
            ))
        })
        .collect()
}

/// Register the recipe value constructors.
pub fn register(engine: &mut Engine) {
    engine
        .register_type_with_name::<ScriptProfileKind>("BlendProfile")
        .register_fn("circular", || ScriptProfileKind(BlendProfileKind::Circular))
        .register_fn("curvature_continuous", || {
            ScriptProfileKind(BlendProfileKind::CurvatureContinuous)
        })
        .register_fn("conic", |rho: Dynamic| {
            let rho = args::float("rho", &rho)?;
            if !(rho > 0.0 && rho < 1.0) {
                return Err(argument_error(
                    "rho",
                    format!("conic shoulder fullness must be in (0, 1), got {rho}"),
                ));
            }
            Ok::<_, Box<EvalAltResult>>(ScriptProfileKind(BlendProfileKind::Conic { rho }))
        });
    engine
        .register_type_with_name::<ScriptRadiusLaw>("RadiusLaw")
        .register_fn("constant", || ScriptRadiusLaw(RadiusLaw::Constant))
        .register_fn("linear", |end: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(ScriptRadiusLaw(RadiusLaw::Linear {
                end: args::scalar("end", &end)?,
            }))
        })
        .register_fn("by_stations", |list: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(ScriptRadiusLaw(RadiusLaw::ByStations(stations(list)?)))
        });
    engine
        .register_type_with_name::<ScriptChamferKind>("ChamferKind")
        .register_fn("equal", || ScriptChamferKind(ChamferSpec::Equal))
        .register_fn("two_distance", |face: FaceSelector, d2: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(ScriptChamferKind(ChamferSpec::TwoDistance {
                face,
                distance_2: args::scalar("distance_2", &d2)?,
            }))
        })
        .register_fn("distance_angle", |face: FaceSelector, angle: Dynamic| {
            Ok::<_, Box<EvalAltResult>>(ScriptChamferKind(ChamferSpec::DistanceAngle {
                face,
                angle: args::scalar("angle", &angle)?,
            }))
        });
}
