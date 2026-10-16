//! Argument parsing from script values: scalars, units, enum
//! spellings, and trailing option maps. Every rejection is a
//! spec-§10 `ArgumentError` raised before any kernel work runs.

use popkern_history::feature::{JoinStyleKind, OperationMode, ShellDirectionKind, TrimPolicyKind};
use popkern_history::feature_blend::{
    CapKindDesc, ContinuityKind, CornerConstructionDesc, CornerModeDesc, InterferencePolicyKind,
};
use popkern_history::unit::Unit;
use popkern_part::{Dir, SurfaceKind};
use popkern_script::Scalar;
use rhai::{Dynamic, EvalAltResult, Map};

use crate::error::argument_error;

/// A script value as a numeric argument: a float or integer literal,
/// a [`Scalar`] handle (`p("wall")`, `expr("2*wall")`), or a string
/// treated as an expression (the Rust skin's rule).
pub fn scalar(argument: &str, value: &Dynamic) -> Result<Scalar, Box<EvalAltResult>> {
    if let Some(s) = value.clone().try_cast::<Scalar>() {
        return Ok(s);
    }
    if let Ok(f) = value.as_float() {
        return Ok(Scalar::Literal(f));
    }
    if let Ok(i) = value.as_int() {
        return Ok(Scalar::Literal(i as f64));
    }
    if let Ok(text) = value.clone().into_string() {
        return Ok(Scalar::Expr(text));
    }
    Err(argument_error(
        argument,
        format!(
            "expected a number, p(name), expr(text) or an expression string, got {}",
            value.type_name()
        ),
    ))
}

/// A plain float argument (heights, chord lengths, station parameters).
pub fn float(argument: &str, value: &Dynamic) -> Result<f64, Box<EvalAltResult>> {
    if let Ok(f) = value.as_float() {
        return Ok(f);
    }
    if let Ok(i) = value.as_int() {
        return Ok(i as f64);
    }
    Err(argument_error(
        argument,
        format!("expected a number, got {}", value.type_name()),
    ))
}

/// A `[x, y]` or `[x, y, z]` literal.
pub fn point(
    argument: &str,
    value: &Dynamic,
    arity: usize,
) -> Result<Vec<f64>, Box<EvalAltResult>> {
    let items = value.clone().into_array().map_err(|_| {
        argument_error(
            argument,
            format!(
                "expected a [..] array of {arity} numbers, got {}",
                value.type_name()
            ),
        )
    })?;
    if items.len() != arity {
        return Err(argument_error(
            argument,
            format!("expected {arity} coordinates, got {}", items.len()),
        ));
    }
    items.iter().map(|v| float(argument, v)).collect()
}

/// A `[x, y, z]` literal as a fixed array.
pub fn vec3(argument: &str, value: &Dynamic) -> Result<[f64; 3], Box<EvalAltResult>> {
    let v = point(argument, value, 3)?;
    Ok([v[0], v[1], v[2]])
}

/// A string argument.
/// A boolean argument.
pub fn boolean(argument: &str, value: &Dynamic) -> Result<bool, Box<EvalAltResult>> {
    value.as_bool().map_err(|_| {
        argument_error(
            argument,
            format!("expected true or false, got {}", value.type_name()),
        )
    })
}

pub fn text(argument: &str, value: &Dynamic) -> Result<String, Box<EvalAltResult>> {
    value.clone().into_string().map_err(|_| {
        argument_error(
            argument,
            format!("expected a string, got {}", value.type_name()),
        )
    })
}

pub(crate) fn spelling<T: Copy>(
    argument: &str,
    value: &Dynamic,
    table: &[(&str, T)],
) -> Result<T, Box<EvalAltResult>> {
    let word = text(argument, value)?;
    let key = word.trim().to_ascii_lowercase();
    table
        .iter()
        .find(|(name, _)| *name == key)
        .map(|(_, v)| *v)
        .ok_or_else(|| {
            let options: Vec<&str> = table.iter().map(|(n, _)| *n).collect();
            argument_error(
                argument,
                format!(
                    "unknown spelling '{word}'; expected one of {}",
                    options.join(" | ")
                ),
            )
        })
}

/// `plane | cylinder | cone | sphere | torus | nurbs | gregory`.
pub fn surface_kind(argument: &str, value: &Dynamic) -> Result<SurfaceKind, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("plane", SurfaceKind::Plane),
            ("cylinder", SurfaceKind::Cylinder),
            ("cone", SurfaceKind::Cone),
            ("sphere", SurfaceKind::Sphere),
            ("torus", SurfaceKind::Torus),
            ("nurbs", SurfaceKind::Nurbs),
            ("gregory", SurfaceKind::Gregory),
        ],
    )
}

/// `x | y | z`.
pub fn dir(argument: &str, value: &Dynamic) -> Result<Dir, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[("x", Dir::X), ("y", Dir::Y), ("z", Dir::Z)],
    )
}

/// `g0 | g1 | g2`.
pub fn continuity(argument: &str, value: &Dynamic) -> Result<ContinuityKind, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("g0", ContinuityKind::G0),
            ("g1", ContinuityKind::G1),
            ("g2", ContinuityKind::G2),
        ],
    )
}

/// `open | flat | rounded | tangent_match`.
pub fn cap(argument: &str, value: &Dynamic) -> Result<CapKindDesc, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("open", CapKindDesc::Open),
            ("flat", CapKindDesc::Flat),
            ("rounded", CapKindDesc::Rounded),
            ("tangent_match", CapKindDesc::TangentMatch),
        ],
    )
}

/// `smooth | mitre`.
pub fn corner(argument: &str, value: &Dynamic) -> Result<CornerModeDesc, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("smooth", CornerModeDesc::Smooth),
            ("mitre", CornerModeDesc::Mitre),
        ],
    )
}

/// `auto | setback | fitted | classic | rolling_ball` (PK-603, PK-861).
pub fn corner_cap(
    argument: &str,
    value: &Dynamic,
) -> Result<CornerConstructionDesc, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("auto", CornerConstructionDesc::Auto),
            ("setback", CornerConstructionDesc::Setback),
            ("fitted", CornerConstructionDesc::Fitted),
            ("classic", CornerConstructionDesc::Classic),
            ("rolling_ball", CornerConstructionDesc::RollingBall),
        ],
    )
}

/// `decline | serve_non_interfering`.
pub fn interference(
    argument: &str,
    value: &Dynamic,
) -> Result<InterferencePolicyKind, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("decline", InterferencePolicyKind::Decline),
            (
                "serve_non_interfering",
                InterferencePolicyKind::ServeNonInterfering,
            ),
        ],
    )
}

/// `inward | outward | midplane`.
pub fn shell_direction(
    argument: &str,
    value: &Dynamic,
) -> Result<ShellDirectionKind, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("inward", ShellDirectionKind::Inward),
            ("outward", ShellDirectionKind::Outward),
            ("midplane", ShellDirectionKind::Midplane),
        ],
    )
}

/// `arc | intersection`.
pub fn join_style(argument: &str, value: &Dynamic) -> Result<JoinStyleKind, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("arc", JoinStyleKind::Arc),
            ("intersection", JoinStyleKind::Intersection),
        ],
    )
}

/// `trim | decline`.
pub fn trim_policy(argument: &str, value: &Dynamic) -> Result<TrimPolicyKind, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("trim", TrimPolicyKind::Trim),
            ("decline", TrimPolicyKind::Decline),
        ],
    )
}

/// `new_body | boss | cut`.
pub fn mode(argument: &str, value: &Dynamic) -> Result<OperationMode, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("new_body", OperationMode::NewBody),
            ("boss", OperationMode::Boss),
            ("cut", OperationMode::Cut),
        ],
    )
}

/// Length/angle units by symbol or name.
pub fn unit(argument: &str, value: &Dynamic) -> Result<Unit, Box<EvalAltResult>> {
    spelling(
        argument,
        value,
        &[
            ("m", Unit::Meters),
            ("meters", Unit::Meters),
            ("mm", Unit::Millimeters),
            ("millimeters", Unit::Millimeters),
            ("cm", Unit::Centimeters),
            ("centimeters", Unit::Centimeters),
            ("in", Unit::Inches),
            ("inches", Unit::Inches),
            ("rad", Unit::Radians),
            ("radians", Unit::Radians),
            ("deg", Unit::Degrees),
            ("degrees", Unit::Degrees),
            ("ratio", Unit::Ratio),
        ],
    )
}

/// The script-visible spelling of a unit (inverse of [`unit`]).
pub fn unit_name(unit: Unit) -> &'static str {
    match unit {
        Unit::Meters => "m",
        Unit::Millimeters => "mm",
        Unit::Centimeters => "cm",
        Unit::Inches => "in",
        Unit::Radians => "rad",
        Unit::Degrees => "deg",
        Unit::Ratio => "ratio",
        Unit::SquareMeters => "m2",
        Unit::CubicMeters => "m3",
        Unit::PerMeter => "1/m",
        Unit::MetersToFourth => "m4",
    }
}

/// A trailing options map (`#{ label: "rounds", continuity: "g2" }`).
/// Unknown keys are argument errors so misspellings never pass
/// silently.
pub struct Options {
    map: Map,
    op: &'static str,
    allowed: &'static [&'static str],
}

impl Options {
    /// Parse `value` as the option map of `op`; `()` means no options.
    pub fn parse(
        op: &'static str,
        value: &Dynamic,
        allowed: &'static [&'static str],
    ) -> Result<Options, Box<EvalAltResult>> {
        let map = if value.is_unit() {
            Map::new()
        } else {
            value.clone().try_cast::<Map>().ok_or_else(|| {
                argument_error(
                    "options",
                    format!(
                        "{op}: expected an options map #{{..}}, got {}",
                        value.type_name()
                    ),
                )
            })?
        };
        for key in map.keys() {
            if !allowed.contains(&key.as_str()) {
                return Err(argument_error(
                    key.as_str(),
                    format!(
                        "{op} has no option '{key}'; options are {}",
                        allowed.join(", ")
                    ),
                ));
            }
        }
        Ok(Options { map, op, allowed })
    }

    /// The raw value of an option, if given.
    pub fn get(&self, key: &str) -> Option<&Dynamic> {
        debug_assert!(self.allowed.contains(&key), "{}: {key}", self.op);
        self.map.get(key)
    }

    /// A string option, defaulting to `""`.
    pub fn label(&self) -> Result<String, Box<EvalAltResult>> {
        match self.get("label") {
            Some(v) => text("label", v),
            None => Ok(String::new()),
        }
    }

    /// Apply `parse` to an option when present.
    pub fn with<T>(
        &self,
        key: &str,
        parse: impl FnOnce(&str, &Dynamic) -> Result<T, Box<EvalAltResult>>,
    ) -> Result<Option<T>, Box<EvalAltResult>> {
        self.get(key).map(|v| parse(key, v)).transpose()
    }

    /// [`Options::with`] for parsers shared across verbs whose allowed
    /// lists differ (a thread spec reads `left_hand`/`starts`, which a
    /// clearance hole's option list rightly omits — no thread is
    /// modeled). A key this verb does not allow reads as absent:
    /// `parse` already rejected it if the user passed it, so the
    /// allowed list stays the single authority on the verb's surface.
    pub fn with_if_allowed<T>(
        &self,
        key: &str,
        parse: impl FnOnce(&str, &Dynamic) -> Result<T, Box<EvalAltResult>>,
    ) -> Result<Option<T>, Box<EvalAltResult>> {
        if !self.allowed.contains(&key) {
            return Ok(None);
        }
        self.with(key, parse)
    }
}
