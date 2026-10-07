//! Fastener verbs (PK-460): `bolt`, `threaded_rod`, `nut` on the
//! model; `tapped_hole`, `clearance_hole` on a body.
//!
//! Thread specs are designation strings (`"M6"`, `"M6x0.75"`,
//! `"1/4-20"`, `"#8-32"` — millimetre-magnitude numbers, bound in
//! model units) or maps `#{ major: .., pitch: .. }`. Author fastener
//! models in a default-unit session with millimetre-magnitude
//! numbers (an M6 bolt is 6.0 units across). Fastener bodies build in the insertion frame (axis +Z,
//! origin at the bearing plane) and place with `translate`/`rotate`;
//! hole cuts take an `axis_at(origin, direction)` whose direction
//! points into the material, so a fastener on the same axis mates by
//! construction.

use popkern_history::feature_fastener::FastenerKind;
use popkern_history::unit::Unit;
use popkern_script::body_fastener::{
    ClearanceFitArg, ClearanceHoleArgs, FastenerOptions, TapHoleOptions, ThreadSpecArg,
};
use popkern_script::scalar::Scalar;
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::{IntoScriptResult, argument_error};
use crate::model_api::ScriptModel;
use crate::profile_api::ScriptEnd;

const BOLT_OPTIONS: &[&str] = &[
    "label",
    "head",
    "thread_length",
    "left_hand",
    "starts",
    "tip_chamfer",
];
const ROD_OPTIONS: &[&str] = &["label", "left_hand", "starts", "tip_chamfer"];
const NUT_OPTIONS: &[&str] = &[
    "label",
    "clearance",
    "thickness",
    "left_hand",
    "starts",
    "corner_chamfer",
    "mouth_chamfer",
];
const TAPPED_OPTIONS: &[&str] = &["label", "clearance", "left_hand", "starts", "entry_chamfer"];
const CLEARANCE_OPTIONS: &[&str] = &[
    "label",
    "fit",
    "diameter",
    "countersink",
    "counterbore",
    "entry_chamfer",
];

/// Parse a thread spec argument: a designation string (millimetre
/// dimensions, converted into the model's length `unit`) or a
/// `#{ major, pitch }` map (model units). A metre (default-unit)
/// model keeps the convention that fastener sizes are millimetre
/// MAGNITUDES (an "M6" records 6.0); a mm, cm or inch model records
/// the true length — an M6 in an inch model is 0.236 in, not 6 in
/// (PK-462).
fn thread_spec(
    value: &Dynamic,
    options: &Options,
    unit: Unit,
) -> Result<ThreadSpecArg, Box<EvalAltResult>> {
    let per_mm = match unit {
        Unit::Millimeters | Unit::Centimeters | Unit::Inches => 1.0e-3 / unit.to_si(1.0),
        _ => 1.0,
    };
    let mut arg = if let Ok(text) = value.clone().into_string() {
        let parsed = popkern_fastener::ThreadSpec::parse(&text)
            .map_err(|e| argument_error("spec", &e.to_string()[..]))?;
        ThreadSpecArg {
            major: Scalar::Literal(parsed.major_diameter() * per_mm),
            pitch: Scalar::Literal(parsed.pitch() * per_mm),
            left_handed: false,
            starts: 1,
        }
    } else if let Some(map) = value.read_lock::<rhai::Map>() {
        let get = |key: &str| -> Result<Scalar, Box<EvalAltResult>> {
            map.get(key)
                .map(|v| args::scalar(key, v))
                .transpose()?
                .ok_or_else(|| {
                    argument_error(
                        "spec",
                        &format!(
                            "a thread spec map needs '{key}' (e.g. #{{ major: 6.0, pitch: 1.0 }})"
                        )[..],
                    )
                })
        };
        for key in map.keys() {
            if key.as_str() != "major" && key.as_str() != "pitch" {
                return Err(argument_error(
                    "spec",
                    &format!("thread spec maps take only 'major' and 'pitch' (got '{key}')")[..],
                ));
            }
        }
        ThreadSpecArg {
            major: get("major")?,
            pitch: get("pitch")?,
            left_handed: false,
            starts: 1,
        }
    } else {
        return Err(argument_error(
            "spec",
            "a thread spec is a designation string ('M6', 'M6x0.75', '1/4-20', '#8-32') \
             or a #{ major, pitch } map",
        ));
    };
    if let Some(true) = options.with_if_allowed("left_hand", args::boolean)? {
        arg.left_handed = true;
    }
    if let Some(starts) = options.with_if_allowed("starts", |a, v| {
        v.as_int().map_err(|_| argument_error(a, "an integer"))
    })? {
        if !(1..=8).contains(&starts) {
            return Err(argument_error("starts", "between 1 and 8"));
        }
        arg.starts = starts as u32;
    }
    Ok(arg)
}

/// A blind depth or through-all, from the shared end tokens.
fn hole_depth(value: &Dynamic) -> Result<Option<Scalar>, Box<EvalAltResult>> {
    let end = value.clone().try_cast::<ScriptEnd>().ok_or_else(|| {
        argument_error("depth", "hole depths are blind(distance) or through_all()")
    })?;
    match end {
        ScriptEnd::Extrude(popkern_script::End::Blind(d)) => Ok(Some(d)),
        ScriptEnd::Extrude(popkern_script::End::ThroughAll) => Ok(None),
        _ => Err(argument_error(
            "depth",
            "hole depths are blind(distance) or through_all()",
        )),
    }
}

fn fastener_body(
    m: &mut ScriptModel,
    kind: FastenerKind,
    op: &'static str,
    allowed: &'static [&'static str],
    spec: &Dynamic,
    length: Option<&Dynamic>,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse(op, options, allowed)?;
    let unit = m.borrow_mut().length_unit();
    let thread = thread_spec(spec, &options, unit)?;
    let length = match length {
        Some(v) => args::scalar("length", v)?,
        // A nut's `length` drives its thickness; 0 keeps the
        // standard nominal.
        None => options
            .with("thickness", args::scalar)?
            .unwrap_or(Scalar::Literal(0.0)),
    };
    // Only query keys this verb's option list allows (the Options
    // guard treats an off-list query as a programming error).
    let fastener_options = FastenerOptions {
        thread_length: if allowed.contains(&"thread_length") {
            options.with("thread_length", args::scalar)?
        } else {
            None
        },
        clearance: if allowed.contains(&"clearance") {
            options.with("clearance", args::scalar)?
        } else {
            None
        },
        no_tip_chamfer: allowed.contains(&"tip_chamfer")
            && matches!(options.with("tip_chamfer", args::boolean)?, Some(false)),
        no_corner_chamfer: allowed.contains(&"corner_chamfer")
            && matches!(options.with("corner_chamfer", args::boolean)?, Some(false)),
        no_mouth_chamfer: allowed.contains(&"mouth_chamfer")
            && matches!(options.with("mouth_chamfer", args::boolean)?, Some(false)),
    };
    let body = m
        .borrow_mut()
        .fastener(kind, thread, length, fastener_options, &options.label()?)
        .script()?;
    Ok(ScriptBody::new(m.clone(), body))
}

fn bolt(
    m: &mut ScriptModel,
    spec: Dynamic,
    length: Dynamic,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let parsed = Options::parse("bolt", options, BOLT_OPTIONS)?;
    let kind = parsed
        .with("head", |a, v| {
            args::spelling(
                a,
                v,
                &[
                    ("hex", FastenerKind::HexBolt),
                    ("socket", FastenerKind::SocketCapScrew),
                    ("countersunk", FastenerKind::CountersunkScrew),
                ],
            )
        })?
        .unwrap_or(FastenerKind::HexBolt);
    fastener_body(m, kind, "bolt", BOLT_OPTIONS, &spec, Some(&length), options)
}

fn threaded_rod(
    m: &mut ScriptModel,
    spec: Dynamic,
    length: Dynamic,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    fastener_body(
        m,
        FastenerKind::ThreadedRod,
        "threaded_rod",
        ROD_OPTIONS,
        &spec,
        Some(&length),
        options,
    )
}

fn nut(
    m: &mut ScriptModel,
    spec: Dynamic,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    fastener_body(
        m,
        FastenerKind::HexNut,
        "nut",
        NUT_OPTIONS,
        &spec,
        None,
        options,
    )
}

fn tapped_hole(
    sb: &mut ScriptBody,
    axis: popkern_script::AxisSpec,
    spec: Dynamic,
    depth: Dynamic,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("tapped_hole", options, TAPPED_OPTIONS)?;
    let unit = sb.model.borrow_mut().length_unit();
    let thread = thread_spec(&spec, &options, unit)?;
    let depth = hole_depth(&depth)?;
    let tap = TapHoleOptions {
        clearance: options.with("clearance", args::scalar)?,
        no_entry_chamfer: matches!(options.with("entry_chamfer", args::boolean)?, Some(false)),
    };
    let model = sb.model.clone();
    let body = sb
        .body
        .tapped_hole(
            &mut model.borrow_mut(),
            axis,
            thread,
            depth,
            tap,
            &options.label()?,
        )
        .script()?;
    Ok(ScriptBody::new(model, body))
}

fn clearance_hole(
    sb: &mut ScriptBody,
    axis: popkern_script::AxisSpec,
    spec: Dynamic,
    depth: Dynamic,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("clearance_hole", options, CLEARANCE_OPTIONS)?;
    let unit = sb.model.borrow_mut().length_unit();
    let thread = thread_spec(&spec, &options, unit)?;
    let depth = hole_depth(&depth)?;
    #[derive(Clone, Copy)]
    enum FitSeries {
        Close,
        Normal,
        Loose,
    }
    let fit = match (
        options.with("fit", |a, v| {
            args::spelling(
                a,
                v,
                &[
                    ("close", FitSeries::Close),
                    ("normal", FitSeries::Normal),
                    ("loose", FitSeries::Loose),
                ],
            )
        })?,
        options.with("diameter", args::scalar)?,
    ) {
        (Some(_), Some(_)) => {
            return Err(argument_error(
                "fit",
                "give either a fit series or an explicit diameter, not both",
            ));
        }
        (_, Some(d)) => ClearanceFitArg::Custom(d),
        (Some(FitSeries::Close), None) => ClearanceFitArg::Close,
        (Some(FitSeries::Normal), None) => ClearanceFitArg::Normal,
        (Some(FitSeries::Loose), None) => ClearanceFitArg::Loose,
        (None, None) => ClearanceFitArg::Normal,
    };
    let hole = ClearanceHoleArgs {
        fit,
        countersink: options.with("countersink", args::scalar)?,
        counterbore: options.with("counterbore", args::scalar)?,
        no_entry_chamfer: matches!(options.with("entry_chamfer", args::boolean)?, Some(false)),
    };
    let model = sb.model.clone();
    let body = sb
        .body
        .clearance_hole(
            &mut model.borrow_mut(),
            axis,
            thread,
            depth,
            hole,
            &options.label()?,
        )
        .script()?;
    Ok(ScriptBody::new(model, body))
}

/// Register the fastener verbs.
pub fn register(engine: &mut Engine) {
    engine
        .register_fn(
            "bolt",
            |m: &mut ScriptModel, spec: Dynamic, length: Dynamic| {
                bolt(m, spec, length, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "bolt",
            |m: &mut ScriptModel, spec: Dynamic, length: Dynamic, o: Dynamic| {
                bolt(m, spec, length, &o)
            },
        )
        .register_fn(
            "threaded_rod",
            |m: &mut ScriptModel, spec: Dynamic, length: Dynamic| {
                threaded_rod(m, spec, length, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "threaded_rod",
            |m: &mut ScriptModel, spec: Dynamic, length: Dynamic, o: Dynamic| {
                threaded_rod(m, spec, length, &o)
            },
        )
        .register_fn("nut", |m: &mut ScriptModel, spec: Dynamic| {
            nut(m, spec, &Dynamic::UNIT)
        })
        .register_fn("nut", |m: &mut ScriptModel, spec: Dynamic, o: Dynamic| {
            nut(m, spec, &o)
        })
        .register_fn(
            "tapped_hole",
            |sb: &mut ScriptBody,
             axis: crate::profile_api::ScriptAxis,
             spec: Dynamic,
             depth: Dynamic| {
                tapped_hole(sb, axis.spec(), spec, depth, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "tapped_hole",
            |sb: &mut ScriptBody,
             axis: crate::profile_api::ScriptAxis,
             spec: Dynamic,
             depth: Dynamic,
             o: Dynamic| { tapped_hole(sb, axis.spec(), spec, depth, &o) },
        )
        .register_fn(
            "clearance_hole",
            |sb: &mut ScriptBody,
             axis: crate::profile_api::ScriptAxis,
             spec: Dynamic,
             depth: Dynamic| {
                clearance_hole(sb, axis.spec(), spec, depth, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "clearance_hole",
            |sb: &mut ScriptBody,
             axis: crate::profile_api::ScriptAxis,
             spec: Dynamic,
             depth: Dynamic,
             o: Dynamic| { clearance_hole(sb, axis.spec(), spec, depth, &o) },
        );
}
