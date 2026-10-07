//! Hinge verbs (PK-693): `hinge`, `hinge_component` and `living_hinge`
//! on the model.
//!
//! `m.hinge("butt" | "print_in_place", length, #{…})` returns
//! `#{ leaf_a, leaf_b, pin }` (`pin` is `()` for a print-in-place hinge,
//! and every component not asked for with `only:` is `()`). Dimensions
//! left out follow the length; a screw designation (`screw: "M3"`) adds
//! countersunk holes (`holes:` per leaf, default 2). Hinge dimensions are
//! model units; a screw designation follows the fastener convention
//! (millimetre magnitudes in a default-unit model).

use popkern_history::feature_hinge::{HingeComponentDesc, HingeKindDesc, LivingHingeStyleDesc};
use popkern_script::body_hinge::{HingeHolesArg, HingeOptions, LivingHingeOptions};
use popkern_script::scalar::Scalar;
use rhai::{Dynamic, Engine, EvalAltResult, Map};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::{IntoScriptResult, argument_error};
use crate::model_api::ScriptModel;

const HINGE_OPTIONS: &[&str] = &[
    "label",
    "leaf_width",
    "thickness",
    "knuckles",
    "pin",
    "knuckle_diameter",
    "clearance",
    "gap",
    "angle",
    "holes",
    "screw",
    "countersink",
    "only",
];
const COMPONENT_OPTIONS: &[&str] = &["label"];
const LIVING_OPTIONS: &[&str] = &[
    "label",
    "panel_width",
    "thickness",
    "web",
    "root_radius",
    "slot_pitch",
    "slot_width",
    "bridge",
];

fn kind(text: &str) -> Result<HingeKindDesc, Box<EvalAltResult>> {
    match text {
        "butt" => Ok(HingeKindDesc::Butt),
        "print_in_place" => Ok(HingeKindDesc::PrintInPlace),
        other => Err(argument_error(
            "kind",
            format!("'{other}' is not a hinge kind; kinds are \"butt\" and \"print_in_place\""),
        )),
    }
}

fn component(argument: &str, value: &Dynamic) -> Result<HingeComponentDesc, Box<EvalAltResult>> {
    match args::text(argument, value)?.as_str() {
        "leaf_a" => Ok(HingeComponentDesc::LeafA),
        "leaf_b" => Ok(HingeComponentDesc::LeafB),
        "pin" => Ok(HingeComponentDesc::Pin),
        other => Err(argument_error(
            argument,
            format!(
                "'{other}' is not a hinge component; components are \"leaf_a\", \"leaf_b\" and \"pin\""
            ),
        )),
    }
}

fn hinge(
    m: &mut ScriptModel,
    kind_text: &str,
    length: &Dynamic,
    options: &Dynamic,
) -> Result<Map, Box<EvalAltResult>> {
    let kind = kind(kind_text)?;
    let o = Options::parse("hinge", options, HINGE_OPTIONS)?;
    let unit = m.borrow_mut().length_unit();
    let holes = match o.get("screw") {
        Some(screw) => Some(HingeHolesArg {
            per_leaf: o
                .with("holes", args::scalar)?
                .unwrap_or(Scalar::Literal(2.0)),
            thread: crate::fastener_api::thread_spec(screw, &o, unit)?,
            countersink: o.with("countersink", args::boolean)?.unwrap_or(true),
        }),
        None => {
            if o.get("holes").is_some() || o.get("countersink").is_some() {
                return Err(argument_error(
                    "holes",
                    "screw holes need a screw: add screw: \"M3\" (or #{ major, pitch })",
                ));
            }
            None
        }
    };
    let options = HingeOptions {
        leaf_width: o.with("leaf_width", args::scalar)?,
        thickness: o.with("thickness", args::scalar)?,
        knuckles: o.with("knuckles", args::scalar)?,
        pin: o.with("pin", args::scalar)?,
        knuckle_diameter: o.with("knuckle_diameter", args::scalar)?,
        clearance: o.with("clearance", args::scalar)?,
        gap: o.with("gap", args::scalar)?,
        angle: o.with("angle", args::scalar)?,
        holes,
        only: o.with("only", component)?,
    };
    let length = args::scalar("length", length)?;
    let bodies = m
        .borrow_mut()
        .hinge(kind, length, &options, &o.label()?)
        .script()?;
    let wrap = |b: Option<popkern_script::body::Body>| match b {
        Some(b) => Dynamic::from(ScriptBody::new(m.clone(), b)),
        None => Dynamic::UNIT,
    };
    let mut map = Map::new();
    map.insert("leaf_a".into(), wrap(bodies.leaf_a));
    map.insert("leaf_b".into(), wrap(bodies.leaf_b));
    map.insert("pin".into(), wrap(bodies.pin));
    Ok(map)
}

fn hinge_component(
    m: &mut ScriptModel,
    owner: &str,
    which: &Dynamic,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let o = Options::parse("hinge_component", options, COMPONENT_OPTIONS)?;
    let which = component("component", which)?;
    let body = m
        .borrow_mut()
        .hinge_component(owner, which, &o.label()?)
        .script()?;
    Ok(ScriptBody::new(m.clone(), body))
}

fn living_hinge(
    m: &mut ScriptModel,
    style_text: &str,
    span: &Dynamic,
    length: &Dynamic,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let style = match style_text {
        "straight" => LivingHingeStyleDesc::Straight,
        "slotted" => LivingHingeStyleDesc::Slotted,
        other => {
            return Err(argument_error(
                "style",
                format!(
                    "'{other}' is not a living hinge style; styles are \"straight\" and \"slotted\""
                ),
            ));
        }
    };
    let o = Options::parse("living_hinge", options, LIVING_OPTIONS)?;
    let (own, other): (&[&str], &[&str]) = match style {
        LivingHingeStyleDesc::Straight => (
            &["web", "root_radius"],
            &["slot_pitch", "slot_width", "bridge"],
        ),
        LivingHingeStyleDesc::Slotted => (
            &["slot_pitch", "slot_width", "bridge"],
            &["web", "root_radius"],
        ),
    };
    if let Some(key) = other.iter().copied().find(|k| o.get(k).is_some()) {
        return Err(argument_error(
            key,
            format!(
                "a {style_text} living hinge takes {} (not {key})",
                own.join(", ")
            ),
        ));
    }
    let options = LivingHingeOptions {
        panel_width: o.with("panel_width", args::scalar)?,
        thickness: o.with("thickness", args::scalar)?,
        web: o.with("web", args::scalar)?,
        root_radius: o.with("root_radius", args::scalar)?,
        slot_pitch: o.with("slot_pitch", args::scalar)?,
        slot_width: o.with("slot_width", args::scalar)?,
        bridge: o.with("bridge", args::scalar)?,
    };
    let span = args::scalar("span", span)?;
    let length = args::scalar("length", length)?;
    let body = m
        .borrow_mut()
        .living_hinge(style, span, length, &options, &o.label()?)
        .script()?;
    Ok(ScriptBody::new(m.clone(), body))
}

/// Register the hinge verbs.
pub fn register(engine: &mut Engine) {
    engine
        .register_fn(
            "hinge",
            |m: &mut ScriptModel, kind: &str, length: Dynamic| {
                hinge(m, kind, &length, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "hinge",
            |m: &mut ScriptModel, kind: &str, length: Dynamic, o: Dynamic| {
                hinge(m, kind, &length, &o)
            },
        )
        .register_fn(
            "hinge_component",
            |m: &mut ScriptModel, owner: &str, which: Dynamic| {
                hinge_component(m, owner, &which, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "hinge_component",
            |m: &mut ScriptModel, owner: &str, which: Dynamic, o: Dynamic| {
                hinge_component(m, owner, &which, &o)
            },
        )
        .register_fn(
            "living_hinge",
            |m: &mut ScriptModel, style: &str, span: Dynamic, length: Dynamic| {
                living_hinge(m, style, &span, &length, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "living_hinge",
            |m: &mut ScriptModel, style: &str, span: Dynamic, length: Dynamic, o: Dynamic| {
                living_hinge(m, style, &span, &length, &o)
            },
        );
}
