//! Seam hinge verbs (PK-702): a hinge grown on two existing bodies.
//!
//! ```rhai
//! let h = m.hinge_bodies("print_in_place", box, lid, axis(edges().at_height(20.0).direction("x")), #{ knuckles: 5 });
//! // h.a, h.b (the bodies, now knuckled), h.pin (a butt hinge's pin, else ())
//! ```
//!
//! With `axis(edges()…)` (an edge of the first body) the seam runs
//! along that edge and the hinge is `length` long (default 80% of the
//! edge), centred on it; with `axis_at(origin, direction)` it starts at
//! `origin` and `length` is required. Show-as-code spells the parts:
//! `a.seam_hinge(kind, axis_at(…), length, #{ outward, up, … })` on body
//! A, `b.seam_hinge_component(owner, "leaf_b")`, and
//! `m.seam_hinge_component(owner, "pin")`.

use popkern_core::geometry::{Point3, Vec3};
use popkern_history::feature_hinge::{HingeComponentDesc, HingeKindDesc, SeamDesc};
use popkern_script::AxisSpec;
use popkern_script::body_seam_hinge::SeamHingeOptions;
use popkern_script::scalar::Scalar;
use rhai::{Array, Dynamic, Engine, EvalAltResult, Map};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::{IntoScriptResult, argument_error};
use crate::model_api::ScriptModel;
use crate::profile_api::ScriptAxis;

const SIZE_KEYS: [&str; 6] = [
    "knuckles",
    "pin",
    "wall",
    "knuckle_diameter",
    "clearance",
    "gap",
];
const BODIES_OPTIONS: &[&str] = &[
    "label",
    "length",
    "knuckles",
    "pin",
    "wall",
    "knuckle_diameter",
    "clearance",
    "gap",
];
const SEAM_OPTIONS: &[&str] = &[
    "label",
    "outward",
    "up",
    "knuckles",
    "pin",
    "wall",
    "knuckle_diameter",
    "clearance",
    "gap",
];
const LABEL_ONLY: &[&str] = &["label"];

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

fn sizes(o: &Options) -> Result<SeamHingeOptions, Box<EvalAltResult>> {
    debug_assert!(SIZE_KEYS.len() == 6);
    Ok(SeamHingeOptions {
        knuckles: o.with("knuckles", args::scalar)?,
        pin: o.with("pin", args::scalar)?,
        wall: o.with("wall", args::scalar)?,
        knuckle_diameter: o.with("knuckle_diameter", args::scalar)?,
        clearance: o.with("clearance", args::scalar)?,
        gap: o.with("gap", args::scalar)?,
    })
}

fn same_model(m: &ScriptModel, b: &ScriptBody) -> Result<(), Box<EvalAltResult>> {
    if !m.same_document(&b.model) {
        return Err(argument_error("body", "this body belongs to another model"));
    }
    Ok(())
}

fn hinge_bodies(
    m: &mut ScriptModel,
    kind_text: &str,
    a: ScriptBody,
    b: ScriptBody,
    axis: ScriptAxis,
    options: &Dynamic,
) -> Result<Map, Box<EvalAltResult>> {
    let kind = kind(kind_text)?;
    same_model(m, &a)?;
    same_model(m, &b)?;
    let o = Options::parse("hinge_bodies", options, BODIES_OPTIONS)?;
    let length = o.with("length", args::scalar)?;
    let (origin, direction, length) = match axis.spec() {
        AxisSpec::Literal { origin, direction } => {
            let length = length.ok_or_else(|| {
                argument_error("length", "a seam given by axis_at(...) needs length: …")
            })?;
            (
                Point3::new(origin[0], origin[1], origin[2]),
                Vec3::new(direction[0], direction[1], direction[2]),
                length,
            )
        }
        AxisSpec::Edge(selector) => {
            let (p, q) = m.borrow_mut().seam_edge(a.body, &selector).script()?;
            let span = q - p;
            let edge = span.length();
            let length = length.unwrap_or(Scalar::Literal(0.8 * edge));
            // Centre a literal length on the edge; a parameter starts
            // where a literal of its current value would.
            let l = match &length {
                Scalar::Literal(v) => *v,
                _ => 0.8 * edge,
            };
            let unit = span * (1.0 / edge.max(f64::MIN_POSITIVE));
            (p + unit * (0.5 * (edge - l)), span, length)
        }
    };
    let bodies = m
        .borrow_mut()
        .hinge_bodies(
            kind,
            a.body,
            b.body,
            origin,
            direction,
            length,
            &sizes(&o)?,
            &o.label()?,
        )
        .script()?;
    let wrap = |x: Option<popkern_script::body::Body>| match x {
        Some(x) => Dynamic::from(ScriptBody::new(m.clone(), x)),
        None => Dynamic::UNIT,
    };
    let mut map = Map::new();
    map.insert("a".into(), wrap(bodies.leaf_a));
    map.insert("b".into(), wrap(bodies.leaf_b));
    map.insert("pin".into(), wrap(bodies.pin));
    Ok(map)
}

fn vec3(argument: &str, value: &Dynamic) -> Result<[f64; 3], Box<EvalAltResult>> {
    let array = value
        .clone()
        .try_cast::<Array>()
        .filter(|a| a.len() == 3)
        .ok_or_else(|| argument_error(argument, "expected [x, y, z]"))?;
    let mut out = [0.0; 3];
    for (k, v) in array.iter().enumerate() {
        out[k] = args::float(argument, v)?;
    }
    Ok(out)
}

/// `a.seam_hinge(kind, axis_at(origin, direction), length, #{ outward, up, … })`:
/// body A's knuckles with the seam spelled out (show-as-code).
fn seam_hinge(
    a: &mut ScriptBody,
    kind_text: &str,
    axis: ScriptAxis,
    length: Dynamic,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let kind = kind(kind_text)?;
    let o = Options::parse("seam_hinge", options, SEAM_OPTIONS)?;
    let AxisSpec::Literal { origin, direction } = axis.spec() else {
        return Err(argument_error(
            "axis",
            "seam_hinge takes axis_at(origin, direction)",
        ));
    };
    let outward = o
        .with("outward", vec3)?
        .ok_or_else(|| argument_error("outward", "seam_hinge needs outward: [x, y, z]"))?;
    let up = o
        .with("up", vec3)?
        .ok_or_else(|| argument_error("up", "seam_hinge needs up: [x, y, z]"))?;
    let seam = SeamDesc {
        origin,
        axis: direction,
        outward,
        up,
    };
    let length = args::scalar("length", &length)?;
    let body = a
        .model
        .borrow_mut()
        .seam_hinge(kind, a.body, seam, length, &sizes(&o)?, &o.label()?)
        .script()?;
    Ok(ScriptBody::new(a.model.clone(), body))
}

fn component(text: &str) -> Result<HingeComponentDesc, Box<EvalAltResult>> {
    match text {
        "leaf_b" => Ok(HingeComponentDesc::LeafB),
        "pin" => Ok(HingeComponentDesc::Pin),
        other => Err(argument_error(
            "component",
            format!(
                "'{other}': a seam hinge's other components are \"leaf_b\" (on a body) and \"pin\""
            ),
        )),
    }
}

/// The two forms of `seam_hinge_component`: on a body (`leaf_b`) or the
/// model (`pin`).
fn seam_component(
    model: &ScriptModel,
    body: Option<popkern_script::body::Body>,
    owner: &str,
    which: &str,
    options: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let o = Options::parse("seam_hinge_component", options, LABEL_ONLY)?;
    let made = model
        .borrow_mut()
        .seam_hinge_component(owner, component(which)?, body, &o.label()?)
        .script()?;
    Ok(ScriptBody::new(model.clone(), made))
}

/// Register the seam hinge verbs.
pub fn register(engine: &mut Engine) {
    engine
        .register_fn(
            "hinge_bodies",
            |m: &mut ScriptModel, k: &str, a: ScriptBody, b: ScriptBody, axis: ScriptAxis| {
                hinge_bodies(m, k, a, b, axis, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "hinge_bodies",
            |m: &mut ScriptModel,
             k: &str,
             a: ScriptBody,
             b: ScriptBody,
             axis: ScriptAxis,
             o: Dynamic| { hinge_bodies(m, k, a, b, axis, &o) },
        )
        .register_fn(
            "seam_hinge",
            |a: &mut ScriptBody, k: &str, axis: ScriptAxis, length: Dynamic, o: Dynamic| {
                seam_hinge(a, k, axis, length, &o)
            },
        )
        .register_fn(
            "seam_hinge_component",
            |b: &mut ScriptBody, owner: &str, which: &str| {
                seam_component(&b.model.clone(), Some(b.body), owner, which, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "seam_hinge_component",
            |b: &mut ScriptBody, owner: &str, which: &str, o: Dynamic| {
                seam_component(&b.model.clone(), Some(b.body), owner, which, &o)
            },
        )
        .register_fn(
            "seam_hinge_component",
            |m: &mut ScriptModel, owner: &str, which: &str| {
                seam_component(m, None, owner, which, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "seam_hinge_component",
            |m: &mut ScriptModel, owner: &str, which: &str, o: Dynamic| {
                seam_component(m, None, owner, which, &o)
            },
        );
}
