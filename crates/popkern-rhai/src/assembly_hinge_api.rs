//! Assembly hinge verbs (PK-701):
//!
//! ```rhai
//! let a = assembly("box");
//! let h = a.hinge("lid", "print_in_place", 60.0, #{ angle: 110.0 });
//! // h.leaf_a, h.leaf_b (instances), h.pin (instance or ()), h.angle == "lid_angle"
//! a.limits("lid_angle", 0.0, 250.0);
//! let swing = a.swing("lid_angle", #{ from: 0.0, to: 250.0, step: 5.0 });
//! ```

use popkern_history::feature_hinge::{HingeComponentDesc, HingeKindDesc};
use popkern_script::body_hinge::{HingeHolesArg, HingeOptions};
use popkern_script::scalar::Scalar;
use rhai::{Dynamic, Engine, EvalAltResult, Map};

use crate::args::{self, Options};
use crate::assembly_api::{IntoAssemblyResult, ScriptAssembly, ScriptInstance};
use crate::error::argument_error;

const HINGE_OPTIONS: &[&str] = &[
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
];
const SWING_OPTIONS: &[&str] = &["from", "to", "step", "resolution", "clearance", "exclude"];

fn hinge(
    a: &mut ScriptAssembly,
    name: &str,
    kind: &str,
    length: &Dynamic,
    options: &Dynamic,
) -> Result<Map, Box<EvalAltResult>> {
    if a.parent != popkern_assembly::Parent::Root {
        return Err(argument_error(
            "hinge",
            "an assembly hinge is placed at the top level of the assembly",
        ));
    }
    let kind = match kind {
        "butt" => HingeKindDesc::Butt,
        "print_in_place" => HingeKindDesc::PrintInPlace,
        other => {
            return Err(argument_error(
                "kind",
                format!("'{other}' is not a hinge kind; kinds are \"butt\" and \"print_in_place\""),
            ));
        }
    };
    let o = Options::parse("hinge", options, HINGE_OPTIONS)?;
    let holes = match o.get("screw") {
        Some(screw) => Some(HingeHolesArg {
            per_leaf: o
                .with("holes", args::scalar)?
                .unwrap_or(Scalar::Literal(2.0)),
            thread: crate::fastener_api::thread_spec(
                screw,
                &o,
                popkern_history::unit::Unit::Meters,
            )?,
            countersink: o.with("countersink", args::boolean)?.unwrap_or(true),
        }),
        None => None,
    };
    let options = HingeOptions {
        leaf_width: o.with("leaf_width", args::scalar)?,
        thickness: o.with("thickness", args::scalar)?,
        knuckles: o.with("knuckles", args::scalar)?,
        pin: o.with("pin", args::scalar)?,
        knuckle_diameter: o.with("knuckle_diameter", args::scalar)?,
        clearance: o.with("clearance", args::scalar)?,
        gap: o.with("gap", args::scalar)?,
        angle: o.with("angle", args::float)?.map(Scalar::Literal),
        holes,
        only: None,
    };
    let length = args::scalar("length", length)?;
    let parts = a
        .borrow_mut()
        .assembly
        .add_hinge(name, kind, length, &options)
        .assembly()?;
    let mut map = Map::new();
    map.insert("pin".into(), Dynamic::UNIT);
    for (component, _, _, at) in parts.parts {
        let key = match component {
            HingeComponentDesc::LeafA => "leaf_a",
            HingeComponentDesc::LeafB => "leaf_b",
            HingeComponentDesc::Pin => "pin",
        };
        map.insert(
            key.into(),
            Dynamic::from(ScriptInstance {
                owner: a.clone(),
                at,
            }),
        );
    }
    map.insert("angle".into(), Dynamic::from(parts.angle));
    Ok(map)
}

fn swing(
    a: &mut ScriptAssembly,
    param: &str,
    options: &Dynamic,
) -> Result<Dynamic, Box<EvalAltResult>> {
    let o = Options::parse("swing", options, SWING_OPTIONS)?;
    let (from, to) = {
        let state = a.borrow_mut();
        let p = state
            .assembly
            .params()
            .iter()
            .find(|p| p.name == param)
            .cloned()
            .ok_or_else(|| argument_error("param", format!("no assembly parameter '{param}'")))?;
        let [lo, hi] = p.bounds.unwrap_or([p.display_value(), p.display_value()]);
        (lo, hi)
    };
    let from = o.with("from", args::float)?.unwrap_or(from);
    let to = o.with("to", args::float)?.unwrap_or(to);
    let step = o.with("step", args::float)?.unwrap_or(5.0);
    let resolution = o.with("resolution", args::float)?.unwrap_or(0.1);
    let clearance = o.with("clearance", args::float)?;
    let exclude = o
        .with("exclude", |a, v| {
            let pairs = v
                .clone()
                .try_cast::<rhai::Array>()
                .ok_or_else(|| argument_error(a, "expected [[prefix, prefix], ...]"))?;
            pairs
                .into_iter()
                .map(|p| {
                    let p = p
                        .try_cast::<rhai::Array>()
                        .filter(|p| p.len() == 2)
                        .ok_or_else(|| argument_error(a, "each entry is [prefix, prefix]"))?;
                    let s = |d: &Dynamic| {
                        d.clone()
                            .into_immutable_string()
                            .map(|s| s.to_string())
                            .map_err(|_| argument_error(a, "prefixes are strings"))
                    };
                    Ok([s(&p[0])?, s(&p[1])?])
                })
                .collect::<Result<Vec<_>, Box<EvalAltResult>>>()
        })?
        .unwrap_or_default();
    let report = a
        .borrow_mut()
        .assembly
        .swing_with(
            param,
            &popkern_assembly::SwingOptions {
                from,
                to,
                step,
                resolution,
                clearance,
                exclude,
            },
        )
        .assembly()?;
    rhai::serde::to_dynamic(&report).map_err(|e| argument_error("swing", e.to_string()))
}

fn limits(
    a: &mut ScriptAssembly,
    param: &str,
    min: Dynamic,
    max: Dynamic,
) -> Result<(), Box<EvalAltResult>> {
    let (min, max) = (args::float("min", &min)?, args::float("max", &max)?);
    a.borrow_mut()
        .assembly
        .set_param_bounds(param, min, max)
        .assembly()
}

const SEAM_OPTIONS: &[&str] = &[
    "length",
    "knuckles",
    "pin",
    "wall",
    "knuckle_diameter",
    "clearance",
    "gap",
    "angle",
];

/// `a.hinge_bodies(name, kind, box, lid, axis_at(origin, direction),
/// #{ length, … })` (PK-712): two placed parts hinged along their seam,
/// the lid opening by `<name>_angle`.
fn hinge_bodies(
    a: &mut ScriptAssembly,
    name: &str,
    kind: &str,
    box_: ScriptInstance,
    lid: ScriptInstance,
    axis: crate::profile_api::ScriptAxis,
    options: &Dynamic,
) -> Result<Map, Box<EvalAltResult>> {
    let kind = match kind {
        "butt" => HingeKindDesc::Butt,
        "print_in_place" => HingeKindDesc::PrintInPlace,
        other => {
            return Err(argument_error(
                "kind",
                format!("'{other}' is not a hinge kind; kinds are \"butt\" and \"print_in_place\""),
            ));
        }
    };
    let o = Options::parse("hinge_bodies", options, SEAM_OPTIONS)?;
    let popkern_script::AxisSpec::Literal { origin, direction } = axis.spec() else {
        return Err(argument_error(
            "axis",
            "an assembly seam hinge takes axis_at(origin, direction) in assembly coordinates",
        ));
    };
    let length = o
        .with("length", args::scalar)?
        .ok_or_else(|| argument_error("length", "an assembly seam hinge needs length: …"))?;
    let sizes = popkern_script::body_seam_hinge::SeamHingeOptions {
        knuckles: o.with("knuckles", args::scalar)?,
        pin: o.with("pin", args::scalar)?,
        wall: o.with("wall", args::scalar)?,
        knuckle_diameter: o.with("knuckle_diameter", args::scalar)?,
        clearance: o.with("clearance", args::scalar)?,
        gap: o.with("gap", args::scalar)?,
    };
    let start = o.with("angle", args::float)?.unwrap_or(0.0);
    let parts = a
        .borrow_mut()
        .assembly
        .add_seam_hinge(
            name,
            kind,
            box_.at,
            lid.at,
            popkern_core::geometry::Point3::new(origin[0], origin[1], origin[2]),
            popkern_core::geometry::Vec3::new(direction[0], direction[1], direction[2]),
            length,
            &sizes,
            start,
        )
        .assembly()?;
    let mut map = Map::new();
    map.insert("angle".into(), Dynamic::from(parts.angle));
    map.insert(
        "pin".into(),
        match parts.pin {
            Some(at) => Dynamic::from(ScriptInstance {
                owner: a.clone(),
                at,
            }),
            None => Dynamic::UNIT,
        },
    );
    Ok(map)
}

/// Register the assembly hinge verbs.
pub fn register(engine: &mut Engine) {
    engine
        .register_fn(
            "hinge",
            |a: &mut ScriptAssembly, name: &str, kind: &str, length: Dynamic| {
                hinge(a, name, kind, &length, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "hinge",
            |a: &mut ScriptAssembly, name: &str, kind: &str, length: Dynamic, o: Dynamic| {
                hinge(a, name, kind, &length, &o)
            },
        )
        .register_fn("swing", |a: &mut ScriptAssembly, param: &str| {
            swing(a, param, &Dynamic::UNIT)
        })
        .register_fn(
            "swing",
            |a: &mut ScriptAssembly, param: &str, o: Dynamic| swing(a, param, &o),
        )
        .register_fn("limits", limits)
        .register_fn(
            "hinge_bodies",
            |a: &mut ScriptAssembly,
             name: &str,
             kind: &str,
             b: ScriptInstance,
             l: ScriptInstance,
             axis: crate::profile_api::ScriptAxis,
             o: Dynamic| { hinge_bodies(a, name, kind, b, l, axis, &o) },
        );
}
