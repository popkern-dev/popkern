//! Mates in the Rhai skin (design §8.2, PK-631):
//!
//! ```rhai
//! a.ground(ring);
//! a.concentric(planet.datum("z"), carrier.face("pin_0"));
//! a.coincident(planet.datum("xy"), ring.datum("xy"), #{ aligned: true });
//! a.angle(planet.datum("zx"), ring.datum("zx"), 40.0, #{ about: ring.datum("z") });
//! ```
//!
//! An entity is `instance.face(label_or_spelled_name)`,
//! `instance.face(faces()…)` (a selection that must match exactly one
//! face of the part), or `instance.datum("xy" | "yz" | "zx" | "x" | "y" |
//! "z")`; `a.occurrence("stage/planet_1")` names an occurrence inside a
//! sub-assembly. Each mate solves at once and returns the solve's
//! report; a conflicting mate raises and moves nothing.

use popkern_assembly::{Datum, Mate, MateEntity, MateKind, MateRef, Parent};
use popkern_script::FaceSelector;
use popkern_script::scalar::Scalar;
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::assembly_api::{IntoAssemblyResult, ScriptAssembly, ScriptInstance};
use crate::error::argument_error;

/// An occurrence a mate can reference (a root instance or a path).
#[derive(Clone)]
pub struct ScriptOccurrence {
    owner: ScriptAssembly,
    path: String,
}

impl std::fmt::Debug for ScriptOccurrence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Occurrence({})", self.path)
    }
}

/// An entity reference value.
#[derive(Debug, Clone)]
pub struct ScriptMateRef(pub MateRef);

fn occurrence_of(i: &ScriptInstance) -> Result<ScriptOccurrence, Box<EvalAltResult>> {
    if i.at.parent != Parent::Root {
        return Err(argument_error(
            "instance",
            "an instance inside a sub-assembly is named by its path: a.occurrence(\"stage/planet_1\")",
        ));
    }
    let name = i
        .owner
        .borrow_mut()
        .assembly
        .instance(i.at)
        .assembly()?
        .name
        .clone();
    Ok(ScriptOccurrence {
        owner: i.owner.clone(),
        path: name,
    })
}

fn datum(o: ScriptOccurrence, which: &str) -> Result<ScriptMateRef, Box<EvalAltResult>> {
    let datum = Datum::parse(which).ok_or_else(|| {
        argument_error(
            "datum",
            format!("expected xy, yz, zx, x, y or z, got '{which}'"),
        )
    })?;
    Ok(ScriptMateRef(MateRef {
        instance: o.path,
        entity: MateEntity::Datum { datum },
    }))
}

/// `instance.axis([o], [d])` / `instance.plane([o], [n])`: a line or a
/// plane in the part's own frame (PK-712).
fn literal(
    o: ScriptOccurrence,
    point: &Dynamic,
    direction: &Dynamic,
    axis: bool,
) -> Result<ScriptMateRef, Box<EvalAltResult>> {
    let triple = |name: &str, v: &Dynamic| -> Result<[f64; 3], Box<EvalAltResult>> {
        let a = v
            .clone()
            .try_cast::<rhai::Array>()
            .filter(|a| a.len() == 3)
            .ok_or_else(|| argument_error(name, "expected [x, y, z]"))?;
        Ok([
            args::float(name, &a[0])?,
            args::float(name, &a[1])?,
            args::float(name, &a[2])?,
        ])
    };
    let origin = triple("origin", point)?;
    let dir = triple("direction", direction)?;
    let entity = if axis {
        MateEntity::Axis {
            origin,
            direction: dir,
        }
    } else {
        MateEntity::Plane {
            origin,
            normal: dir,
        }
    };
    Ok(ScriptMateRef(MateRef {
        instance: o.path,
        entity,
    }))
}

fn face_named(o: ScriptOccurrence, name: &str) -> ScriptMateRef {
    ScriptMateRef(MateRef {
        instance: o.path,
        entity: MateEntity::Face {
            name: name.to_string(),
        },
    })
}

/// The one face of the occurrence's part `selector` matches, by its
/// spelled name.
fn face_selected(
    o: ScriptOccurrence,
    selector: FaceSelector,
) -> Result<ScriptMateRef, Box<EvalAltResult>> {
    let mut state = o.owner.borrow_mut();
    let rows = state.assembly.tree();
    let row = rows
        .iter()
        .find(|r| r.path.to_string() == o.path)
        .ok_or_else(|| argument_error("instance", format!("no occurrence '{}'", o.path)))?;
    let definition = row.definition;
    let model = state.assembly.part_model_mut(definition).assembly()?;
    let mut names = Vec::new();
    for body in model.alive_bodies() {
        let brep = model.brep_of(body).map_err(crate::error::script_error)?;
        let snapshot = model
            .resolve_faces_in(&brep, &selector)
            .map_err(crate::error::script_error)?;
        names.extend(snapshot.names);
    }
    if names.len() != 1 {
        return Err(argument_error(
            "face",
            format!(
                "the selection must match exactly one face of '{}', matched {}",
                o.path,
                names.len()
            ),
        ));
    }
    let spelled = model
        .spell_entity(&names[0])
        .map_err(crate::error::script_error)?;
    drop(state);
    Ok(face_named(o, &spelled))
}

fn add(
    a: &mut ScriptAssembly,
    kind: MateKind,
    x: ScriptMateRef,
    y: ScriptMateRef,
) -> Result<Dynamic, Box<EvalAltResult>> {
    let report = a
        .borrow_mut()
        .assembly
        .add_mate(Mate {
            kind,
            a: x.0,
            b: y.0,
            drive: None,
        })
        .assembly()?;
    rhai::serde::to_dynamic(&report).map_err(|e| argument_error("mate", e.to_string()))
}

fn aligned(name: &'static str, options: &Dynamic) -> Result<bool, Box<EvalAltResult>> {
    let options = Options::parse(name, options, &["aligned"])?;
    Ok(options.with("aligned", args::boolean)?.unwrap_or(false))
}

fn angle(
    a: &mut ScriptAssembly,
    x: ScriptMateRef,
    y: ScriptMateRef,
    degrees: Dynamic,
    options: Dynamic,
) -> Result<Dynamic, Box<EvalAltResult>> {
    let about = if options.is_unit() {
        None
    } else {
        let map = options
            .clone()
            .try_cast::<rhai::Map>()
            .ok_or_else(|| argument_error("options", "expected #{ about: entity }"))?;
        for key in map.keys() {
            if !["about", "scale", "offset"].contains(&key.as_str()) {
                return Err(argument_error(
                    key,
                    "unknown option (angle takes: about, scale, offset)",
                ));
            }
        }
        match map.get("about") {
            Some(v) => Some(
                v.clone()
                    .try_cast::<ScriptMateRef>()
                    .ok_or_else(|| {
                        argument_error(
                            "about",
                            "expected an entity: instance.datum(\"z\") or instance.face(…)",
                        )
                    })?
                    .0,
            ),
            None => None,
        }
    };
    let (scale, offset) = drive_map(&options)?;
    match args::scalar("degrees", &degrees)? {
        Scalar::Literal(v) => {
            no_drive(scale, offset)?;
            add(a, MateKind::Angle { degrees: v, about }, x, y)
        }
        // `p("name")`: the assembly parameter drives the angle (PK-701),
        // as offset + scale × parameter.
        Scalar::Param(param) => {
            let v = offset + scale * param_value(a, &param)?;
            let report = add(a, MateKind::Angle { degrees: v, about }, x, y)?;
            drive_last(a, &param, scale, offset)?;
            Ok(report)
        }
        Scalar::Expr(_) => Err(argument_error(
            "degrees",
            "a mate angle is a number or p(\"assembly parameter\")",
        )),
    }
}

/// A distance mate, its value a number or `p("assembly parameter")`.
fn distance(
    a: &mut ScriptAssembly,
    x: ScriptMateRef,
    y: ScriptMateRef,
    d: &Dynamic,
    options: &Dynamic,
) -> Result<Dynamic, Box<EvalAltResult>> {
    let o = Options::parse("distance", options, &["aligned", "scale", "offset"])?;
    let aligned = o.with("aligned", args::boolean)?.unwrap_or(false);
    let scale = o.with("scale", args::float)?.unwrap_or(1.0);
    let offset = o.with("offset", args::float)?.unwrap_or(0.0);
    match args::scalar("distance", d)? {
        Scalar::Literal(distance) => {
            no_drive(scale, offset)?;
            add(a, MateKind::Distance { distance, aligned }, x, y)
        }
        Scalar::Param(param) => {
            let distance = offset + scale * param_value(a, &param)?;
            let report = add(a, MateKind::Distance { distance, aligned }, x, y)?;
            drive_last(a, &param, scale, offset)?;
            Ok(report)
        }
        Scalar::Expr(_) => Err(argument_error(
            "distance",
            "a mate distance is a number or p(\"assembly parameter\")",
        )),
    }
}

/// The display value of assembly parameter `name`.
fn param_value(a: &ScriptAssembly, name: &str) -> Result<f64, Box<EvalAltResult>> {
    a.borrow_mut()
        .assembly
        .params()
        .iter()
        .find(|p| p.name == name)
        .map(|p| p.display_value())
        .ok_or_else(|| argument_error("p", format!("no assembly parameter '{name}'")))
}

/// Bind the mate just added to assembly parameter `name`.
fn drive_last(
    a: &mut ScriptAssembly,
    name: &str,
    scale: f64,
    offset: f64,
) -> Result<(), Box<EvalAltResult>> {
    let mut state = a.borrow_mut();
    let index = state.assembly.mates().len().saturating_sub(1);
    state
        .assembly
        .drive_mate(name, index, scale, offset)
        .assembly()
}

/// `scale` / `offset` from an angle's option map (defaults 1, 0).
fn drive_map(options: &Dynamic) -> Result<(f64, f64), Box<EvalAltResult>> {
    let Some(map) = options.clone().try_cast::<rhai::Map>() else {
        return Ok((1.0, 0.0));
    };
    let get = |key: &str, default: f64| match map.get(key) {
        Some(v) => args::float(key, v),
        None => Ok(default),
    };
    Ok((get("scale", 1.0)?, get("offset", 0.0)?))
}

/// `scale` / `offset` only make sense on a mate a parameter drives.
fn no_drive(scale: f64, offset: f64) -> Result<(), Box<EvalAltResult>> {
    if scale != 1.0 || offset != 0.0 {
        return Err(argument_error(
            "scale",
            "scale and offset map a driving parameter: give the value as p(\"name\")",
        ));
    }
    Ok(())
}

/// Register the mate verbs.
pub fn register(engine: &mut Engine) {
    engine
        .register_type_with_name::<ScriptOccurrence>("Occurrence")
        .register_type_with_name::<ScriptMateRef>("MateEntity")
        .register_fn("occurrence", |a: &mut ScriptAssembly, path: &str| {
            ScriptOccurrence {
                owner: a.clone(),
                path: path.to_string(),
            }
        })
        .register_fn("datum", |i: &mut ScriptInstance, d: &str| {
            datum(occurrence_of(i)?, d)
        })
        .register_fn("datum", |o: &mut ScriptOccurrence, d: &str| {
            datum(o.clone(), d)
        })
        .register_fn("face", |i: &mut ScriptInstance, n: &str| {
            Ok::<_, Box<EvalAltResult>>(face_named(occurrence_of(i)?, n))
        })
        .register_fn("face", |o: &mut ScriptOccurrence, n: &str| {
            face_named(o.clone(), n)
        })
        .register_fn("face", |i: &mut ScriptInstance, s: FaceSelector| {
            face_selected(occurrence_of(i)?, s)
        })
        .register_fn("face", |o: &mut ScriptOccurrence, s: FaceSelector| {
            face_selected(o.clone(), s)
        })
        .register_fn("axis", |i: &mut ScriptInstance, p: Dynamic, d: Dynamic| {
            literal(occurrence_of(i)?, &p, &d, true)
        })
        .register_fn(
            "axis",
            |o: &mut ScriptOccurrence, p: Dynamic, d: Dynamic| literal(o.clone(), &p, &d, true),
        )
        .register_fn("plane", |i: &mut ScriptInstance, p: Dynamic, n: Dynamic| {
            literal(occurrence_of(i)?, &p, &n, false)
        })
        .register_fn(
            "plane",
            |o: &mut ScriptOccurrence, p: Dynamic, n: Dynamic| literal(o.clone(), &p, &n, false),
        )
        .register_fn("ground", |a: &mut ScriptAssembly, i: ScriptInstance| {
            let name = occurrence_of(&i)?.path;
            a.borrow_mut().assembly.set_ground(&name).assembly()
        })
        .register_fn("ground", |a: &mut ScriptAssembly, o: ScriptOccurrence| {
            a.borrow_mut().assembly.set_ground(&o.path).assembly()
        })
        .register_fn(
            "concentric",
            |a: &mut ScriptAssembly, x: ScriptMateRef, y: ScriptMateRef| {
                add(a, MateKind::Concentric, x, y)
            },
        )
        .register_fn(
            "coincident",
            |a: &mut ScriptAssembly, x: ScriptMateRef, y: ScriptMateRef| {
                add(a, MateKind::Coincident { aligned: false }, x, y)
            },
        )
        .register_fn(
            "coincident",
            |a: &mut ScriptAssembly, x: ScriptMateRef, y: ScriptMateRef, o: Dynamic| {
                add(
                    a,
                    MateKind::Coincident {
                        aligned: aligned("coincident", &o)?,
                    },
                    x,
                    y,
                )
            },
        )
        .register_fn(
            "distance",
            |a: &mut ScriptAssembly, x: ScriptMateRef, y: ScriptMateRef, d: Dynamic| {
                distance(a, x, y, &d, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "distance",
            |a: &mut ScriptAssembly, x: ScriptMateRef, y: ScriptMateRef, d: Dynamic, o: Dynamic| {
                distance(a, x, y, &d, &o)
            },
        )
        .register_fn(
            "angle",
            |a: &mut ScriptAssembly, x: ScriptMateRef, y: ScriptMateRef, d: Dynamic| {
                angle(a, x, y, d, Dynamic::UNIT)
            },
        )
        .register_fn("angle", angle)
        .register_fn("solve", |a: &mut ScriptAssembly| {
            let report = a.borrow_mut().assembly.solve_mates().assembly()?;
            rhai::serde::to_dynamic(&report).map_err(|e| argument_error("solve", e.to_string()))
        });
}
