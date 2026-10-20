//! Revolute joints in the Rhai skin (design
//! `docs/architecture/six-axis-arm.md` §4.3):
//!
//! ```rhai
//! // link_2 turns about the shoulder axis; a slider "j2" drives it.
//! a.revolute("j2", link_1.axis([0,0,0],[0,1,0]), link_2.axis([0,0,0],[0,1,0]),
//!            #{ limits: [-360.0, 360.0] });
//! // The rotor follows j2 at −29×, measured against the housing.
//! a.revolute("j2_rotor", link_1.axis(…), rotor.datum("z"),
//!            #{ follow: "j2", ratio: -29.0 });
//! ```
//!
//! Options: `zero: [entity, entity]` (directions whose signed angle
//! about the first axis is the joint angle; default perpendicular to
//! each axis), `angle` (start, degrees; default the current angle),
//! `limits: [min, max]`, `offset` (axial), `follow` + `ratio` (+
//! `phase`). The joint solves at once and returns its report.

use popkern_assembly::RevoluteOptions;
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::assembly_api::{IntoAssemblyResult, ScriptAssembly};
use crate::assembly_mate_api::ScriptMateRef;
use crate::error::argument_error;

const REVOLUTE_OPTIONS: &[&str] = &[
    "zero", "angle", "limits", "offset", "follow", "ratio", "phase",
];

fn entity(name: &str, v: &Dynamic) -> Result<popkern_assembly::MateRef, Box<EvalAltResult>> {
    v.clone()
        .try_cast::<ScriptMateRef>()
        .map(|r| r.0)
        .ok_or_else(|| {
            argument_error(
                name,
                "expected an entity: instance.axis(…), .datum(…) or .face(…)",
            )
        })
}

fn pair(name: &str, v: &Dynamic) -> Result<rhai::Array, Box<EvalAltResult>> {
    v.clone()
        .try_cast::<rhai::Array>()
        .filter(|a| a.len() == 2)
        .ok_or_else(|| argument_error(name, "expected a pair [first, second]"))
}

fn options(map: &Dynamic) -> Result<RevoluteOptions, Box<EvalAltResult>> {
    let o = Options::parse("revolute", map, REVOLUTE_OPTIONS)?;
    let zero = o.with("zero", |a, v| {
        let p = pair(a, v)?;
        Ok((entity(a, &p[0])?, entity(a, &p[1])?))
    })?;
    let limits = o.with("limits", |a, v| {
        let p = pair(a, v)?;
        Ok([args::float(a, &p[0])?, args::float(a, &p[1])?])
    })?;
    let follow = o.with("follow", |a, v| {
        v.clone()
            .into_immutable_string()
            .map(|s| s.to_string())
            .map_err(|_| argument_error(a, "expected the leading parameter's name"))
    })?;
    let ratio = o.with("ratio", args::float)?;
    let follow = match (follow, ratio) {
        (Some(leader), Some(ratio)) => Some((leader, ratio)),
        (Some(leader), None) => Some((leader, 1.0)),
        (None, Some(_)) => {
            return Err(argument_error(
                "ratio",
                "a ratio needs follow: \"leading parameter\"",
            ));
        }
        (None, None) => None,
    };
    Ok(RevoluteOptions {
        zero,
        angle: o.with("angle", args::float)?,
        limits,
        offset: o.with("offset", args::float)?.unwrap_or(0.0),
        follow,
        phase: o.with("phase", args::float)?,
    })
}

fn revolute(
    a: &mut ScriptAssembly,
    name: &str,
    x: ScriptMateRef,
    y: ScriptMateRef,
    map: &Dynamic,
) -> Result<Dynamic, Box<EvalAltResult>> {
    let options = options(map)?;
    let report = a
        .borrow_mut()
        .assembly
        .add_revolute(name, x.0, y.0, &options)
        .assembly()?;
    rhai::serde::to_dynamic(&report).map_err(|e| argument_error("revolute", e.to_string()))
}

/// Register `revolute`, `release_joint` and `drive_joint` (PK-808: a
/// released joint turns with the other mates within its limits).
pub fn register(engine: &mut Engine) {
    engine
        .register_fn(
            "release_joint",
            |a: &mut ScriptAssembly, name: &str| -> Result<(), Box<EvalAltResult>> {
                a.borrow_mut().assembly.release_joint(name).assembly()
            },
        )
        .register_fn(
            "drive_joint",
            |a: &mut ScriptAssembly, name: &str| -> Result<(), Box<EvalAltResult>> {
                a.borrow_mut().assembly.drive_joint(name).assembly()
            },
        )
        .register_fn(
            "revolute",
            |a: &mut ScriptAssembly, name: &str, x: ScriptMateRef, y: ScriptMateRef| {
                revolute(a, name, x, y, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "revolute",
            |a: &mut ScriptAssembly, name: &str, x: ScriptMateRef, y: ScriptMateRef, o: Dynamic| {
                revolute(a, name, x, y, &o)
            },
        );
}
