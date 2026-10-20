//! Fastener insertion in films (PK-854): `s.insert(path, #{ distance })`
//! slides a fastener occurrence in along its own axis, and
//! `s.screw_in(path, #{ distance })` turns it one thread lead a turn as
//! it goes. A fastener's axis is its insertion frame — a fastener placed
//! on a hole (`seat`, `frame`) has its local +z down the hole from the
//! seat — and its lead the thread's pitch times starts, both read from
//! the built assembly into the film's document snapshot. The tracks are
//! the plain moves and turns a film could write by hand: the fastener
//! ends seated, so the film's last frame is the real assembly.

use std::collections::BTreeMap;

use popkern_assembly::{Assembly, Placement};
use popkern_history::feature::Feature;
use rhai::{Dynamic, EvalAltResult, NativeCallContext};
use serde::{Deserialize, Serialize};

use super::api::{ScriptStep, target};
use super::model::Channel;
use crate::args::{self, Options};
use crate::error::argument_error;

/// A fastener occurrence's axis, in document units.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FastenerAxis {
    /// The seat: the occurrence's world origin.
    pub origin: [f64; 3],
    /// Unit direction out of the hole, toward the head (the insertion
    /// frame's -z).
    pub out: [f64; 3],
    /// Axial advance per turn.
    pub lead: f64,
    /// Whether the thread is left-handed.
    pub left_handed: bool,
}

/// The axis of every leaf occurrence of `assembly` whose part is one
/// fastener (a part built by exactly one fastener feature).
pub(super) fn fastener_axes(assembly: &Assembly) -> BTreeMap<String, FastenerAxis> {
    let mut out = BTreeMap::new();
    for row in assembly.leaves() {
        let Some(model) = assembly
            .definition(row.definition)
            .ok()
            .and_then(|d| d.model())
        else {
            continue;
        };
        let tree = model.tree();
        let threads: Vec<_> = tree
            .timeline()
            .iter()
            .filter_map(|&id| match &tree.feature(id)?.feature {
                Feature::Fastener { thread, .. } => Some(thread.clone()),
                _ => None,
            })
            .collect();
        let [thread] = threads.as_slice() else {
            continue;
        };
        let Some(pitch) = tree.parameters().get(thread.pitch).map(|p| p.value) else {
            continue;
        };
        // Parameters are SI; the assembly places parts in their own
        // document units.
        let unit = model.length_unit().to_si(1.0);
        let lead = pitch * f64::from(thread.starts.max(1)) / unit;
        if !(lead.is_finite() && lead > 0.0) {
            continue;
        }
        out.insert(
            row.path.to_string(),
            axis_of(&row.world, lead, thread.left_handed),
        );
    }
    out
}

fn axis_of(world: &Placement, lead: f64, left_handed: bool) -> FastenerAxis {
    let (o, z) = (world.origin(), world.z_axis());
    FastenerAxis {
        origin: [o.x, o.y, o.z],
        out: [-z.x, -z.y, -z.z],
        lead,
        left_handed,
    }
}

const INSERT_OPTIONS: &[&str] = &["distance", "at", "duration", "ease"];

/// `s.insert(path, #{ distance })`: slide in from `distance` out along
/// the fastener's axis to its seat.
pub(super) fn insert(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    path: Dynamic,
    options: Dynamic,
) -> Result<(), Box<EvalAltResult>> {
    insertion(ctx, s, path, options, "insert", false)
}

/// `s.screw_in(path, #{ distance })`: [`insert`], turning one lead a
/// turn about the axis (the thread's hand decides the sense).
pub(super) fn screw_in(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    path: Dynamic,
    options: Dynamic,
) -> Result<(), Box<EvalAltResult>> {
    insertion(ctx, s, path, options, "screw_in", true)
}

fn insertion(
    ctx: NativeCallContext,
    s: &mut ScriptStep,
    path: Dynamic,
    options: Dynamic,
    action: &'static str,
    turning: bool,
) -> Result<(), Box<EvalAltResult>> {
    let context = s.context();
    let target = target(&context, "target", &path)?;
    let Some(axis) = context.fasteners.get(&target) else {
        let known: Vec<&str> = context
            .fasteners
            .keys()
            .map(String::as_str)
            .take(24)
            .collect();
        return Err(argument_error(
            "target",
            format!(
                "'{target}' is not a fastener occurrence with an axis{}",
                if known.is_empty() {
                    String::new()
                } else {
                    format!("; the fasteners are {}", known.join(", "))
                }
            ),
        ));
    };
    let options = Options::parse(action, &options, INSERT_OPTIONS)?;
    let distance = options
        .with("distance", args::float)?
        .ok_or_else(|| argument_error("options", format!("{action} needs distance: a length")))?;
    if !(distance.is_finite() && distance >= 0.0) {
        return Err(argument_error("distance", "a finite length >= 0"));
    }
    let timing = s.timing(&options)?;
    let si = context.unit_scale;
    let from = axis.out.map(|c| c * distance * si);
    s.push(
        &ctx,
        timing,
        Channel::Move {
            target: target.clone(),
            from,
            to: [0.0; 3],
        },
    );
    if turning {
        // Backed out by `distance`, a right-hand thread sits turned
        // positively about its out-axis by one turn a lead.
        let turns = distance / axis.lead * 360.0;
        s.push(
            &ctx,
            timing,
            Channel::Turn {
                target,
                origin: axis.origin.map(|c| c * si),
                direction: axis.out,
                from: if axis.left_handed { -turns } else { turns },
                to: 0.0,
            },
        );
    }
    Ok(())
}
