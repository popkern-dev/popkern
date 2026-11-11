//! An operand thinner than the linear tolerance (PK-1131).
//!
//! The tolerance is the resolution: two faces within it of each other
//! are one place. A solid every one of whose faces lies within the
//! tolerance of another of its faces has no inside at that resolution
//! (a slab a millionth thick: its walls coincide, its sides are strips
//! no wider than the tolerance). Cut from a body it removes nothing,
//! added it adds nothing, and what it shares with a body is empty; the
//! boolean resolves it so, deterministically, instead of building an
//! arrangement of faces the tolerance cannot tell apart (a cube cut by
//! a saw 1e-6 wide declined with a broken loop).
//!
//! The test: every probe point inside each face (its interior point and
//! a grid of stations inside its trim, else the middle of its boundary
//! samples) lies within the tolerance of some other face of the solid, measured to that face's trim (the foot on
//! its surface when on its trim, else the nearest of its boundary
//! samples, which can only overstate the distance). One face farther
//! than the tolerance from every other face gives the solid an inside.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::BoolOperation;

/// What the boolean absorbed under the linear tolerance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SubToleranceKind {
    /// An operand thinner than the tolerance everywhere: it removed or
    /// added nothing, and what it shared with the other was empty.
    OperandVanished {
        /// Which operand (0 the target, 1 the tool).
        operand: usize,
    },
    /// The result was thinner than the tolerance everywhere (a lens a
    /// tolerance deep): it is empty.
    ResultVanished,
    /// Facing faces of the operands lay a gap at or under the tolerance
    /// apart: the union fused them across it.
    SlotClosed,
}

/// A feature at or under the linear tolerance the boolean resolved as
/// a coincidence instead of building it (PK-1131): reported, so the
/// caller knows to tighten the tolerance to keep it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SubToleranceFeature {
    /// What was absorbed.
    pub kind: SubToleranceKind,
    /// A point of it.
    pub location: Point3,
    /// Its measured width (at most the tolerance).
    pub measured_width: f64,
    /// The linear tolerance it was resolved at.
    pub tolerance: f64,
}

impl std::fmt::Display for SubToleranceFeature {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.kind {
            SubToleranceKind::OperandVanished { operand } => write!(
                f,
                "operand {operand} is {:.1e} thick at {:?}, under the linear tolerance {:.1e}: \
                 it has no inside and contributed nothing (tighten the tolerance to keep it)",
                self.measured_width, self.location, self.tolerance
            ),
            SubToleranceKind::ResultVanished => write!(
                f,
                "the result is {:.1e} thick at {:?}, under the linear tolerance {:.1e}: \
                 it has no inside and is empty (tighten the tolerance to keep it)",
                self.measured_width, self.location, self.tolerance
            ),
            SubToleranceKind::SlotClosed => write!(
                f,
                "the operands face each other {:.1e} apart at {:?}, under the linear \
                 tolerance {:.1e}: the union closed the gap (tighten the tolerance to keep it)",
                self.measured_width, self.location, self.tolerance
            ),
        }
    }
}

static EVENTS: std::sync::Mutex<Vec<SubToleranceFeature>> = std::sync::Mutex::new(Vec::new());

/// Record a resolved sub-tolerance feature (process registry, drained
/// by [`take_sub_tolerance_events`], and a shadow line).
pub(crate) fn record(feature: SubToleranceFeature) {
    crate::complex_shadow::append_log_line(&format!("SubTolerance {feature}"));
    if let Ok(mut events) = EVENTS.lock() {
        events.push(feature);
    }
}

/// Drain every sub-tolerance feature recorded since the last call.
pub fn take_sub_tolerance_events() -> Vec<SubToleranceFeature> {
    EVENTS
        .lock()
        .map(|mut events| std::mem::take(&mut *events))
        .unwrap_or_default()
}

/// Boundary samples per edge for a face's trimmed distance.
const EDGE_SAMPLES: usize = 32;

/// Where an operand thinner than the tolerance was found, and how thick
/// it is there at most.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ThinOperand {
    /// A point of the operand.
    pub location: Point3,
    /// The largest face-to-face distance found (at most the tolerance).
    pub width: f64,
}

/// Whether `brep` (a solid) is thinner than the linear tolerance
/// everywhere, by the face test above.
pub(crate) fn thin_operand(brep: &BRep, tolerance: &Tolerance) -> Option<ThinOperand> {
    // At the tolerance, allowing its rounding (walls placed 1e-6 apart
    // at x = 5 measure 1.00000000014e-6, faces 1e-6 apart on a turned
    // block 1.0000000018e-6): `Tolerance::within_linear`.
    let within = |d: f64| tolerance.within_linear(d);
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    if faces.len() < 2 {
        return None;
    }
    let boxes: std::collections::HashMap<FaceHandle, _> =
        crate::face_bounding_boxes(brep).into_iter().collect();
    let mut widest: f64 = 0.0;
    let mut location = None;
    for &f in &faces {
        let mut probes = probes::FaceProbes::new(brep, f);
        if probes.is_empty() {
            return None;
        }
        // A probe with no other face within the tolerance gives the
        // solid an inside; boxes farther than it from the probe cannot
        // hold one.
        let clear = probes.find_map(
            |_| true,
            |p| {
                let near = faces
                    .iter()
                    .filter(|&&g| g != f && boxes.get(&g).is_some_and(|b| within(box_gap(b, &p))))
                    .map(|&g| trimmed_distance(brep, g, &p, tolerance))
                    .find(|&d| within(d));
                match near {
                    Some(d) => {
                        widest = widest.max(d);
                        location.get_or_insert(p);
                        None
                    }
                    None => Some(()),
                }
            },
        );
        if clear.is_some() {
            return None;
        }
    }
    Some(ThinOperand {
        location: location?,
        width: widest,
    })
}

/// The distance from `p` to box `b` (zero inside it).
fn box_gap(b: &popkern_core::geometry::intersection::surface_subdivision::Aabb, p: &Point3) -> f64 {
    let gap = |lo: f64, hi: f64, v: f64| (lo - v).max(v - hi).max(0.0);
    let (dx, dy, dz) = (
        gap(b.min.x, b.max.x, p.x),
        gap(b.min.y, b.max.y, p.y),
        gap(b.min.z, b.max.z, p.z),
    );
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// The distance from `p` to face `g` within its trim (an upper bound
/// off the trim: the nearest boundary sample).
fn trimmed_distance(brep: &BRep, g: FaceHandle, p: &Point3, tolerance: &Tolerance) -> f64 {
    let Some(surface) = brep.faces.get(g).and_then(|f| f.surface.as_ref()) else {
        return f64::INFINITY;
    };
    if let Some((q, uv)) = crate::postprocess::shell_contact::foot(&surface.geometry, p)
        && crate::postprocess::shell_contact::on_trim(
            brep,
            g,
            &surface.geometry,
            uv,
            tolerance.linear_reach(),
        )
    {
        return p.distance_to(&q);
    }
    brep.face_outer_loop_sample_positions(g, EDGE_SAMPLES)
        .iter()
        .map(|q| p.distance_to(q))
        .fold(f64::INFINITY, f64::min)
}

/// The result of `op` when operand `thin` (0 = a, 1 = b) has no inside:
/// the other operand for a union or a cut tool, empty otherwise.
pub(crate) fn resolve(op: BoolOperation, a: &BRep, b: &BRep, thin: usize) -> BRep {
    let empty = || BRep::new(a.tolerance);
    match (op, thin) {
        (BoolOperation::Union, 0) => b.clone(),
        (BoolOperation::Union, _) | (BoolOperation::Subtraction, 1) => a.clone(),
        _ => empty(),
    }
}

#[cfg(test)]
#[path = "sub_tolerance_tests.rs"]
mod tests;

#[path = "sub_tolerance_probes.rs"]
mod probes;

#[path = "sub_tolerance_slot.rs"]
mod slot;
pub(crate) use slot::closed_slot;
