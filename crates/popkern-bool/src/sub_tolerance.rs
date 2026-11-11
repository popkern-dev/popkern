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
    /// Two faces of the result facing the same way lay a step at or
    /// under the tolerance apart: the boolean made them one face.
    StepFlattened,
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
            SubToleranceKind::StepFlattened => write!(
                f,
                "two faces lay a step {:.1e} apart at {:?}, under the linear tolerance \
                 {:.1e}: the boolean made them one face (tighten the tolerance to keep it)",
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

/// Record each slot a selection closed as
/// [`SubToleranceKind::SlotClosed`].
pub(crate) fn record_closed_slots(slots: &[slot::Slot], tolerance: &Tolerance) {
    for s in slots {
        record(SubToleranceFeature {
            kind: SubToleranceKind::SlotClosed,
            location: s.location,
            measured_width: s.width,
            tolerance: tolerance.linear(),
        });
    }
}

/// The steps a pass of face merges flattened, gathered so that one
/// step is reported once: two merges across one step (the tops of
/// blocks overlapping, the coincident piece between them merged into
/// each side in turn) are one feature, the faces they joined one face
/// (PK-1163).
#[derive(Default)]
pub(crate) struct FlattenedSteps {
    /// Each merge: the two faces it joined and the step between them,
    /// if one above rounding.
    merges: Vec<(
        [popkern_core::topology::FaceHandle; 2],
        Option<SubToleranceFeature>,
    )>,
}

impl FlattenedSteps {
    /// Note the merge about to join the faces of `edge`, measured before
    /// it (`commit` once it has succeeded).
    pub(crate) fn measure(
        brep: &BRep,
        edge: popkern_core::topology::EdgeHandle,
        tolerance: &Tolerance,
    ) -> Option<(
        [popkern_core::topology::FaceHandle; 2],
        Option<SubToleranceFeature>,
    )> {
        let he = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
        let twin = brep.half_edges.get(he.mate()?)?;
        let (fa, fb) = (
            brep.loops.get(he.loop_)?.face,
            brep.loops.get(twin.loop_)?.face,
        );
        let rounding = tolerance.linear_reach() - tolerance.linear();
        let step = crate::postprocess::merge_rechart::faces_offset(brep, fa, fb)
            .filter(|&(d, _)| d > rounding && tolerance.within_linear(d))
            .map(|(d, at)| SubToleranceFeature {
                kind: SubToleranceKind::StepFlattened,
                location: at,
                measured_width: d,
                tolerance: tolerance.linear(),
            });
        Some(([fa, fb], step))
    }

    /// Keep a measured merge that succeeded.
    pub(crate) fn commit(
        &mut self,
        merge: Option<(
            [popkern_core::topology::FaceHandle; 2],
            Option<SubToleranceFeature>,
        )>,
    ) {
        if let Some(m) = merge {
            self.merges.push(m);
        }
    }

    /// Record the widest step of each set of faces the merges joined.
    pub(crate) fn record(self) {
        use popkern_core::topology::FaceHandle;
        let mut parent: std::collections::HashMap<FaceHandle, FaceHandle> = Default::default();
        fn find(
            parent: &mut std::collections::HashMap<
                popkern_core::topology::FaceHandle,
                popkern_core::topology::FaceHandle,
            >,
            f: popkern_core::topology::FaceHandle,
        ) -> popkern_core::topology::FaceHandle {
            let p = *parent.entry(f).or_insert(f);
            if p == f {
                return f;
            }
            let root = find(parent, p);
            parent.insert(f, root);
            root
        }
        for ([a, b], _) in &self.merges {
            let (ra, rb) = (find(&mut parent, *a), find(&mut parent, *b));
            if ra != rb {
                parent.insert(ra, rb);
            }
        }
        let mut widest: std::collections::BTreeMap<FaceHandle, SubToleranceFeature> =
            Default::default();
        for ([a, _], step) in self.merges {
            let Some(step) = step else {
                continue;
            };
            let root = find(&mut parent, a);
            let slot = widest.entry(root).or_insert(step);
            if step.measured_width > slot.measured_width {
                *slot = step;
            }
        }
        for step in widest.into_values() {
            record(step);
        }
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
    thin_operand_in(brep, &[], tolerance)
}

/// [`thin_operand`] with whatever of `brep`'s face boxes are built (the
/// boolean builds each operand's once for its sub-tolerance tests and
/// intersection stage, PK-1161); a face on a surface with no closed
/// form whose box is not given has it built when a probe asks.
pub(crate) fn thin_operand_in(
    brep: &BRep,
    face_boxes: &[(
        FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    tolerance: &Tolerance,
) -> Option<ThinOperand> {
    // At the tolerance, allowing its rounding (walls placed 1e-6 apart
    // at x = 5 measure 1.00000000014e-6, faces 1e-6 apart on a turned
    // block 1.0000000018e-6): `Tolerance::within_linear`.
    let within = |d: f64| tolerance.within_linear(d);
    let mut faces: Vec<FaceHandle> = brep.faces.keys().collect();
    if faces.len() < 2 {
        return None;
    }
    // A body with an inside ends the test at its first face's interior
    // point, whose cost grows with the face's loops (a 48-gon cap's
    // took 0.05 ms, the whole test on a prism): faces of fewest
    // half-edges first (PK-1161). A thin body tests every face, so the
    // order changes only which point reports it.
    let mut sides: std::collections::HashMap<FaceHandle, usize> = std::collections::HashMap::new();
    for he in brep.half_edges.values() {
        if let Some(l) = brep.loops.get(he.loop_) {
            *sides.entry(l.face).or_default() += 1;
        }
    }
    faces.sort_by_key(|f| sides.get(f).copied().unwrap_or(0));
    let near_faces = near::NearFaces::new(brep, face_boxes);
    let reach = tolerance.linear_reach();
    let mut widest: f64 = 0.0;
    let mut location = None;
    for &f in &faces {
        let mut probes = probes::FaceProbes::new(brep, f);
        if probes.is_empty() {
            return None;
        }
        // A probe with no other face within the tolerance gives the
        // solid an inside; a face whose surface or box is farther than
        // it from the probe cannot be within it.
        let clear = probes.find_map(
            |_| true,
            |p| {
                let near = faces
                    .iter()
                    .filter(|&&g| g != f && near_faces.may_reach(g, &p, reach))
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

#[path = "sub_tolerance_near.rs"]
mod near;
use near::box_gap;

#[path = "sub_tolerance_slot.rs"]
mod slot;
pub(crate) use slot::{closed_slot_in, closed_slots_among};
