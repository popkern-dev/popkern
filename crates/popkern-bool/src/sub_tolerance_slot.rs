//! A slot at or under the linear tolerance closed by a union (PK-1131).
//!
//! Two operands whose facing faces lie a gap at or under the tolerance
//! apart are one body at its resolution: the union fuses them across
//! the gap (the plane × plane lines a tolerance off each other's edges
//! lie on those edges), adding up to the faces' area times the gap. A
//! gap measured above rounding was set by the caller, and its closing
//! is reported, as an absorbed operand is.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::coincidence::classify::classify_face_orientation;
use crate::coincidence::types::CoincidenceOrientation;

/// How far past the facing face the side test steps, in tolerances: a
/// slot's far side is the operand's material for at least this depth
/// unless the operand is that thin, when the closing goes unreported.
const SIDE_STEP: f64 = 100.0;

/// A slot between two operands: where it lies, how wide it is, and a
/// point of each operand's material past its facing face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Slot {
    /// The middle of the gap at the measured probe.
    pub location: Point3,
    /// The gap measured there.
    pub width: f64,
    /// [`SIDE_STEP`] tolerances past `a`'s face, away from `b`.
    pub into_a: Point3,
    /// [`SIDE_STEP`] tolerances past `b`'s face, away from `a`: inside
    /// `b`.
    pub into_b: Point3,
}

/// Where a face of `a` faces a face of `b` across a gap above rounding
/// and at or under the tolerance, and the gap measured there.
pub(crate) fn closed_slot(a: &BRep, b: &BRep, tolerance: &Tolerance) -> Option<Slot> {
    closed_slot_in(
        (a, &crate::intersect::compute_face_aabbs(a)),
        (b, &crate::intersect::compute_face_aabbs(b)),
        tolerance,
    )
}

/// [`closed_slot`] with each operand's face boxes given (PK-1161).
pub(crate) fn closed_slot_in(
    (a, boxes_a): (&BRep, &[(FaceHandle, Aabb)]),
    (b, boxes_b): (&BRep, &[(FaceHandle, Aabb)]),
    tolerance: &Tolerance,
) -> Option<Slot> {
    let reach = tolerance.linear_reach();
    let rounding = reach - tolerance.linear();
    for (fa, box_a) in boxes_a {
        let mut probes = super::probes::FaceProbes::new(a, *fa);
        for (fb, box_b) in boxes_b {
            if !box_a.overlaps(box_b, reach) {
                continue;
            }
            let Some(surface) = b.faces.get(*fb).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let geometry = &surface.geometry;
            // A slot's probe lies within the tolerance and above
            // rounding of the other face's surface: a probe farther
            // than the tolerance from its box, or off that band, is
            // rejected before its own trim is tested.
            let foot = |p: &Point3| {
                crate::postprocess::shell_contact::foot(geometry, p).map(|(q, uv)| {
                    let d = p.distance_to(&q);
                    (q, uv, d)
                })
            };
            let in_band = |p: &Point3| {
                super::box_gap(box_b, p) <= reach
                    && foot(p).is_some_and(|(_, _, d)| d > rounding && d <= reach)
            };
            let gap = probes.find_map(in_band, |p| {
                let (q, uv, d) = foot(&p)?;
                crate::postprocess::shell_contact::on_trim(b, *fb, geometry, uv, reach)
                    .then_some((p, q, d))
            });
            let Some((p, q, d)) = gap else {
                continue;
            };
            let facing = classify_face_orientation(a, *fa, boxes_a, b, *fb, boxes_b, tolerance);
            // Across a slot, past `b`'s face from `a`'s lies `b`'s
            // material; where `b` overlaps `a` by as little, it does
            // not (the step leaves `b` back through `a`).
            let step = SIDE_STEP * tolerance.linear() / d;
            let past = q + (q - p) * step;
            if facing == CoincidenceOrientation::Opposite
                && crate::classify::point_in_solid_final(&past, b, boxes_b, tolerance)
                    .is_ok_and(|c| c == crate::error::PointClassification::Inside)
            {
                return Some(Slot {
                    location: p + (q - p) * 0.5,
                    width: d,
                    into_a: p + (p - q) * step,
                    into_b: past,
                });
            }
        }
    }
    None
}

/// The slots a selection of `bodies` closes: for each pair of solids
/// facing across a gap at or under the tolerance, the slot whose two
/// sides `keep` both keeps (a union-like expression fuses across it; a
/// selection keeping one side leaves a face where it was). Each side's
/// membership is read by classifying its point against every solid.
pub(crate) fn closed_slots_among(
    bodies: &[&BRep],
    keep: impl Fn(&[usize]) -> bool,
    tolerance: &Tolerance,
) -> Vec<Slot> {
    use popkern_core::topology::BodyKind;
    let solids: Vec<usize> = (0..bodies.len())
        .filter(|&i| bodies[i].body_kind() == BodyKind::Solid)
        .collect();
    if solids.len() < 2 {
        return Vec::new();
    }
    let boxes: Vec<_> = bodies
        .iter()
        .map(|b| crate::intersect::compute_face_aabbs(b))
        .collect();
    let members = |p: &Point3| -> Vec<usize> {
        solids
            .iter()
            .copied()
            .filter(|&k| {
                crate::classify::point_in_solid_final(p, bodies[k], &boxes[k], tolerance)
                    .is_ok_and(|c| c == crate::error::PointClassification::Inside)
            })
            .collect()
    };
    let mut out = Vec::new();
    for (n, &i) in solids.iter().enumerate() {
        for &j in &solids[n + 1..] {
            let Some(slot) = closed_slot(bodies[i], bodies[j], tolerance) else {
                continue;
            };
            if keep(&members(&slot.into_a)) && keep(&members(&slot.into_b)) {
                out.push(slot);
            }
        }
    }
    out
}
