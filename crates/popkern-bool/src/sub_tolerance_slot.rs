//! A slot at or under the linear tolerance closed by a union (PK-1131).
//!
//! Two operands whose facing faces lie a gap at or under the tolerance
//! apart are one body at its resolution: the union fuses them across
//! the gap (the plane × plane lines a tolerance off each other's edges
//! lie on those edges), adding up to the faces' area times the gap. A
//! gap measured above rounding was set by the caller, and its closing
//! is reported, as an absorbed operand is.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::coincidence::classify::classify_face_orientation;
use crate::coincidence::types::CoincidenceOrientation;

/// How far past the facing face the side test steps, in tolerances: a
/// slot's far side is the operand's material for at least this depth
/// unless the operand is that thin, when the closing goes unreported.
const SIDE_STEP: f64 = 100.0;

/// Where a face of `a` faces a face of `b` across a gap above rounding
/// and at or under the tolerance, and the gap measured there.
pub(crate) fn closed_slot(a: &BRep, b: &BRep, tolerance: &Tolerance) -> Option<(Point3, f64)> {
    let reach = tolerance.linear_reach();
    let rounding = reach - tolerance.linear();
    let (boxes_a, boxes_b) = (
        crate::intersect::compute_face_aabbs(a),
        crate::intersect::compute_face_aabbs(b),
    );
    for (fa, box_a) in &boxes_a {
        let mut probes: Option<Vec<Point3>> = None;
        for (fb, box_b) in &boxes_b {
            if !box_a.overlaps(box_b, reach) {
                continue;
            }
            let Some(surface) = b.faces.get(*fb).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let probes =
                probes.get_or_insert_with(|| super::probe_points(a, *fa).unwrap_or_default());
            // A probe farther than the tolerance from the other face's
            // box cannot lie within it of the face: no projection.
            let gap = probes
                .iter()
                .filter(|p| super::box_gap(box_b, p) <= reach)
                .find_map(|&p| {
                    let (q, uv) = crate::postprocess::shell_contact::foot(&surface.geometry, &p)?;
                    let d = p.distance_to(&q);
                    let on_face = crate::postprocess::shell_contact::on_trim(
                        b,
                        *fb,
                        &surface.geometry,
                        uv,
                        reach,
                    );
                    (d > rounding && d <= reach && on_face).then_some((p, q, d))
                });
            let Some((p, q, d)) = gap else {
                continue;
            };
            let facing = classify_face_orientation(a, *fa, &boxes_a, b, *fb, &boxes_b, tolerance);
            // Across a slot, past `b`'s face from `a`'s lies `b`'s
            // material; where `b` overlaps `a` by as little, it does
            // not (the step leaves `b` back through `a`).
            let past = q + (q - p) * (SIDE_STEP * tolerance.linear() / d);
            if facing == CoincidenceOrientation::Opposite
                && crate::classify::point_in_solid_final(&past, b, &boxes_b, tolerance)
                    .is_ok_and(|c| c == crate::error::PointClassification::Inside)
            {
                return Some((p + (q - p) * 0.5, d));
            }
        }
    }
    None
}
