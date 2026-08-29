//! Search hints for the PK-255 pole-placement search (see
//! docs/architecture/exact-boundary-pole-placement.md): the blend-ward
//! pull direction and scale, derived from the prepared fillet-end
//! arcs alone. Each arc spans its strip's end-section plane (three
//! arc samples give the plane); the OTHER arcs' midpoints give the
//! blend side; the pull is the negated sum of the confident
//! into-strip normals and the scale is the largest arc chord.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

/// Orientation is trusted only when the reference clears the plane by
/// this fraction of the chord.
const ORIENTATION_CONFIDENCE_FRACTION: f64 = 1.0e-2;

/// `Some((direction, scale))` when at least one plane orientation is
/// confident; `None` disables the pole search for the corner.
pub(super) fn pole_search_hints(fillet_arcs: &[&NurbsCurve]) -> Option<(Vec3, f64)> {
    let mids: Vec<Point3> = fillet_arcs
        .iter()
        .map(|arc| {
            let (t0, t1) = arc.parameter_range();
            arc.point_at(0.5 * (t0 + t1))
        })
        .collect();
    let mut pull = Vec3::new(0.0, 0.0, 0.0);
    let mut scale = 0.0_f64;
    let mut confident = 0usize;
    for (arc_index, arc) in fillet_arcs.iter().enumerate() {
        let (t0, t1) = arc.parameter_range();
        let start = arc.point_at(t0);
        let mid = arc.point_at(0.5 * (t0 + t1));
        let end = arc.point_at(t1);
        let chord = start.distance_to(&end);
        if chord <= 1.0e-12 {
            continue;
        }
        scale = scale.max(chord);
        let Some(normal) = (mid - start).cross(&(end - start)).normalized() else {
            continue;
        };
        let others: Vec<&Point3> = mids
            .iter()
            .enumerate()
            .filter(|(k, _)| *k != arc_index)
            .map(|(_, p)| p)
            .collect();
        if others.is_empty() {
            continue;
        }
        let denom = others.len() as f64;
        let reference = Point3::new(
            others.iter().map(|p| p.x).sum::<f64>() / denom,
            others.iter().map(|p| p.y).sum::<f64>() / denom,
            others.iter().map(|p| p.z).sum::<f64>() / denom,
        );
        let side = (reference - mid).dot(&normal);
        if side.abs() < chord * ORIENTATION_CONFIDENCE_FRACTION {
            continue;
        }
        let into_strip = if side > 0.0 { -normal } else { normal };
        pull = pull - into_strip;
        confident += 1;
    }
    if confident == 0 || scale <= 0.0 {
        return None;
    }
    pull.normalized().map(|direction| (direction, scale))
}
