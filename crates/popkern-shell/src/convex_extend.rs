//! Faces deleted from a convex body, its other faces extended until they
//! meet (PK-1220).
//!
//! Deleting faces and healing the gap runs the neighbouring faces on
//! until they close it. For a convex body whose deleted faces are its
//! convex blends (rounds and corner balls) that healed body is exactly
//! the intersection of the kept faces' supporting regions: each region
//! holds the body (measured as for the shell's cavity, at its vertices,
//! along its edges and inside every face's trim), and the blends lie
//! inside every kept region, so taking them away leaves the regions'
//! intersection. A box rounded on every edge with every round and
//! corner ball deleted is the box again; the corner-by-corner regions
//! declined the cellular assembly where three edge regions and a corner
//! ball met.
//!
//! The intersection must close near the body: kept regions that leave a
//! direction open (a face deleted with no neighbour across it) bound
//! nothing, and the extension declines, as does a body not convex in
//! this sense or a boolean that declines.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::convex_cavity::{body_samples, consumed_convex_blend, region_of};

/// `body` with `deleted` faces taken away and its other faces extended
/// until they meet, when the body is convex and every deleted face a
/// convex blend; `None` otherwise.
pub fn extend_convex_faces(
    body: &BRep,
    deleted: &[FaceHandle],
    tolerance: Tolerance,
) -> Option<BRep> {
    use popkern_bool::error::BoolOperation;
    let lin = tolerance.linear();
    if deleted.is_empty() {
        return None;
    }
    let samples = body_samples(body);
    let (lo, hi) = bounds(&samples)?;
    let size = lo.distance_to(&hi);
    if !(size.is_finite() && size > lin) {
        return None;
    }
    let centre = Point3::new(
        0.5 * (lo.x + hi.x),
        0.5 * (lo.y + hi.y),
        0.5 * (lo.z + hi.z),
    );
    let reach = 2.0 * size;
    let band = 10.0 * lin;
    let mut healed: Option<BRep> = None;
    for (face, data) in &body.faces {
        let surface = &data.surface.as_ref()?.geometry;
        if deleted.contains(&face) {
            // A deleted face must be a convex blend: its tangent planes
            // support the body with its axis or centre inside.
            if !consumed_convex_blend(body, face, surface, &samples, f64::INFINITY, band) {
                return None;
            }
            continue;
        }
        let solid = region_of(surface, &samples, band)?.solid(centre, reach, tolerance)?;
        healed = Some(match healed {
            None => solid,
            Some(h) => {
                popkern_bool::boolean(BoolOperation::Intersection, h, solid, tolerance, None)
                    .ok()?
            }
        });
    }
    let healed = healed?;
    // Closed near the body: an open direction runs out to the regions'
    // own reach.
    let samples_out: Vec<Point3> = healed.vertices.values().map(|v| v.position).collect();
    let (hlo, hhi) = bounds(&samples_out)?;
    if hlo.distance_to(&hhi) > 1.5 * size {
        return None;
    }
    healed
        .validate_collect_structural()
        .is_empty()
        .then_some(healed)
}

/// The axis-aligned bounds of `points`.
fn bounds(points: &[Point3]) -> Option<(Point3, Point3)> {
    let first = *points.first()?;
    Some(points.iter().fold((first, first), |(lo, hi), p| {
        (
            Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z)),
            Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z)),
        )
    }))
}
