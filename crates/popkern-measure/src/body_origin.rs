//! Where to take a body's divergence terms about (PK-1135).

use popkern_core::geometry::Point3;
use popkern_core::topology::BRep;

/// The centre of `brep`'s vertex box when the origin lies farther from
/// it than the box's diagonal (`None` otherwise, and for a body whose
/// vertices all coincide: nothing to move it by). The moved body's box
/// centre is at the origin, so this answers `None` for it.
pub(crate) fn far_from_origin(brep: &BRep) -> Option<Point3> {
    let mut points = brep.vertices.values().map(|v| v.position);
    let first = points.next()?;
    let (lo, hi) = points.fold((first, first), |(lo, hi), p| {
        (
            Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z)),
            Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z)),
        )
    });
    let centre = Point3::new(
        0.5 * (lo.x + hi.x),
        0.5 * (lo.y + hi.y),
        0.5 * (lo.z + hi.z),
    );
    let size = lo.distance_to(&hi);
    let far = centre.distance_to(&Point3::new(0.0, 0.0, 0.0));
    (size > 0.0 && far > size).then_some(centre)
}
