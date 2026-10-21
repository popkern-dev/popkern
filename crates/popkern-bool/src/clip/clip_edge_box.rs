//! Conservative edge boxes for the clip's per-boundary-edge loops
//! (PK-815): exact prunings that skip edges which cannot matter, so the
//! loops cost what the record reaches rather than the face's boundary.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle};

/// A box holding edge `eh`: its two vertices plus its curve's
/// conservative extent (`face_box::conservative_curve_extent`).
pub(crate) fn edge_box(brep: &BRep, eh: EdgeHandle) -> Option<Aabb> {
    let edge = brep.edges.get(eh)?;
    let he = brep.half_edges.get(edge.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    let mut points = vec![
        brep.vertices.get(he.origin)?.position,
        brep.vertices.get(twin.origin)?.position,
    ];
    crate::face_box::conservative_curve_extent(edge.curve.as_ref()?, &mut points);
    Aabb::from_points(&points)
}

/// The distance from `p` to box `b` (zero inside): no point of the box,
/// so none of the geometry it holds, is nearer.
pub(crate) fn box_distance(b: &Aabb, p: &Point3) -> f64 {
    let d = |x: f64, lo: f64, hi: f64| (lo - x).max(0.0).max(x - hi);
    let (dx, dy, dz) = (
        d(p.x, b.min.x, b.max.x),
        d(p.y, b.min.y, b.max.y),
        d(p.z, b.min.z, b.max.z),
    );
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Whether edge `eh` lies wholly on one side of `surface`, a plane, by
/// more than `margin` (every corner of its box): it cannot meet it.
/// `false` for any other surface, or an edge without a box.
pub(crate) fn edge_clear_of_plane(
    brep: &BRep,
    eh: EdgeHandle,
    surface: &SurfaceGeometry,
    margin: f64,
) -> bool {
    let SurfaceGeometry::Plane { origin, normal, .. } = surface else {
        return false;
    };
    let (Some(n), Some(b)) = (normal.normalized(), edge_box(brep, eh)) else {
        return false;
    };
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for x in [b.min.x, b.max.x] {
        for y in [b.min.y, b.max.y] {
            for z in [b.min.z, b.max.z] {
                let s = (Point3::new(x, y, z) - *origin).dot(&n);
                (lo, hi) = (lo.min(s), hi.max(s));
            }
        }
    }
    lo > margin || hi < -margin
}
