//! Facts about a complex edge's faces the radial order reads.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

use super::ids::FinId;
use super::model::CellComplex;

/// Whether a fin's face lies on a plane (its primary owner's surface).
/// A straight edge between planes is ordered at one station: each
/// fin's inward direction n × t is the same along it (on a 48-gon
/// prism's union with a box, 0.24 of the complex's 0.93 ms).
pub(super) fn on_plane(surface: Option<&SurfaceOnFace>) -> bool {
    surface.is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Plane { .. }))
}

/// How far `fin`'s face reaches from `p` before another of its edges:
/// the distance from `p` to the nearest edge of the face other than the
/// fin's own (infinite for a face bounded by that edge alone). Faces
/// tied at first order are ordered within this reach: a lens between
/// two faces crossing at a small angle is exactly this wide, and beyond
/// it their order is the other crossing's.
pub(super) fn face_clearance(c: &CellComplex, fin: FinId, p: Point3) -> f64 {
    use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};
    let Some(f) = c.fins.get(fin.index()) else {
        return f64::INFINITY;
    };
    let Some(face) = c.faces.get(f.face.index()) else {
        return f64::INFINITY;
    };
    let own = f.edge;
    face.loops
        .iter()
        .flat_map(|&l| c.loop_fins(l))
        .filter_map(|g| {
            let e = c.fins[g.index()].edge;
            (e != own).then(|| &c.edges[e.index()])
        })
        .map(|edge| match &edge.curve {
            Some(cv) => {
                let t = closest_point_in_window(&cv.geometry, &p, cv.t_start, cv.t_end);
                cv.geometry.point_at(t).distance_to(&p)
            }
            None => {
                let (a, b) = (
                    c.vertices[edge.start.index()].position,
                    c.vertices[edge.end.index()].position,
                );
                let ab = b - a;
                let k = ((p - a).dot(&ab) / ab.dot(&ab).max(1e-300)).clamp(0.0, 1.0);
                (a + ab * k).distance_to(&p)
            }
        })
        .fold(f64::INFINITY, f64::min)
}
