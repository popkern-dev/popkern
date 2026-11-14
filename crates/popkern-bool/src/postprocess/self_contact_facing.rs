//! Whether two faces sharing a vertex face each other there, and
//! whether a touch found beside the vertex is its own (split from
//! `self_contact.rs` under the file-size rule).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;

/// Whether the faces on `a` and `b` face each other at `p`: their
/// outward unit normals there opposed within a hundredth of a radian's
/// cosine (a tangent pinch, the faces' interiors on either side of one
/// tangent plane). The chart normal `∂u × ∂v` is scaled by the chart's
/// speed: a ball of radius 3 read it at 9 away from its poles, and 1e-4
/// from a pole at 3e-4, which never passed for opposed (PK-1187).
pub(super) fn facing_at(
    a: &popkern_core::geometry::surface::SurfaceOnFace,
    b: &popkern_core::geometry::surface::SurfaceOnFace,
    p: Point3,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    let outward = |s: &popkern_core::geometry::surface::SurfaceOnFace| {
        let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(&s.geometry, &p)?;
        let n = s.geometry.normal_at(u, v).normalized()?;
        Some(if s.orientation { n } else { n * -1.0 })
    };
    match (outward(a), outward(b)) {
        (Some(x), Some(y)) => x.dot(&y) < -0.99995,
        _ => false,
    }
}

/// Whether the gap between surfaces `a` and `b` at `at` keeps falling
/// toward the nearest of the `shared` vertices: the touch found there is
/// the vertex's, reached along a band under the tolerance, not a touch
/// of its own. A dimple whose rim touches a plate's edge pulls away from
/// it as y²/8: scaled by 1e-3 the wall stays under the tolerance 9e-5
/// from the pinch, and the search reported a contact on that slope
/// 3.6e-5 from it, where a ball resting beside a shared vertex has a
/// minimum of its own (PK-1187).
pub(super) fn falls_to_shared(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    at: Point3,
    shared: &[Point3],
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let gap = |q: Point3| -> Option<f64> {
        let (u, v) = project_point_to_uv(a, &q)?;
        let on_a = a.point_at(u, v);
        let (s, t) = project_point_to_uv(b, &on_a)?;
        Some(b.point_at(s, t).distance_to(&on_a))
    };
    let Some(s) = shared
        .iter()
        .min_by(|p, q| p.distance_to(&at).total_cmp(&q.distance_to(&at)))
    else {
        return false;
    };
    let toward = at + (*s - at) * 0.25;
    matches!((gap(at), gap(toward)), (Some(here), Some(nearer)) if nearer < here)
}
