//! Geometry primitives for the N-quad ring decomposition: the
//! central-vertex normal, its saddle fallback, and the straight
//! spoke builder. Split from `gregory_three_quad.rs` for the
//! 500-line structural cap.

use super::boundary::{BlendBoundary, BoundarySegment};
use super::gregory_parallelogram::sub3;
use crate::fillet_surface::FilletSurfaceResult;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

/// Approximate outward normal of the hex plane at the central
/// vertex `C`. Sums cross products of consecutive
/// centroid-to-hex-vertex vectors — for a convex hex in CCW
/// order the sum points along the plane normal. Returns `None`
/// if the 6 cross products sum to a near-zero vector
/// (degenerate hex, e.g., all 6 vertices on a line, or a saddle
/// where the positive + negative curvature contributions cancel).
pub(super) fn central_vertex_normal(hex_vertices: &[Point3], c: Point3) -> Option<Vec3> {
    if hex_vertices.len() != 6 {
        return None;
    }
    let mut nx = 0.0_f64;
    let mut ny = 0.0_f64;
    let mut nz = 0.0_f64;
    for k in 0..6 {
        let a = sub3(hex_vertices[k], c);
        let b = sub3(hex_vertices[(k + 1) % 6], c);
        let cross = a.cross(&b);
        nx += cross.x;
        ny += cross.y;
        nz += cross.z;
    }
    let n = Vec3::new(nx, ny, nz);
    n.normalized()
}

/// §5.3 item #5 — saddle central-vertex normal fallback.
///
/// When [`central_vertex_normal`] returns `None` (saddle-case
/// cancellation of the cross-product sum), pick the normal
/// direction from the fillet-surface normals at each fillet-end
/// arc midpoint. Averaging the absolute-valued normals (to avoid
/// sign cancellation across the ridge) gives a direction
/// perpendicular to the saddle ridge — the correct common
/// tangent-plane normal at C.
pub(super) fn saddle_normal_fallback(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<Vec3> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    let mut n = 0usize;
    let mut seed: Option<Vec3> = None;
    for (i, seg) in boundary.segments.iter().enumerate() {
        let BoundarySegment::FilletEnd { curve, .. } = seg else {
            continue;
        };
        // Pair each FilletEnd arc (even index 2k) with the fillet
        // strip at the same index.
        let strip_idx = i / 2;
        let Some(strip) = fillet_surfaces.get(strip_idx) else {
            continue;
        };
        let (t0, t1) = curve.parameter_range();
        let mid = curve.point_at(0.5 * (t0 + t1));
        let Some((u, v)) = strip.surface.closest_point(&mid) else {
            continue;
        };
        let normal = strip.surface.normal_at(u, v);
        let Some(unit) = normal.normalized() else {
            continue;
        };
        // Flip the contribution if it's anti-aligned with the
        // first sampled normal — this ensures the average is in
        // one consistent half-space across the saddle ridge.
        let aligned = match seed {
            None => {
                seed = Some(unit);
                unit
            }
            Some(ref s) => {
                if unit.dot(s) >= 0.0 {
                    unit
                } else {
                    Vec3::new(-unit.x, -unit.y, -unit.z)
                }
            }
        };
        sum = Vec3::new(sum.x + aligned.x, sum.y + aligned.y, sum.z + aligned.z);
        n += 1;
    }
    if n == 0 {
        return None;
    }
    sum.normalized()
}

/// Build a straight-line spoke as a degree-1 NURBS from `a` to `b`
/// through a midpoint. Returns `None` if the two points coincide.
pub(super) fn build_spoke(a: Point3, b: Point3) -> Option<NurbsCurve> {
    if a.distance_to(&b) < 1e-15 {
        return None;
    }
    let mid = Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5);
    crate::center_curve::fit_nurbs_through_points(&[a, mid, b]).ok()
}
