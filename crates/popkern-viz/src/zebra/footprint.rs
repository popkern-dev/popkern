//! Where a pixel's located surface point sits relative to its facet
//! and to the surface's parameter domain (PK-578): closed-direction
//! seams that are not domain boundaries, and hits outside the visible
//! facet's 3D footprint (silhouette slivers).

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::dot;
use super::locate::active_bounds;

/// Whether the bound `uv` sits on is a SEAM of a geometrically closed
/// direction rather than a domain edge (PK-578): a clamped NURBS ring
/// reports no periodicity, so its `u = u₁` edge looks like a boundary,
/// but the surface continues across it onto `u₀`. A grazing silhouette
/// ray on coil_spring's tube defeated the ray solves there, and the
/// boundary solve returned a seam point 0.075 away in `u` with a
/// normal 33° off, while the mesh under the pixel matched the true
/// surface. On a seam, the pixel's surface point is not on the bound;
/// the closest-point fallback is.
pub(super) fn on_closed_seam(surface: &SurfaceGeometry, uv: (f64, f64), reach: f64) -> bool {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let [bu, bv] = active_bounds(surface, uv);
    let tol = 1e-7 * reach.max(1e-12);
    let here = surface.point_at(uv.0, uv.1);
    let closed_u = bu.is_some_and(|b| {
        let other = if b == u0 { u1 } else { u0 };
        surface.point_at(other, uv.1).distance_to(&here) <= tol
    });
    let closed_v = bv.is_some_and(|b| {
        let other = if b == v0 { v1 } else { v0 };
        surface.point_at(uv.0, other).distance_to(&here) <= tol
    });
    closed_u || closed_v
}

/// How far `point` lies outside the triangle `tri`, in barycentric
/// units of the triangle's own plane: 0 over the facet, otherwise the
/// most negative barycentric coordinate's magnitude of `point`'s
/// projection. A degenerate triangle reports 0.
///
/// Measured in 3D, not in the chart: a chord along a curving iso-line
/// lies ON the surface but at another parameter (a coil_spring seam
/// chord sits 0.016 in `u` off its vertices' iso-line, 2.7 facet
/// widths), so a parameter-footprint test classified genuine,
/// directly-underneath hits as silhouette slivers and hid a 5.2°
/// display error from the view statistics (PK-587). A silhouette
/// sliver is a hit displaced along the grazing ray, which is outside
/// the facet in 3D as well.
pub(super) fn footprint_excess(tri: &[[f64; 3]; 3], point: [f64; 3]) -> f64 {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let (e1, e2, r) = (sub(tri[1], tri[0]), sub(tri[2], tri[0]), sub(point, tri[0]));
    let (a, b, c) = (dot(e1, e1), dot(e1, e2), dot(e2, e2));
    let det = a * c - b * b;
    if det.is_nan() || det <= 1e-24 * a * c {
        return 0.0;
    }
    let (p, q) = (dot(r, e1), dot(r, e2));
    let l1 = (c * p - b * q) / det;
    let l2 = (a * q - b * p) / det;
    let l0 = 1.0 - l1 - l2;
    (-l0.min(l1).min(l2)).max(0.0)
}
