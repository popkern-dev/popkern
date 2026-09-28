//! Where a pixel's located surface point sits relative to its facet
//! and to the surface's parameter domain (PK-578): closed-direction
//! seams that are not domain boundaries, and hits outside the visible
//! facet's parameter footprint (silhouette slivers).

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::locate::{active_bounds, unwrap_near};

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

/// How far `uv` lies outside the triangle with parameters `uvs`
/// (already unwrapped onto one sheet), in barycentric units: 0 inside,
/// otherwise the most negative barycentric coordinate's magnitude.
/// A degenerate parameter triangle reports 0.
pub(super) fn footprint_excess(
    surface: Option<&SurfaceGeometry>,
    uvs: &[(f64, f64); 3],
    uv: (f64, f64),
) -> f64 {
    let (pu, pv) = surface.map_or((None, None), |s| s.parameter_periodicity());
    let (x, y) = (
        unwrap_near(uv.0, uvs[0].0, pu),
        unwrap_near(uv.1, uvs[0].1, pv),
    );
    let (a, b, c) = (uvs[0], uvs[1], uvs[2]);
    let det = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
    if det.abs() < 1e-300 {
        return 0.0;
    }
    let l1 = ((x - a.0) * (c.1 - a.1) - (c.0 - a.0) * (y - a.1)) / det;
    let l2 = ((b.0 - a.0) * (y - a.1) - (x - a.0) * (b.1 - a.1)) / det;
    let l0 = 1.0 - l1 - l2;
    (-l0.min(l1).min(l2)).max(0.0)
}
