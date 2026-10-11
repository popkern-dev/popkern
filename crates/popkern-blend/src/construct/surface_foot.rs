//! Exact support normals for the setback rim passes (PK-776).
//!
//! The rim jet, the curvature stations and the band refit each took
//! the support's normal by central differences (`h = 1e-6` on point
//! evaluations). The noise of that estimate, about 1e-10, reached the
//! cap's corner pieces, which are ~0.01 across. Curvature there scales
//! with control-point noise over the squared piece size, so the rim's
//! measured G2 jump moved by 25-46% under a one-ulp change of the
//! input. The surface's own partials repeat to roundoff.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// The unit normal of `surface` at `uv` from its own first partials
/// (`S_u × S_v`); `None` where the chart degenerates (one partial under
/// 1e-3 of the other, or a vanishing cross product).
pub(crate) fn unit_normal(surface: &SurfaceGeometry, (u, v): (f64, f64)) -> Option<Vec3> {
    let (su, sv) = (surface.partial_u(u, v), surface.partial_v(u, v));
    let (lu, lv) = (su.length(), sv.length());
    if lu.min(lv) < 1e-3 * lu.max(lv) {
        return None;
    }
    su.cross(&sv).normalized()
}
