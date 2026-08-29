//! §5.3 item #7 — tangent-corner planar cap with G1-tangent-match.
//!
//! At a **tangent corner** (dihedral ≈ π at one or more incident
//! edges — typically where a fillet's trim-edge meets its support
//! face), the rolling-ball contact solver degenerates: there's no
//! well-defined "outside" direction from which a sphere can sit
//! on the fillet strips, and polar-NURBS fits produce curvature
//! artefacts. The usual commercial approach is a **planar cap** whose
//! normal is chosen to align with the fillet surfaces at the
//! boundary arcs — delivering G1-tangent-match at every fillet-
//! end curve without pretending to be curved.
//!
//! # Algorithm
//!
//! 1. For each `FilletEnd` boundary segment, evaluate the incident
//!    fillet surface's unit normal at the arc midpoint.
//! 2. Average those normals as the cap plane's normal.
//! 3. Place the plane origin at the centroid of the boundary arc
//!    endpoints.
//!
//! The resulting plane is tangent to each fillet surface at the
//! midpoint of the corresponding arc (within tangent-approximation
//! error). For a pure "cylinder tangent to plane along a circle"
//! configuration this match is exact; on general tangent
//! configurations the residual angle is bounded by the dihedral
//! deviation from π (which is small by definition of the Tangent
//! classification).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::surface::BlendSurfaceResult;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod};

/// Try to build a planar-cap blend surface for a tangent-classified
/// corner. Returns `None` if the boundary has no fillet-end arcs or
/// the averaged normal is degenerate (zero-length sum, e.g. the
/// incident fillet normals cancel out).
///
/// The returned `BlendSurfaceResult` reports
/// [`CornerMethod::TangentPlane`] and [`ContinuityLevel::G1`] —
/// the caller is expected to emit a
/// [`crate::error::BlendWarning::TangentCornerCap`] alongside it.
pub(super) fn try_tangent_plane_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    points: &[Point3],
    boundary_curves: &[NurbsCurve],
) -> Option<BlendSurfaceResult> {
    if points.is_empty() {
        return None;
    }

    // Average the fillet-surface normal at each arc midpoint. The
    // resulting direction is the "best G1-match" plane normal.
    let mut normal_sum = Vec3::new(0.0, 0.0, 0.0);
    let mut normal_count = 0;
    for seg in &boundary.segments {
        let BoundarySegment::FilletEnd {
            fillet_index,
            curve,
        } = seg
        else {
            continue;
        };
        let Some(fillet) = fillet_surfaces.get(*fillet_index) else {
            continue;
        };
        let (t0, t1) = curve.parameter_range();
        let mid = curve.point_at(0.5 * (t0 + t1));
        let Some((u, v)) = project_point_to_uv(&fillet.surface, &mid) else {
            continue;
        };
        let su = fillet.surface.partial_u(u, v);
        let sv = fillet.surface.partial_v(u, v);
        let n = su.cross(&sv);
        // Skip degenerate partials (zero-length cross product).
        let n_mag = n.length();
        if n_mag < 1e-15 {
            continue;
        }
        let n_unit = Vec3::new(n.x / n_mag, n.y / n_mag, n.z / n_mag);
        normal_sum = Vec3::new(
            normal_sum.x + n_unit.x,
            normal_sum.y + n_unit.y,
            normal_sum.z + n_unit.z,
        );
        normal_count += 1;
    }
    if normal_count == 0 {
        return None;
    }
    let normal = normal_sum.normalized()?;

    // Plane origin = centroid of the boundary points.
    let inv_n = 1.0 / (points.len() as f64);
    let origin = Point3::new(
        points.iter().map(|p| p.x).sum::<f64>() * inv_n,
        points.iter().map(|p| p.y).sum::<f64>() * inv_n,
        points.iter().map(|p| p.z).sum::<f64>() * inv_n,
    );

    // Orthonormal u-axis perpendicular to normal.
    let u_axis = pick_perpendicular(&normal);
    let v_axis = normal.cross(&u_axis).normalized()?;

    Some(BlendSurfaceResult {
        surface: SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        },
        boundary_curves: boundary_curves.to_vec(),
        method: CornerMethod::TangentPlane,
        continuity: ContinuityLevel::G1,
    })
}

/// Pick a unit vector perpendicular to `n`, biased toward the
/// world axis least aligned with `n` so the result is
/// numerically stable.
fn pick_perpendicular(n: &Vec3) -> Vec3 {
    let (ax, ay, az) = (n.x.abs(), n.y.abs(), n.z.abs());
    let seed = if ax <= ay && ax <= az {
        Vec3::new(1.0, 0.0, 0.0)
    } else if ay <= az {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let dot = seed.dot(n);
    let raw = Vec3::new(seed.x - dot * n.x, seed.y - dot * n.y, seed.z - dot * n.z);
    raw.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0))
}

#[cfg(test)]
#[path = "tangent_corner_tests.rs"]
mod tests;
