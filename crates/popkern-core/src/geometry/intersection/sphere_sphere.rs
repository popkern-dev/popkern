//! Exact sphere × sphere intersection (PK-1134).
//!
//! Two spheres meet in a circle, touch at a point, coincide, or miss.
//! The circle lies in their radical plane, at distance
//! `a = (d² + r₁² − r₂²) / 2d` from the first centre along the line of
//! centres, with radius `√(r₁² − a²)`. Without this arm the pair went to
//! the marcher, which on two balls whose circle runs through a pole of
//! each returned 196,069 curves.
//!
//! A touch (one sphere passing into the other by less than the tangency
//! band) emits nothing, as
//! every arm does: contacts are typed by contact discovery beside the
//! curves. Coincident spheres fall through to the coincidence detector.

use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::plane_sphere::{fit_sphere_uv_curve, nurbs_circle_3d};
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// The sphere × sphere answer, or `None` when the pair is not two
/// spheres or they coincide (the coincidence detector's case).
pub(crate) fn try_sphere_sphere(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (
        SurfaceGeometry::Sphere {
            center: c1,
            radius: r1,
            axis: a1,
        },
        SurfaceGeometry::Sphere {
            center: c2,
            radius: r2,
            axis: a2,
        },
    ) = (surface1, surface2)
    else {
        return None;
    };
    let lin = tolerance.linear();
    let line = *c2 - *c1;
    let d = line.length();
    if d <= lin {
        // Concentric: coincident (equal radii) or nested apart.
        return ((r1 - r2).abs() > lin).then(Vec::new);
    }
    if d > r1 + r2 + lin || d < (r1 - r2).abs() - lin {
        return Some(Vec::new());
    }
    // A touch: the spheres pass into each other by less than the
    // tangency band, so their surfaces stay within the tolerance of
    // each other over the whole lens between them, however wide its
    // circle (balls overlapping by half a tolerance meet in a circle of
    // radius 1.2e-3; read as a crossing it was cut across a seam it
    // grazed and a three-ball union lost both domes, PK-1130). Outside
    // (`r₁ + r₂ − d`) or inside (`d − |r₁ − r₂|`), as the other arms
    // judge a touch by `Tolerance::tangency_band`.
    let band = tolerance.tangency_band(super::coordinate_scale(&[*c1, *c2]) + r1.max(*r2));
    if r1 + r2 - d < band || d - (r1 - r2).abs() < band {
        return Some(Vec::new());
    }
    let axis = line * (1.0 / d);
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let rho_sq = r1 * r1 - a * a;
    let rho = rho_sq.sqrt();
    let centre: Point3 = *c1 + axis * a;
    let (e1, e2) = frame(axis);
    let curve_3d = nurbs_circle_3d(centre, e1, e2, rho);
    let curve_uv1 = fit_sphere_uv_curve(&curve_3d, (*c1, *r1, *a1), tolerance);
    let curve_uv2 = fit_sphere_uv_curve(&curve_3d, (*c2, *r2, *a2), tolerance);
    Some(vec![SurfaceSurfaceIntersection::Curve(CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed: true,
    })])
}

/// An orthonormal pair perpendicular to unit `n`.
fn frame(n: Vec3) -> (Vec3, Vec3) {
    let pick = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let e1 = (pick - n * pick.dot(&n))
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    (e1, n.cross(&e1))
}

#[cfg(test)]
#[path = "sphere_sphere_tests.rs"]
mod tests;
