//! Exact recompute arms for curved analytic skins (PK-322).
//!
//! Scope: the analytic pairs natural-quadric bodies produce —
//! axis-aligned plane × (cylinder | cone) circles, straight-ruling seam
//! edges, and vertex positions on mixed plane/quadric intersections.
//! Every arm verifies its candidate against the offset surfaces by
//! sampled residuals and returns `None` on any doubt, so unhandled
//! configurations keep their (stale) clone curves and decline honestly
//! at the serve gate. General curved pairs (numeric SSI marching,
//! valence>3 least-squares vertices) are PK-324.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Exact closest point on an analytic surface, `None` for classes
/// without a closed form here or at parametrically ambiguous inputs
/// (points on an axis, at a center, or at a cone apex).
pub(crate) fn project_to_analytic_surface(surface: &SurfaceGeometry, p: &Point3) -> Option<Point3> {
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let n = normal.normalized()?;
            Some(*p - n * (*p - *origin).dot(&n))
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let a = axis.normalized()?;
            let rel = *p - *origin;
            let ax = rel.dot(&a);
            let radial = rel - a * ax;
            let radial_n = radial.normalized()?;
            Some(*origin + a * ax + radial_n * *radius)
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let a = axis.normalized()?;
            let rel = *p - *apex;
            let t = rel.dot(&a);
            let radial = rel - a * t;
            let r = radial.length();
            let radial_n = radial.normalized()?;
            // In (axial, radial) coordinates the ruling is the ray
            // s·(cosα, sinα), s ≥ 0; project onto it.
            let (sin_a, cos_a) = half_angle.psin_cos();
            let s = t * cos_a + r * sin_a;
            if s <= 0.0 {
                // Nearest surface point is the apex — parametrically
                // degenerate; the apex arm owns that case.
                return None;
            }
            Some(*apex + a * (s * cos_a) + radial_n * (s * sin_a))
        }
        SurfaceGeometry::Sphere { center, radius, .. } => {
            let rel_n = (*p - *center).normalized()?;
            Some(*center + rel_n * *radius)
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let a = axis.normalized()?;
            let rel = *p - *center;
            let ring_dir = (rel - a * rel.dot(&a)).normalized()?;
            let ring_center = *center + ring_dir * *major_radius;
            let tube_n = (*p - ring_center).normalized()?;
            Some(ring_center + tube_n * *minor_radius)
        }
        _ => None,
    }
}

/// Exact distance from `p` to the analytic surface.
///
/// Computed directly (not via [`project_to_analytic_surface`]) so that
/// parametrically ambiguous points — the cone apex, points on an axis —
/// still get a well-defined residual: a seam edge's end vertex IS the
/// offset cone's apex, and verification must accept it.
pub(crate) fn analytic_residual(surface: &SurfaceGeometry, p: &Point3) -> Option<f64> {
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let n = normal.normalized()?;
            Some((*p - *origin).dot(&n).abs())
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let a = axis.normalized()?;
            let rel = *p - *origin;
            let radial = rel - a * rel.dot(&a);
            Some((radial.length() - radius).abs())
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let a = axis.normalized()?;
            let rel = *p - *apex;
            let t = rel.dot(&a);
            let r = (rel - a * t).length();
            let (sin_a, cos_a) = half_angle.psin_cos();
            let s = t * cos_a + r * sin_a;
            if s <= 0.0 {
                // Nearest surface point is the apex itself.
                Some((t * t + r * r).sqrt())
            } else {
                Some((r * cos_a - t * sin_a).abs())
            }
        }
        SurfaceGeometry::Sphere { center, radius, .. } => {
            Some(((*p - *center).length() - radius).abs())
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let a = axis.normalized()?;
            let rel = *p - *center;
            let z = rel.dot(&a);
            let rho = (rel - a * z).length();
            let d_ring = ((rho - major_radius).powi(2) + z * z).sqrt();
            Some((d_ring - minor_radius).abs())
        }
        _ => None,
    }
}

/// Solve a vertex position on the intersection of analytic surfaces by
/// sequential exact projection, seeded from the original vertex (branch
/// selection: nearest intersection to the seed).
///
/// Returns `Some` only when the iterate converges onto every surface
/// within half the linear tolerance — a non-converged position must not
/// replace the clone (the serve gate declines the stale clone; a
/// half-moved vertex could pass validation while being wrong).
pub(crate) fn solve_vertex_on_surfaces(
    surfaces: &[&SurfaceGeometry],
    seed: Point3,
    tolerance: &Tolerance,
) -> Option<Point3> {
    if surfaces.is_empty() {
        return None;
    }
    let accept = tolerance.linear() * 0.5;
    let mut pos = seed;
    for _ in 0..200 {
        let mut moved: f64 = 0.0;
        for surf in surfaces {
            let q = project_to_surface(surf, &pos)?;
            moved = moved.max(pos.distance_to(&q));
            pos = q;
        }
        if moved < 1e-14 {
            break;
        }
    }
    let converged = surfaces
        .iter()
        .all(|s| surface_distance(s, &pos).is_some_and(|r| r <= accept));
    converged.then_some(pos)
}

/// Closest point on any surface: exact analytic arms first, iterative
/// UV projection for the rest (fitted NURBS offsets).
pub(crate) fn project_to_surface(surface: &SurfaceGeometry, p: &Point3) -> Option<Point3> {
    use popkern_core::geometry::surface::SurfaceEval;
    if let Some(q) = project_to_analytic_surface(surface, p) {
        return Some(q);
    }
    if matches!(
        surface,
        SurfaceGeometry::Plane { .. }
            | SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Cone { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Torus { .. }
    ) {
        // Analytic arm declined for a REASON (parametric ambiguity —
        // axis point, cone apex); the iterative fallback would pick an
        // arbitrary branch. Keep the decline.
        return None;
    }
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, p)?;
    Some(surface.point_at(u, v))
}

/// Distance to any surface: exact analytic residuals first, iterative
/// projection distance for the rest.
pub(crate) fn surface_distance(surface: &SurfaceGeometry, p: &Point3) -> Option<f64> {
    if let Some(r) = analytic_residual(surface, p) {
        return Some(r);
    }
    Some(project_to_surface(surface, p)?.distance_to(p))
}

/// Quadric axis data for the coaxial-circle arm.
struct AxisForm {
    origin: Point3,
    axis: Vec3,
    kind: AxisKind,
}

enum AxisKind {
    Cylinder { radius: f64 },
    Cone { half_angle: f64 },
}

fn axis_form(surface: &SurfaceGeometry) -> Option<AxisForm> {
    match surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some(AxisForm {
            origin: *origin,
            axis: axis.normalized()?,
            kind: AxisKind::Cylinder { radius: *radius },
        }),
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => Some(AxisForm {
            origin: *apex,
            axis: axis.normalized()?,
            kind: AxisKind::Cone {
                half_angle: *half_angle,
            },
        }),
        _ => None,
    }
}

/// Apply the coaxial affine map `p ↦ c1 + (p − c0)·s` to a curve
/// geometry, preserving parametrization and direction exactly.
fn map_coaxial_circle_geometry(
    geometry: &CurveGeometry,
    c0: Point3,
    c1: Point3,
    s: f64,
) -> Option<CurveGeometry> {
    match geometry {
        CurveGeometry::Nurbs(n) => {
            // Rational NURBS evaluation commutes with affine maps of
            // the control points.
            let mut mapped = n.clone();
            for cp in &mut mapped.control_points {
                *cp = c1 + (*cp - c0) * s;
            }
            Some(CurveGeometry::Nurbs(mapped))
        }
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => Some(CurveGeometry::Circle {
            center: c1 + (*center - c0) * s,
            normal: *normal,
            radius: radius * s,
        }),
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        } => Some(CurveGeometry::Arc {
            center: c1 + (*center - c0) * s,
            normal: *normal,
            radius: radius * s,
            start_angle: *start_angle,
            end_angle: *end_angle,
        }),
        _ => None,
    }
}

/// Recompute the curve of an edge whose adjacent offset surfaces are an
/// axis-perpendicular plane and a coaxial cylinder or cone: the
/// intersection is the exact circle at the plane's station, produced by
/// the radial-scale + axial-translate affine map of the clone curve
/// (parametrization and direction preserved by construction — the
/// PK-303 lesson made direction a contract, not a convention).
///
/// `start_pos` is the recomputed start-vertex position; the candidate
/// must pass through it and lie on both surfaces at sampled parameters,
/// else `None` (keep the clone, decline at the serve gate).
pub(crate) fn coaxial_circle_edge_curve(
    old: &CurveOnEdge,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    start_pos: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    let lin = tolerance.linear();
    let (plane_origin, plane_normal, quadric) = match (surf1, surf2) {
        (SurfaceGeometry::Plane { origin, normal, .. }, q) => {
            (*origin, normal.normalized()?, axis_form(q)?)
        }
        (q, SurfaceGeometry::Plane { origin, normal, .. }) => {
            (*origin, normal.normalized()?, axis_form(q)?)
        }
        _ => return None,
    };
    // Coaxial: plane normal parallel to the quadric axis (structural
    // coaxiality guard, not a geometric tolerance — oblique sections
    // are conics and belong to the general SSI arm, PK-324).
    if plane_normal.cross(&quadric.axis).length() > 1e-9 {
        return None;
    }

    // New circle: axis ∩ plane, radius from the quadric at that station.
    let denom = quadric.axis.dot(&plane_normal);
    if denom.abs() < 1e-12 {
        return None;
    }
    let t_axis = (plane_origin - quadric.origin).dot(&plane_normal) / denom;
    let c1 = quadric.origin + quadric.axis * t_axis;
    let r1 = match quadric.kind {
        AxisKind::Cylinder { radius } => radius,
        AxisKind::Cone { half_angle } => {
            if t_axis <= lin {
                return None; // plane at/behind the apex
            }
            t_axis * half_angle.ptan()
        }
    };
    if r1 <= lin {
        return None;
    }

    // Old circle from the clone curve: sampled center/radius about the
    // same axis, with a 3-sample coherence check (heights and radii
    // must agree — anything else is not an axis-centered circle).
    let (t0, t1) = (old.t_start, old.t_end);
    let samples: Vec<Point3> = [0.0, 0.37, 0.71]
        .iter()
        .map(|f| old.geometry.point_at(t0 + (t1 - t0) * f))
        .collect();
    let heights: Vec<f64> = samples
        .iter()
        .map(|p| (*p - quadric.origin).dot(&quadric.axis))
        .collect();
    let radii: Vec<f64> = samples
        .iter()
        .zip(&heights)
        .map(|(p, h)| (*p - (quadric.origin + quadric.axis * *h)).length())
        .collect();
    if heights.iter().any(|h| (h - heights[0]).abs() > lin)
        || radii.iter().any(|r| (r - radii[0]).abs() > lin)
        || radii[0] <= lin
    {
        return None;
    }
    let c0 = quadric.origin + quadric.axis * heights[0];
    let scale = r1 / radii[0];

    let geometry = map_coaxial_circle_geometry(&old.geometry, c0, c1, scale)?;
    let candidate = CurveOnEdge {
        geometry,
        t_start: t0,
        t_end: t1,
    };

    // Verification: the candidate must pass through the recomputed
    // start vertex and lie on both offset surfaces across the span.
    if candidate
        .geometry
        .point_at(candidate.t_start)
        .distance_to(&start_pos)
        > lin
    {
        return None;
    }
    for f in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let p = candidate.geometry.point_at(t0 + (t1 - t0) * f);
        for surf in [surf1, surf2] {
            if analytic_residual(surf, &p).is_none_or(|r| r > lin) {
                return None;
            }
        }
    }
    Some(candidate)
}

/// Recompute a seam edge (both half-edges on the same face) whose clone
/// curve is a straight ruling of a cylinder or cone: the recomputed
/// seam is the chord between the recomputed endpoint vertices, verified
/// to lie on the offset surface at its midpoint (only a genuine ruling
/// passes that on a quadric). Direction follows the vertex order, the
/// same convention as the plane×plane arm.
pub(crate) fn straight_ruling_seam_curve(
    old: &CurveOnEdge,
    surf: &SurfaceGeometry,
    start_pos: Point3,
    end_pos: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    if !matches!(old.geometry, CurveGeometry::Line { .. }) {
        return None;
    }
    if !matches!(
        surf,
        SurfaceGeometry::Cylinder { .. } | SurfaceGeometry::Cone { .. }
    ) {
        return None;
    }
    let lin = tolerance.linear();
    let chord = end_pos - start_pos;
    if chord.length() <= lin {
        return None;
    }
    let mid = start_pos + chord * 0.5;
    if analytic_residual(surf, &mid).is_none_or(|r| r > lin) {
        return None;
    }
    if analytic_residual(surf, &start_pos).is_none_or(|r| r > lin)
        || analytic_residual(surf, &end_pos).is_none_or(|r| r > lin)
    {
        return None;
    }
    Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: start_pos,
            direction: chord,
        },
        t_start: 0.0,
        t_end: 1.0,
    })
}
