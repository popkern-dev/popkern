//! Exact geometry authoring for arc joins (design §5 Stage 3a, PK-325).
//!
//! The rolling ball centered on a gap edge touches both offset surfaces
//! at the offset distance, so every join surface is anchored on ORIGINAL
//! geometry: strips are canal surfaces around the original edge
//! (Cylinder for straight spines, Torus for circular ones, a fillet
//! rational sweep for general spines), rails are the edge curve mapped
//! through the adjacent face's normal field, and corner arcs are great
//! circles on the corner sphere.

use crate::curved_recompute::project_to_surface;
use crate::curved_recompute_similarity::map_similarity_geometry;
use popkern_core::geometry::analytic_curves::circle_basis;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;

/// Sample fractions used for rail fitting and verification.
const RAIL_FRACTIONS: [f64; 9] = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0];

/// Map the clone edge curve through the offset surface's normal field:
/// the rail is the exact projection image of the spine. Tries the exact
/// similarity form first (planes translate; sphere/cylinder/cone circle
/// rails are center/axis-anchored similarities); falls back to a NURBS
/// fit through the projections. Returns the rail and its measured
/// worst-sample deviation from the offset surface.
pub(crate) fn offset_rail(
    clone: &CurveOnEdge,
    offset: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Option<(CurveOnEdge, f64)> {
    let ts: Vec<f64> = RAIL_FRACTIONS
        .iter()
        .map(|f| clone.t_start + f * (clone.t_end - clone.t_start))
        .collect();
    let ps: Vec<Point3> = ts.iter().map(|t| clone.geometry.point_at(*t)).collect();
    let qs: Vec<Point3> = ps
        .iter()
        .map(|p| project_to_surface(offset, p))
        .collect::<Option<Vec<_>>>()?;

    // Exact similarity attempt: anchor + scale from the widest pair,
    // witnessed on every sample.
    let gate = tolerance.linear() * 0.5;
    if let Some(rail) = similarity_rail(clone, &ps, &qs, gate) {
        return Some((rail, 0.0));
    }

    // General rail: NURBS through the projections, deviation measured
    // at fit midpoints against the offset surface.
    let nurbs = popkern_fillet::center_curve::fit_nurbs_through_points(&qs).ok()?;
    let (t0, t1) = nurbs.parameter_range();
    let mut worst = 0.0f64;
    for i in 0..16 {
        let t = t0 + (i as f64 + 0.5) / 16.0 * (t1 - t0);
        let p = nurbs.point_at(t);
        let d = crate::curved_recompute::surface_distance(offset, &p)?;
        worst = worst.max(d);
    }
    Some((
        CurveOnEdge {
            geometry: CurveGeometry::Nurbs(nurbs),
            t_start: t0,
            t_end: t1,
        },
        worst,
    ))
}

fn similarity_rail(
    clone: &CurveOnEdge,
    ps: &[Point3],
    qs: &[Point3],
    gate: f64,
) -> Option<CurveOnEdge> {
    // Widest sample pair sets the scale (PK-291 lesson: never a
    // near-coincident pair).
    let mut best = (0usize, 1usize, 0.0f64);
    for i in 0..ps.len() {
        for j in (i + 1)..ps.len() {
            let d = ps[i].distance_to(&ps[j]);
            if d > best.2 {
                best = (i, j, d);
            }
        }
    }
    if best.2 < 1.0e-12 {
        return None;
    }
    let s = qs[best.0].distance_to(&qs[best.1]) / best.2;
    let anchor = Vec3::new(
        qs[best.0].x - s * ps[best.0].x,
        qs[best.0].y - s * ps[best.0].y,
        qs[best.0].z - s * ps[best.0].z,
    );
    for (p, q) in ps.iter().zip(qs.iter()) {
        let mapped = Point3::new(anchor.x + s * p.x, anchor.y + s * p.y, anchor.z + s * p.z);
        if mapped.distance_to(q) > gate {
            return None;
        }
    }
    let geometry = map_similarity_geometry(&clone.geometry, anchor, s)?;
    Some(CurveOnEdge {
        geometry,
        t_start: clone.t_start,
        t_end: clone.t_end,
    })
}

/// The canal (rolling-ball envelope) surface around the spine at the
/// given radius. Analytic promotion: straight spine → Cylinder,
/// circular spine → Torus (declining spindle tori); general spine →
/// the fillet substrate's rational sweep driven with synthetic inputs.
pub(crate) fn strip_surface(
    spine: &CurveOnEdge,
    radius: f64,
    off1: &SurfaceGeometry,
    off2: &SurfaceGeometry,
    rail1: &CurveOnEdge,
    rail2: &CurveOnEdge,
) -> Result<SurfaceGeometry, String> {
    match &spine.geometry {
        CurveGeometry::Line { origin, direction } => {
            let axis = direction
                .normalized()
                .ok_or("degenerate straight spine direction")?;
            Ok(SurfaceGeometry::Cylinder {
                origin: *origin,
                axis,
                radius,
            })
        }
        CurveGeometry::Circle {
            center,
            normal,
            radius: major,
        }
        | CurveGeometry::Arc {
            center,
            normal,
            radius: major,
            ..
        } => {
            if radius >= *major - 1.0e-12 {
                return Err(format!(
                    "spindle torus: ball radius {radius} >= spine radius {major}"
                ));
            }
            Ok(SurfaceGeometry::Torus {
                center: *center,
                axis: *normal,
                major_radius: *major,
                minor_radius: radius,
            })
        }
        _ => general_canal_surface(spine, radius, off1, off2, rail1, rail2),
    }
}

/// General spine: drive `popkern_fillet::build_fillet_surface` with a
/// synthetic center-curve result — the spine IS the rolling-ball center
/// curve (no solve), the rails are the tangency footprints.
fn general_canal_surface(
    spine: &CurveOnEdge,
    radius: f64,
    off1: &SurfaceGeometry,
    off2: &SurfaceGeometry,
    rail1: &CurveOnEdge,
    rail2: &CurveOnEdge,
) -> Result<SurfaceGeometry, String> {
    const N: usize = 20;
    let sample = |c: &CurveOnEdge, i: usize| {
        let t = c.t_start + i as f64 / (N - 1) as f64 * (c.t_end - c.t_start);
        c.geometry.point_at(t)
    };
    let center_pts: Vec<Point3> = (0..N).map(|i| sample(spine, i)).collect();
    let footprint_1: Vec<Point3> = (0..N).map(|i| sample(rail1, i)).collect();
    let footprint_2: Vec<Point3> = (0..N).map(|i| sample(rail2, i)).collect();
    let center_curve = popkern_fillet::center_curve::fit_nurbs_through_points(&center_pts)?;
    let center = popkern_fillet::center_curve::CenterCurveResult {
        center_curve,
        sweep_center_curve: None,
        sweep_footprint_1: None,
        sweep_footprint_2: None,
        footprint_1,
        footprint_2,
        contact_solutions: None,
    };
    let result = popkern_fillet::fillet_surface::build_fillet_surface(
        off1,
        off2,
        &center,
        radius,
        8,
        true,
        popkern_fillet::ContinuityLevel::G1,
        false,
        false,
        popkern_fillet::types::BlendProfile::Circular,
    )?;
    Ok(result.surface)
}

/// The minor great-circle arc from `a` to `b` on the sphere
/// `(center, radius)`, authored in the canonical [`circle_basis`] frame
/// so the enum evaluation reproduces the endpoints exactly.
pub(crate) fn great_arc(
    center: Point3,
    radius: f64,
    a: Point3,
    b: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    let ra = a - center;
    let rb = b - center;
    let normal = ra.cross(&rb).normalized()?;
    let (u, v) = circle_basis(normal);
    let theta_a = ra.dot(&v).atan2(ra.dot(&u));
    let mut theta_b = rb.dot(&v).atan2(rb.dot(&u));
    while theta_b <= theta_a {
        theta_b += std::f64::consts::TAU;
    }
    if theta_b - theta_a >= std::f64::consts::PI + 1.0e-9 {
        // The cross-product normal guarantees the CCW arc a→b is the
        // minor one; reaching here means a/b are near-antipodal.
        return None;
    }
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle: theta_a,
            end_angle: theta_b,
        },
        t_start: theta_a,
        t_end: theta_b,
    };
    // Basis self-check: the authored angles must reproduce the inputs.
    let gate = (tolerance.linear() * 10.0).max(radius * 1.0e-9);
    if curve.geometry.point_at(theta_a).distance_to(&a) > gate
        || curve.geometry.point_at(theta_b).distance_to(&b) > gate
    {
        return None;
    }
    Some(curve)
}

/// Geometric (∂u×∂v) surface normal at the parametric image of `p`.
/// Defers to the surface's own uv inverse + `normal_at` so the SIGN is
/// the parametrization's, never a hand-derived guess (orientation
/// flags are measured against exactly this normal).
pub(crate) fn surface_geo_normal(surface: &SurfaceGeometry, p: &Point3) -> Option<Vec3> {
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, p)?;
    surface.normal_at(u, v).normalized()
}

#[cfg(test)]
#[path = "arc_join_geometry_tests.rs"]
mod tests;
