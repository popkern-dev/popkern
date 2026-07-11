//! Analytic plane-sphere intersection solver.
//!
//! When the two surfaces are a plane and a sphere the intersection is either
//! empty (no contact), tangent (degenerate — treated as empty), or a circle.
//! Computing this analytically is far cheaper than the general NURBS marching
//! pipeline.

use std::f64::consts::FRAC_1_SQRT_2;

use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Try to analytically intersect a plane and a sphere.
///
/// Returns `Some(Vec<SurfaceSurfaceIntersection>)` when the surface pair is
/// a plane and a sphere (in either order), or `None` when neither surface is
/// a plane-sphere combination and the caller should fall back to the general
/// intersection pipeline.
///
/// Returns `Some(vec![])` for no intersection or tangent (degenerate) cases.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_plane_sphere(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    match (surface1, surface2) {
        (
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            },
            SurfaceGeometry::Sphere { center, radius },
        ) => Some(intersect_plane_sphere(
            *origin, *normal, *u_axis, *v_axis, u1_range, v1_range, *center, *radius, u2_range,
            v2_range, tolerance, false,
        )),
        (
            SurfaceGeometry::Sphere { center, radius },
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            },
        ) => Some(intersect_plane_sphere(
            *origin, *normal, *u_axis, *v_axis, u2_range, v2_range, *center, *radius, u1_range,
            v1_range, tolerance, true,
        )),
        _ => None,
    }
}

/// Core implementation: intersect a plane with a sphere.
///
/// When `swapped` is true the sphere was surface 1 and the plane was surface 2,
/// so we swap the UV curves in the output.
#[allow(clippy::too_many_arguments)]
fn intersect_plane_sphere(
    plane_origin: Point3,
    plane_normal: Vec3,
    plane_u: Vec3,
    plane_v: Vec3,
    plane_u_range: (f64, f64),
    plane_v_range: (f64, f64),
    sphere_center: Point3,
    sphere_radius: f64,
    sphere_u_range: (f64, f64),
    sphere_v_range: (f64, f64),
    tolerance: &Tolerance,
    swapped: bool,
) -> Vec<SurfaceSurfaceIntersection> {
    // Signed distance from sphere centre to the plane.
    let diff = sphere_center - plane_origin;
    let d = diff.dot(&plane_normal);

    // No intersection: sphere is too far from plane.
    if d.abs() > sphere_radius + tolerance.linear() {
        return vec![];
    }

    // Tangent (degenerate): intersection would be a single point.
    if d.abs() > sphere_radius - tolerance.linear() {
        return vec![];
    }

    // Circle radius and projected centre.
    let r_circle = (sphere_radius * sphere_radius - d * d).sqrt();
    let center_proj = sphere_center - plane_normal * d;

    // Build in-plane axes perpendicular to the plane normal.
    // We use the plane's own u/v axes as our circle frame.
    let axis1 = plane_u;
    let axis2 = plane_v;

    // 3D NURBS circle.
    let curve_3d = nurbs_circle_3d(center_proj, axis1, axis2, r_circle);

    // UV on plane: project each 3D control point of the circle into the
    // plane's UV space. This correctly handles non-orthogonal u_axis/v_axis
    // (which produce skewed UV projections, not axis-aligned ellipses).
    let curve_uv_plane = project_curve_to_plane_uv(&curve_3d, plane_origin, plane_u, plane_v);

    // UV on sphere: co-parameterized with `curve_3d` — the pcurve must
    // satisfy S_sphere(uv(t)) = curve_3d(t) for the SAME t. A tilted
    // plane's circle is a transcendental curve in the (longitude,
    // colatitude) chart, so it is fitted: chart the circle at the 3D
    // curve's own parameters, unwrap the longitude by continuity, and
    // interpolate with a degree-1 NURBS whose knots ARE those
    // parameters, refining until mid-knot agreement passes. (The old
    // placeholder was a constant-colatitude u-line: geometrically
    // wrong for any plane not perpendicular to the polar axis, and not
    // co-parameterized even when perpendicular — it only survived
    // because downstream chart embeddings 3D-verify pcurves and fall
    // back to projection.)
    let curve_uv_sphere = fit_sphere_uv_curve(&curve_3d, sphere_center, sphere_radius, tolerance);

    // Respect the UV ranges — if the circle centre is outside the plane's
    // parameter box, the intersection may not be relevant, but we still
    // return it (clipping is the caller's responsibility).
    let _ = (plane_u_range, plane_v_range, sphere_u_range, sphere_v_range);

    let (curve_uv1, curve_uv2) = if swapped {
        (curve_uv_sphere, curve_uv_plane)
    } else {
        (curve_uv_plane, curve_uv_sphere)
    };

    let cos = CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed: true,
    };

    vec![SurfaceSurfaceIntersection::Curve(cos)]
}

/// Build a degree-2 rational NURBS circle with 9 control points.
///
/// The circle lies in the plane spanned by `axis1` and `axis2`, centred at
/// `center`, with the given `radius`.
fn nurbs_circle_3d(center: Point3, axis1: Vec3, axis2: Vec3, radius: f64) -> NurbsCurve {
    nurbs_ellipse_3d(center, axis1, axis2, radius, radius)
}

/// Build a degree-2 rational NURBS ellipse in 3D.
fn nurbs_ellipse_3d(center: Point3, axis1: Vec3, axis2: Vec3, ra: f64, rb: f64) -> NurbsCurve {
    let w = FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            center + axis1 * ra,
            center + axis1 * ra + axis2 * rb,
            center + axis2 * rb,
            center - axis1 * ra + axis2 * rb,
            center - axis1 * ra,
            center - axis1 * ra - axis2 * rb,
            center - axis2 * rb,
            center + axis1 * ra - axis2 * rb,
            center + axis1 * ra,
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Fit the sphere-chart image of `curve_3d` as a co-parameterized
/// degree-1 NURBS: `S_sphere(fit(t)) ≈ curve_3d(t)` for the same `t`.
///
/// Samples the 3D circle at its own parameters, charts each sample
/// (`u = atan2`, `v = acos` — the kernel sphere convention), unwraps
/// the longitude by continuity so seam crossings stay a continuous
/// polyline (period branch ±2π), and interpolates with knots AT the
/// sample parameters. Sample density doubles until every mid-knot
/// deviation passes the linear tolerance or the 4097-point cap is
/// reached (chord sag shrinks quadratically per doubling, so the cap
/// is a safeguard; at default 1e-6 a tilted great-ish circle needs
/// ~4k points — a recorded candidate for a piecewise-quadratic
/// upgrade if profiling ever flags it).
///
/// Pole guard: a sample within the chart singularity keeps its
/// neighbor's longitude (`atan2(0, 0)` is meaningless there and any
/// value evaluates back to the pole).
fn fit_sphere_uv_curve(
    curve_3d: &NurbsCurve,
    sphere_center: Point3,
    sphere_radius: f64,
    tolerance: &Tolerance,
) -> NurbsCurve {
    use crate::geometry::curve::CurveEval;
    use std::f64::consts::TAU;

    let gate = tolerance.linear().max(1.0e-9);
    let chart = |p: Point3, prev_u: Option<f64>| -> (f64, f64) {
        let d = (p - sphere_center) * (1.0 / sphere_radius);
        let v = d.z.clamp(-1.0, 1.0).acos();
        let planar = (d.x * d.x + d.y * d.y).sqrt();
        let mut u = if planar < 1.0e-12 {
            prev_u.unwrap_or(0.0)
        } else {
            let raw = d.y.atan2(d.x);
            if raw < 0.0 { raw + TAU } else { raw }
        };
        if let Some(prev) = prev_u {
            // Continuity unwrap: stay within half a period of the
            // previous sample.
            u -= TAU * ((u - prev) / TAU).round();
        }
        (u, v)
    };

    let mut n = 65usize;
    loop {
        let mut params = Vec::with_capacity(n);
        let mut uvs: Vec<(f64, f64)> = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let prev_u = uvs.last().map(|&(u, _)| u);
            let uv = chart(curve_3d.point_at(t), prev_u);
            params.push(t);
            uvs.push(uv);
        }
        // Mid-knot verification: the polyline evaluated halfway
        // between knots must land back on the 3D circle.
        let mut worst = 0.0f64;
        for i in 0..n - 1 {
            let t_mid = 0.5 * (params[i] + params[i + 1]);
            let (u, v) = (
                0.5 * (uvs[i].0 + uvs[i + 1].0),
                0.5 * (uvs[i].1 + uvs[i + 1].1),
            );
            let on_sphere = sphere_center
                + Vec3::new(
                    sphere_radius * v.sin() * u.cos(),
                    sphere_radius * v.sin() * u.sin(),
                    sphere_radius * v.cos(),
                );
            worst = worst.max(on_sphere.distance_to(&curve_3d.point_at(t_mid)));
        }
        if worst <= gate || n >= 4097 {
            // Clamped degree-1 knot vector with knots at the sample
            // parameters: the curve passes through sample i at t_i,
            // preserving co-parameterization with curve_3d.
            let mut knots = Vec::with_capacity(n + 2);
            knots.push(params[0]);
            knots.extend(params.iter().copied());
            knots.push(params[n - 1]);
            return NurbsCurve {
                degree: 1,
                knots,
                control_points: uvs
                    .into_iter()
                    .map(|(u, v)| Point3::new(u, v, 0.0))
                    .collect(),
                weights: vec![1.0; n],
            };
        }
        n = (n - 1) * 2 + 1;
    }
}

/// Project a 3D NURBS curve onto a plane's UV space by transforming each
/// control point. Preserves degree, knots, and weights.
///
/// Uses the correct inverse parameterization:
/// `u = (p-origin)·u_axis / |u_axis|²`, `v = (p-origin)·v_axis / |v_axis|²`.
fn project_curve_to_plane_uv(
    curve: &NurbsCurve,
    origin: Point3,
    u_axis: Vec3,
    v_axis: Vec3,
) -> NurbsCurve {
    let u_len_sq = u_axis.dot(&u_axis);
    let v_len_sq = v_axis.dot(&v_axis);
    let cps = curve
        .control_points
        .iter()
        .map(|p| {
            let d = *p - origin;
            let u = if u_len_sq > 1e-30 {
                d.dot(&u_axis) / u_len_sq
            } else {
                d.dot(&u_axis)
            };
            let v = if v_len_sq > 1e-30 {
                d.dot(&v_axis) / v_len_sq
            } else {
                d.dot(&v_axis)
            };
            Point3::new(u, v, 0.0)
        })
        .collect();
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: cps,
        weights: curve.weights.clone(),
    }
}

#[cfg(test)]
#[path = "plane_sphere_tests.rs"]
mod tests;
