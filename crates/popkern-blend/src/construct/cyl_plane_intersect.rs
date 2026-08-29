//! Shared cylinder × plane analytic intersection helpers.
//!
//! Used by:
//! * `topology_repair` — to rewrite chord-line edges as analytic
//!   cyl×plane arcs after the fact.
//! * `helpers::compute_perp_anchor_positions_for_face` — to place
//!   anchor vertices ON the cyl×plane intersection at the anchor-split stage,
//!   rather than at perp-walk positions that drift off-cylinder
//!   for non-perpendicular faces (e.g. pyramid slants).  This is
//!   the commercial-grade fix for #136 — anchors are on-surface
//!   from the moment they're created, so the cylinder face's
//!   outer loop never contains chord-line edges through interior.

use popkern_core::geometry::Point3;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, intersect_surfaces,
};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;

/// Compute the analytic cyl×plane intersection, returning the first
/// `Curve` intersection from `intersect_surfaces`. Handles the
/// fillet-specific tangent case (cyl × slant plane = contact line)
/// that `case_parallel` rejects when `d ≈ cr`.
pub(crate) fn analytic_cyl_plane_intersection(
    cyl: &SurfaceGeometry,
    plane: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Option<NurbsCurve> {
    // Cyl UV: u ∈ [0, 2π], v unbounded — give the solver a generous v
    // range. Plane UV: use a reasonably wide square. The analytic
    // plane×cylinder solver doesn't actually clip the result to these
    // ranges, so generous bounds are safe.
    let result = intersect_surfaces(
        cyl,
        (0.0, std::f64::consts::TAU),
        (-10.0, 10.0),
        plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        tolerance,
    )
    .ok()?;
    for inter in &result.intersections {
        if let SurfaceSurfaceIntersection::Curve(c) = inter {
            return Some(c.curve_3d.clone());
        }
    }
    // Empty result — could be a tangent case (cyl perpendicular to
    // plane normal AND axis distance equals radius). Compute the
    // contact line manually.
    cyl_plane_contact_line(cyl, plane, tolerance)
}

/// Compute the contact line where a cylinder is tangent to a plane.
/// Returns a NURBS line along the cylinder axis at the contact angle.
fn cyl_plane_contact_line(
    cyl: &SurfaceGeometry,
    plane: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Option<NurbsCurve> {
    let SurfaceGeometry::Cylinder {
        origin: co,
        axis: ca,
        radius: cr,
    } = cyl
    else {
        return None;
    };
    let SurfaceGeometry::Plane {
        origin: po,
        normal: pn,
        ..
    } = plane
    else {
        return None;
    };
    // Tangent requires axis perpendicular to plane normal (within
    // angular tolerance) AND axial-distance from cyl axis to plane
    // ≈ cyl radius.
    let axis_dot_normal = ca.dot(pn);
    let ca_len = ca.length();
    let pn_len = pn.length();
    if (axis_dot_normal / (ca_len * pn_len)).abs() > tolerance.angular() {
        return None; // not parallel
    }
    let axial_dist = ((*co - *po).dot(pn)) / pn_len;
    if (axial_dist.abs() - cr).abs() > cr * 0.05 {
        return None;
    }
    // Direction from cyl axis toward plane (negative of normal·sign).
    let pn_unit = Vec3::new(pn.x / pn_len, pn.y / pn_len, pn.z / pn_len);
    let toward_plane_sign = -axial_dist.signum();
    let radial = Vec3::new(
        toward_plane_sign * pn_unit.x * *cr,
        toward_plane_sign * pn_unit.y * *cr,
        toward_plane_sign * pn_unit.z * *cr,
    );
    let line_origin = Point3::new(co.x + radial.x, co.y + radial.y, co.z + radial.z);
    let v_lo = -10.0;
    let v_hi = 10.0;
    let a = Point3::new(
        line_origin.x + ca.x * v_lo,
        line_origin.y + ca.y * v_lo,
        line_origin.z + ca.z * v_lo,
    );
    let b = Point3::new(
        line_origin.x + ca.x * v_hi,
        line_origin.y + ca.y * v_hi,
        line_origin.z + ca.z * v_hi,
    );
    Some(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    })
}

/// Find the parameter `s ∈ (eps, 1 − eps)` along the line segment
/// `p0 → p1` where the curve crosses, or `None` if no crossing.
///
/// Implementation: sample the curve densely, project each sample
/// onto the line segment, and look for the segment-side crossing
/// (signed distance from segment changes sign).  Refines via
/// bisection between bracketing samples.
///
/// Returns the 3-D crossing point along with `s` so callers can
/// both place anchors and split BRep edges at the correct fraction.
pub(crate) fn intersect_curve_with_segment(
    curve: &NurbsCurve,
    p0: Point3,
    p1: Point3,
) -> Option<(f64, Point3)> {
    let dir = Vec3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
    let len_sq = dir.dot(&dir);
    if len_sq < 1e-18 {
        return None;
    }
    // Project a 3-D point onto the line and return `s` (the
    // segment-fraction projection) and the perpendicular distance.
    let project = |p: Point3| -> (f64, f64) {
        let rel = Vec3::new(p.x - p0.x, p.y - p0.y, p.z - p0.z);
        let s = rel.dot(&dir) / len_sq;
        let foot = Point3::new(p0.x + s * dir.x, p0.y + s * dir.y, p0.z + s * dir.z);
        let perp = Vec3::new(p.x - foot.x, p.y - foot.y, p.z - foot.z);
        (s, perp.dot(&perp).sqrt())
    };

    // Sample the curve and find the sample with smallest perp
    // distance.  We DON'T require `s ∈ [0, 1]` for the seed — the
    // closest sample may have its projection just outside the
    // segment, yet a near-by t still has a valid in-segment
    // crossing.  We also DON'T reject by absolute perp magnitude
    // here — for a degree-3 NURBS fit through footprint points, a
    // 64-sample grid leaves ~7 mm gaps, so the seed can be that
    // far from the line even though Newton converges to a true
    // crossing.  Newton refines first; final acceptance is
    // perp < 1e-6 m AND s ∈ (0.001, 0.999).
    let (t_lo, t_hi) = curve.parameter_range();
    let n_samples: usize = 256;
    let mut best_perp = f64::INFINITY;
    let mut t_seed = t_lo;
    for i in 0..=n_samples {
        let t = t_lo + (t_hi - t_lo) * (i as f64) / (n_samples as f64);
        let p = curve.point_at(t);
        let (_, perp) = project(p);
        if perp < best_perp {
            best_perp = perp;
            t_seed = t;
        }
    }

    // Newton refine: minimise perp distance starting from t_seed.
    // For analytic NURBS the intersection curve is smooth; a few
    // Newton steps (with bracket clamping) converge quickly.
    let mut t = t_seed;
    let dt = (t_hi - t_lo) * 1e-6;
    let mut converged_perp = best_perp;
    for _ in 0..32 {
        let p_t = curve.point_at(t);
        let (_, perp_t) = project(p_t);
        converged_perp = perp_t;
        if perp_t < 1e-9 {
            break;
        }
        let p_plus = curve.point_at(t + dt);
        let (_, perp_plus) = project(p_plus);
        let dperp = (perp_plus - perp_t) / dt;
        if dperp.abs() < 1e-12 {
            break;
        }
        let step = (perp_t / dperp).clamp(-(t_hi - t_lo) * 0.1, (t_hi - t_lo) * 0.1);
        t = (t - step).clamp(t_lo, t_hi);
    }
    let p_final = curve.point_at(t);
    let (s_final, perp_final) = project(p_final);
    if (0.0..=1.0).contains(&s_final) && perp_final < 1e-6 {
        Some((s_final, p_final))
    } else {
        // Fallback: report final state but only if it's at least
        // a reasonable seed match (within 5% of segment length).
        let seg_len = len_sq.sqrt();
        let _ = converged_perp;
        if (0.0..=1.0).contains(&s_final) && perp_final < seg_len * 0.05 {
            Some((s_final, p_final))
        } else {
            None
        }
    }
}
