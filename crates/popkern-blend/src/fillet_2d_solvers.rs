//! Closed-form and numerical fillet-center solvers used by
//! [`crate::fillet_2d_plan::corner_plan`].
//!
//! Three solvers live here:
//!
//! - [`solve_line_arc_fillet`] — offset-line-offset-circle
//!   intersection. Closed form via a scalar quadratic in the
//!   line parameter.
//! - [`solve_arc_arc_fillet`] — offset-circle-offset-circle
//!   intersection. Closed form via classical two-circle
//!   intersection geometry.
//! - [`solve_generic_fillet`] — numerical Newton solver that
//!   works on **any** two planar curves (Line, Arc, NURBS). Used
//!   whenever at least one adjacent edge is a NURBS, since
//!   offset curves of NURBS aren't NURBS in general. Converges
//!   from a bisector-based initial guess using the fact that
//!   `∂|F − closest_point(F)| / ∂F = outward_unit_normal`, giving
//!   a 2×2 Newton step in the face-plane coordinate frame.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// Solve a Line-Arc fillet via offset-line-offset-circle
/// intersection. Returns `(line_tp, arc_tp, fillet_center)`.
///
/// The fillet centre is on the line offset by `radius` on the
/// bisector side, **and** on the circle concentric with the arc
/// at radius `arc_radius ± radius` (minus for internal tangency
/// when the corner sits on the arc's concave side, plus
/// otherwise). Substituting the offset-line parameterisation
/// into the offset-circle equation gives a scalar quadratic in
/// `t`; the smaller non-negative root is the correct setback.
pub(crate) fn solve_line_arc_fillet(
    corner_pos: Point3,
    line_tangent: Vec3,
    arc_center: Point3,
    arc_radius: f64,
    bisector: Vec3,
    radius: f64,
) -> Result<(Point3, Point3, Point3), String> {
    let proj = line_tangent * bisector.dot(&line_tangent);
    let n_line = (bisector - proj)
        .normalized()
        .ok_or_else(|| "line interior normal degenerate".to_string())?;

    let to_center = arc_center - corner_pos;
    let internal = bisector.dot(&to_center) > 0.0;
    let r_off = if internal {
        arc_radius - radius
    } else {
        arc_radius + radius
    };
    if r_off <= 1.0e-9 {
        return Err(format!(
            "fillet radius {radius:.6} too large for internal tangency on arc radius {arc_radius:.6}"
        ));
    }

    let v = (corner_pos + n_line * radius) - arc_center;
    let b = v.dot(&line_tangent);
    let c = v.length_squared() - r_off * r_off;
    let disc = b * b - c;
    if disc < 0.0 {
        return Err(format!(
            "fillet geometry has no real solution (discriminant = {disc:.6})"
        ));
    }
    let sqrt_d = disc.sqrt();
    let t_plus = -b + sqrt_d;
    let t_minus = -b - sqrt_d;

    let t = if t_minus >= -1.0e-9 {
        t_minus.max(0.0)
    } else if t_plus >= -1.0e-9 {
        t_plus.max(0.0)
    } else {
        return Err(format!(
            "both fillet-center line-parameter roots lie behind the corner (t_minus={t_minus:.6}, t_plus={t_plus:.6})"
        ));
    };

    let fillet_center = corner_pos + n_line * radius + line_tangent * t;
    let line_tangent_point = corner_pos + line_tangent * t;

    let to_fc = fillet_center - arc_center;
    let to_fc_len = to_fc.length();
    if to_fc_len < 1.0e-12 {
        return Err("fillet center coincides with arc center".into());
    }
    let arc_tangent_point = arc_center + to_fc * (arc_radius / to_fc_len);

    Ok((line_tangent_point, arc_tangent_point, fillet_center))
}

/// Solve an Arc-Arc fillet via offset-circle-offset-circle
/// intersection. Returns `(tangent_on_arc_a, tangent_on_arc_b,
/// fillet_center)`. Closed form; picks the intersection root
/// closest to the corner.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_arc_arc_fillet(
    corner_pos: Point3,
    arc_a_center: Point3,
    arc_a_radius: f64,
    arc_b_center: Point3,
    arc_b_radius: f64,
    bisector: Vec3,
    plane_normal: Vec3,
    radius: f64,
) -> Result<(Point3, Point3, Point3), String> {
    let to_ca = arc_a_center - corner_pos;
    let internal_a = bisector.dot(&to_ca) > 0.0;
    let r_a_off = if internal_a {
        arc_a_radius - radius
    } else {
        arc_a_radius + radius
    };
    if r_a_off <= 1.0e-9 {
        return Err(format!(
            "fillet radius {radius:.6} too large for tangency on arc A radius {arc_a_radius:.6}"
        ));
    }
    let to_cb = arc_b_center - corner_pos;
    let internal_b = bisector.dot(&to_cb) > 0.0;
    let r_b_off = if internal_b {
        arc_b_radius - radius
    } else {
        arc_b_radius + radius
    };
    if r_b_off <= 1.0e-9 {
        return Err(format!(
            "fillet radius {radius:.6} too large for tangency on arc B radius {arc_b_radius:.6}"
        ));
    }

    let d_vec = arc_b_center - arc_a_center;
    let d = d_vec.length();
    if d < 1.0e-12 {
        return Err("arc centres coincide; arc-arc fillet requires distinct circles".into());
    }
    if d > r_a_off + r_b_off + 1.0e-9 {
        return Err(format!(
            "offset circles too far apart to intersect (d={d:.6}, r_sum={:.6})",
            r_a_off + r_b_off
        ));
    }
    if d + 1.0e-9 < (r_a_off - r_b_off).abs() {
        return Err(format!(
            "one offset circle strictly contains the other (d={d:.6}, |r_diff|={:.6})",
            (r_a_off - r_b_off).abs()
        ));
    }

    let a = (d * d + r_a_off * r_a_off - r_b_off * r_b_off) / (2.0 * d);
    let h_sq = r_a_off * r_a_off - a * a;
    let h = if h_sq < 0.0 {
        if h_sq > -1.0e-9 {
            0.0
        } else {
            return Err(format!(
                "offset-circle intersection discriminant {h_sq:.6} < 0"
            ));
        }
    } else {
        h_sq.sqrt()
    };

    let d_unit = d_vec * (1.0 / d);
    let p_mid = arc_a_center + d_unit * a;
    let perp = plane_normal
        .cross(&d_unit)
        .normalized()
        .ok_or_else(|| "in-plane perpendicular degenerate".to_string())?;

    let f_plus = p_mid + perp * h;
    let f_minus = p_mid - perp * h;

    let fillet_center =
        if (f_plus - corner_pos).length_squared() <= (f_minus - corner_pos).length_squared() {
            f_plus
        } else {
            f_minus
        };

    let to_fc_a = fillet_center - arc_a_center;
    let len_a = to_fc_a.length();
    if len_a < 1.0e-12 {
        return Err("fillet centre coincides with arc A centre".into());
    }
    let tangent_a = arc_a_center + to_fc_a * (arc_a_radius / len_a);

    let to_fc_b = fillet_center - arc_b_center;
    let len_b = to_fc_b.length();
    if len_b < 1.0e-12 {
        return Err("fillet centre coincides with arc B centre".into());
    }
    let tangent_b = arc_b_center + to_fc_b * (arc_b_radius / len_b);

    Ok((tangent_a, tangent_b, fillet_center))
}

/// Generic numerical fillet-center solver that handles **any**
/// pair of planar curves by Newton iteration in the face plane.
///
/// The fillet centre `F` satisfies `|F − closest_point_a(F)| =
/// radius` and `|F − closest_point_b(F)| = radius`. Differen­tiating
/// the distance-to-curve gives the outward unit normal at the
/// projection — so the 2×2 Jacobian in any face-plane basis
/// `(u, v)` has rows `[n_a · u, n_a · v]` and `[n_b · u, n_b · v]`.
///
/// Initial guess: `corner_pos + bisector · (radius / sin(θ/2))`,
/// matching the Line-Line closed-form result. For small radii on
/// smooth curves Newton converges in a handful of iterations.
///
/// ## Local projection for NURBS
///
/// For non-convex NURBS curves the global closest-point search
/// may return a point far from the corner (a spurious global
/// minimum). This solver uses [`local_closest_point`] for NURBS
/// curves, seeded from the corner's projection parameter, so
/// the projection stays in the neighbourhood of the corner
/// regardless of the curve's global shape.
///
/// ## Wrong-side guard
///
/// The unsigned-distance formulation converges to whichever
/// side of each curve the initial guess lies on. For pathological
/// geometry (e.g. highly curved NURBS with radius close to the
/// curvature limit) Newton can drift across the corner and
/// converge to the exterior-tangency root. After convergence we
/// compute the **interior perpendicular** of each edge at the
/// corner — the in-plane vector perpendicular to `dir_edge`
/// that has positive dot with `bisector` — and verify that
/// `(F − corner_pos) · n_interior > 0` for both edges. This is
/// tighter than a bisector-only check and catches both Newton
/// drift and misaligned initial guesses.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_generic_fillet(
    curve_a: &CurveGeometry,
    curve_b: &CurveGeometry,
    corner_pos: Point3,
    dir_a: Vec3,
    dir_b: Vec3,
    bisector: Vec3,
    plane_normal: Vec3,
    radius: f64,
    corner_angle: f64,
) -> Result<(Point3, Point3, Point3), String> {
    let (u_basis, v_basis) = plane_basis(plane_normal);

    // Pre-compute corner-anchored parameters on each curve. For
    // NURBS these anchor the local projection search; for
    // Line/Arc they're unused (trait `closest_point` is already
    // a local-minimum solver thanks to their convex geometry).
    let t_corner_a = corner_parameter_hint(curve_a, corner_pos);
    let t_corner_b = corner_parameter_hint(curve_b, corner_pos);

    let half_angle = corner_angle * 0.5;
    let sin_half = half_angle.psin().max(1.0e-6);
    let initial_offset = radius / sin_half;
    let mut f = corner_pos + bisector * initial_offset;

    // Running parameter hints: each iteration's projection
    // starts from the previous iteration's result, so NURBS
    // projections track the fillet centre smoothly instead of
    // jumping to unrelated curve minima.
    let mut hint_a = t_corner_a;
    let mut hint_b = t_corner_b;

    for _ in 0..40 {
        let t_a = project_point_local(curve_a, f, hint_a)?;
        let tp_a = curve_a.point_at(t_a);
        let delta_a = f - tp_a;
        let dist_a = delta_a.length();
        hint_a = Some(t_a);

        let t_b = project_point_local(curve_b, f, hint_b)?;
        let tp_b = curve_b.point_at(t_b);
        let delta_b = f - tp_b;
        let dist_b = delta_b.length();
        hint_b = Some(t_b);

        let res_a = dist_a - radius;
        let res_b = dist_b - radius;
        if res_a.abs() < 1.0e-10 && res_b.abs() < 1.0e-10 {
            // Per-edge wrong-side check: F must lie on the
            // bisector-side perpendicular of each edge.
            let n_a_interior = interior_perp_at_corner(dir_a, bisector, plane_normal);
            let n_b_interior = interior_perp_at_corner(dir_b, bisector, plane_normal);
            let side_a = (f - corner_pos).dot(&n_a_interior);
            let side_b = (f - corner_pos).dot(&n_b_interior);
            if side_a <= 1.0e-9 || side_b <= 1.0e-9 {
                return Err(format!(
                    "generic fillet Newton: converged to wrong-side root (side_a = {side_a:.6}, side_b = {side_b:.6})"
                ));
            }
            return Ok((tp_a, tp_b, f));
        }

        if dist_a < 1.0e-12 || dist_b < 1.0e-12 {
            return Err("generic fillet Newton: fillet centre on support curve".into());
        }
        let n_a = delta_a * (1.0 / dist_a);
        let n_b = delta_b * (1.0 / dist_b);

        let j11 = n_a.dot(&u_basis);
        let j12 = n_a.dot(&v_basis);
        let j21 = n_b.dot(&u_basis);
        let j22 = n_b.dot(&v_basis);
        let det = j11 * j22 - j12 * j21;
        if det.abs() < 1.0e-12 {
            return Err(format!(
                "generic fillet Newton: Jacobian singular (det = {det:.6})"
            ));
        }
        let du = (-res_a * j22 + res_b * j12) / det;
        let dv = (res_a * j21 - res_b * j11) / det;
        f = f + u_basis * du + v_basis * dv;
    }
    Err("generic fillet Newton failed to converge in 40 iterations".into())
}

/// For NURBS curves, try to locate the parameter at which the
/// curve passes through (or closest to) `corner_pos`. This is
/// used as the seed for subsequent local projections. For
/// Line/Arc curves, returns `None` — their `closest_point` is
/// already a single-root solver and doesn't need hinting.
fn corner_parameter_hint(curve: &CurveGeometry, corner_pos: Point3) -> Option<f64> {
    match curve {
        CurveGeometry::Nurbs(_) => curve.closest_point(&corner_pos),
        _ => None,
    }
}

/// Project `query` onto `curve`, preferring a local minimum
/// near `hint` (if supplied) for NURBS curves. For Line/Arc
/// curves the trait-provided `closest_point` is a global
/// minimum but also the only minimum, so we delegate to it.
fn project_point_local(
    curve: &CurveGeometry,
    query: Point3,
    hint: Option<f64>,
) -> Result<f64, String> {
    match curve {
        CurveGeometry::Nurbs(nurbs) => {
            if let Some(t0) = hint {
                // Newton-on-parameter-projection, starting from
                // the hint. Falls back to trait `closest_point`
                // if Newton diverges.
                match local_closest_point(nurbs, query, t0) {
                    Ok(t) => Ok(t),
                    Err(_) => curve
                        .closest_point(&query)
                        .ok_or_else(|| "closest_point returned None".to_string()),
                }
            } else {
                curve
                    .closest_point(&query)
                    .ok_or_else(|| "closest_point returned None".to_string())
            }
        }
        _ => curve
            .closest_point(&query)
            .ok_or_else(|| "closest_point returned None".to_string()),
    }
}

/// Find a parameter `t` near `t_initial` at which the NURBS
/// curve's closest-point equation `(query − C(t)) · C'(t) = 0`
/// holds. Newton iteration with curvature-aware derivative,
/// clamped to the curve's parameter range.
fn local_closest_point(nurbs: &NurbsCurve, query: Point3, t_initial: f64) -> Result<f64, String> {
    let (t_lo, t_hi) = nurbs.parameter_range();
    let mut t = t_initial.clamp(t_lo, t_hi);
    for _ in 0..25 {
        let c = nurbs.point_at(t);
        let dc = nurbs.tangent_at(t);
        let diff = query - c;
        let f = diff.dot(&dc);
        if f.abs() < 1.0e-12 {
            return Ok(t);
        }
        // f'(t) = -|C'(t)|² + (query - C(t)) · C''(t). We
        // approximate C''(t) with a central finite difference
        // so we don't require a second-derivative API from the
        // curve. At the solution the second term vanishes, so
        // near convergence we behave like Newton on -|C'|².
        let h = 1.0e-5 * (t_hi - t_lo).max(1.0e-3);
        let tp = (t + h).min(t_hi);
        let tm = (t - h).max(t_lo);
        let dc_p = nurbs.tangent_at(tp);
        let dc_m = nurbs.tangent_at(tm);
        let ddc = if (tp - tm).abs() > 0.0 {
            (dc_p - dc_m) * (1.0 / (tp - tm))
        } else {
            Vec3::new(0.0, 0.0, 0.0)
        };
        let df = -dc.dot(&dc) + diff.dot(&ddc);
        if df.abs() < 1.0e-12 {
            return Err("local projection Newton: zero derivative".into());
        }
        let step = f / df;
        let t_next = (t - step).clamp(t_lo, t_hi);
        if (t_next - t).abs() < 1.0e-14 {
            return Ok(t_next);
        }
        t = t_next;
    }
    Err("local projection Newton failed to converge".into())
}

/// In-plane perpendicular to `dir`, oriented toward the
/// `bisector` side (positive dot product with `bisector`). Used
/// by [`solve_generic_fillet`]'s wrong-side post-check to
/// define each edge's interior half-plane at the corner.
fn interior_perp_at_corner(dir: Vec3, bisector: Vec3, plane_normal: Vec3) -> Vec3 {
    let perp = plane_normal.cross(&dir);
    if perp.dot(&bisector) >= 0.0 {
        perp
    } else {
        -perp
    }
}

/// Pick a right-handed orthonormal basis `(u, v)` spanning the
/// plane perpendicular to `normal`. Used only by
/// [`solve_generic_fillet`]'s Newton step, so the specific
/// choice of axes doesn't matter — any orthonormal pair in the
/// plane works.
fn plane_basis(normal: Vec3) -> (Vec3, Vec3) {
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = normal
        .cross(&axis)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Baseline: the correct bisector must succeed
    /// and place the fillet centre at (r, r, 0).
    #[test]
    fn solve_generic_fillet_accepts_correct_side_root() {
        let curve_a = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let curve_b = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(0.0, 1.0, 0.0),
        };
        let corner = Point3::new(0.0, 0.0, 0.0);
        let plane_normal = Vec3::new(0.0, 0.0, 1.0);
        let dir_a = Vec3::new(1.0, 0.0, 0.0);
        let dir_b = Vec3::new(0.0, 1.0, 0.0);
        let bisector = Vec3::new(1.0, 1.0, 0.0).normalized().unwrap();
        let r = 0.1_f64;
        let angle = std::f64::consts::FRAC_PI_2;
        let (tp_a, tp_b, fc) = solve_generic_fillet(
            &curve_a,
            &curve_b,
            corner,
            dir_a,
            dir_b,
            bisector,
            plane_normal,
            r,
            angle,
        )
        .unwrap();
        assert!((fc.x - 0.1).abs() < 1e-9);
        assert!((fc.y - 0.1).abs() < 1e-9);
        assert!((tp_a.x - 0.1).abs() < 1e-9);
        assert!(tp_a.y.abs() < 1e-9);
        assert!(tp_b.x.abs() < 1e-9);
        assert!((tp_b.y - 0.1).abs() < 1e-9);
    }
}
