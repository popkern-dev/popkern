//! Newton refinement helpers for minimum-distance computations.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;

/// Perform Newton refinement on a curve parameter so that
/// `(C(t) - fixed) · C'(t) = 0`, where `fixed = [x, y, z]`.
///
/// Returns the refined `t`, clamped to `[t_min, t_max]`.
pub(crate) fn refine_curve_param(
    curve: &impl CurveEval,
    mut t: f64,
    t_min: f64,
    t_max: f64,
    fixed: [f64; 3],
) -> f64 {
    for _ in 0..10 {
        let cp = curve.point_at(t);
        let ct = curve.tangent_at(t);
        let dx = cp.x - fixed[0];
        let dy = cp.y - fixed[1];
        let dz = cp.z - fixed[2];
        let dot = dx * ct.x + dy * ct.y + dz * ct.z;
        let ct_sq = ct.dot(&ct);
        if ct_sq <= f64::EPSILON {
            break;
        }
        let t_new = (t - dot / ct_sq).clamp(t_min, t_max);
        if (t_new - t).abs() < 1e-12 {
            t = t_new;
            break;
        }
        t = t_new;
    }
    t
}

/// Perform Newton refinement on surface parameters so that
/// `(S(u,v) - fixed) · Su = 0` and `(S(u,v) - fixed) · Sv = 0`,
/// where `fixed = [x, y, z]` and `bounds = [u_min, u_max, v_min, v_max]`.
///
/// Returns the refined `(u, v)`, each clamped to their respective bounds.
pub(crate) fn refine_surface_params(
    surface: &impl SurfaceEval,
    mut u: f64,
    mut v: f64,
    bounds: [f64; 4],
    fixed: [f64; 3],
) -> (f64, f64) {
    let [u_min, u_max, v_min, v_max] = bounds;
    for _ in 0..10 {
        let sp = surface.point_at(u, v);
        let su = surface.partial_u(u, v);
        let sv = surface.partial_v(u, v);
        let dx = sp.x - fixed[0];
        let dy = sp.y - fixed[1];
        let dz = sp.z - fixed[2];
        let j00 = su.dot(&su);
        let j01 = su.dot(&sv);
        let j11 = sv.dot(&sv);
        let r0 = dx * su.x + dy * su.y + dz * su.z;
        let r1 = dx * sv.x + dy * sv.y + dz * sv.z;
        let det = j00 * j11 - j01 * j01;
        if det.abs() < f64::EPSILON {
            break;
        }
        let du = (j11 * r0 - j01 * r1) / det;
        let dv = (j00 * r1 - j01 * r0) / det;
        let u_new = (u - du).clamp(u_min, u_max);
        let v_new = (v - dv).clamp(v_min, v_max);
        if (u_new - u).abs() < 1e-12 && (v_new - v).abs() < 1e-12 {
            u = u_new;
            v = v_new;
            break;
        }
        u = u_new;
        v = v_new;
    }
    (u, v)
}
