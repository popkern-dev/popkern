//! Normal-curvature sampling along an ambient direction for the
//! KP continuity gates — split from `kp_gate_helpers.rs` for the
//! 500-line structural cap.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::GregoryPatch;

fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

pub fn normal_curvature_along(gp: &GregoryPatch, u: f64, v: f64, t_ambient: Vec3) -> Option<f64> {
    let su = gp.partial_u(u, v);
    let sv = gp.partial_v(u, v);
    let n = su.cross(&sv).normalized()?;
    let e = dot(su, su);
    let f = dot(su, sv);
    let g = dot(sv, sv);
    let det = e * g - f * f;
    if det.abs() < 1e-14 {
        return None;
    }
    let su_dot_t = dot(su, t_ambient);
    let sv_dot_t = dot(sv, t_ambient);
    let a = (g * su_dot_t - f * sv_dot_t) / det;
    let b = (e * sv_dot_t - f * su_dot_t) / det;
    let norm_sq = e * a * a + 2.0 * f * a * b + g * b * b;
    if norm_sq < 1e-14 {
        return None;
    }
    let a = a / norm_sq.sqrt();
    let b = b / norm_sq.sqrt();
    let (suu, suv, svv) = match gp {
        GregoryPatch::RectangularBiNonic(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        GregoryPatch::RectangularTwoStripe(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        GregoryPatch::RectangularBiseptic(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        GregoryPatch::RectangularBiquinticFull(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        _ => return None,
    };
    let l = dot(n, suu);
    let m = dot(n, suv);
    let nn = dot(n, svv);
    Some(l * a * a + 2.0 * m * a * b + nn * b * b)
}
