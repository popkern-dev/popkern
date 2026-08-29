//! Cross-tangent control-point fitting for support-owned Gregory boundaries.

use super::BoundarySide;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};

pub(super) fn selected_cross_tangent_cps<const N: usize>(
    boundary: &[Point3; N],
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    side: BoundarySide,
) -> [Vec3; N] {
    let fitted = fitted_cross_tangent_cps(boundary, support, twin_pcurve, side);
    let local = cp_local_cross_tangent_cps(boundary, support, twin_pcurve, side);
    let fitted_score = score_cross_tangent_cps(boundary, support, twin_pcurve, side, &fitted);
    let local_score = score_cross_tangent_cps(boundary, support, twin_pcurve, side, &local);
    if crate::direct_apex_trace::support_frames_enabled() {
        let trace_t = 0.65;
        let trace_point = bezier_point(boundary, trace_t);
        let trace_tangent = bezier_tangent(boundary, trace_t);
        let trace_normal =
            support_normal_at_boundary_t(support, twin_pcurve, &trace_point, trace_t);
        let trace_cross = bezier_vec(&fitted, trace_t);
        eprintln!(
            "[direct-apex support-frame] side={} fitted_score={fitted_score:.6e} local_score={local_score:.6e} t=.65 tangent=({:.3},{:.3},{:.3}) normal={} fitted_cross=({:.3},{:.3},{:.3})",
            side.name(),
            trace_tangent.x,
            trace_tangent.y,
            trace_tangent.z,
            trace_normal
                .map(|n| format!("({:.3},{:.3},{:.3})", n.x, n.y, n.z))
                .unwrap_or_else(|| "none".to_string()),
            trace_cross.x,
            trace_cross.y,
            trace_cross.z
        );
    }
    if fitted_score <= local_score {
        fitted
    } else {
        local
    }
}

fn fitted_cross_tangent_cps<const N: usize>(
    boundary: &[Point3; N],
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    side: BoundarySide,
) -> [Vec3; N] {
    use popkern_core::geometry::fitting::fit_curve;

    let degree = N - 1;
    let n_samples = 40usize.max(N * 4);
    let mut samples = Vec::with_capacity(n_samples);
    let mut params = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1) as f64;
        let point = bezier_point(boundary, t);
        let tangent = bezier_tangent(boundary, t);
        let Some(normal) = support_normal_at_boundary_t(support, twin_pcurve, &point, t) else {
            return cp_local_cross_tangent_cps(boundary, support, twin_pcurve, side);
        };
        let cross = oriented_cross_tangent(tangent, normal, side)
            .normalized()
            .map(|unit| unit * tangent.length())
            .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
        samples.push(Point3::new(cross.x, cross.y, cross.z));
        params.push(t);
    }
    let Ok(fitted) = fit_curve(&samples, &params, degree, N) else {
        return cp_local_cross_tangent_cps(boundary, support, twin_pcurve, side);
    };
    if fitted.control_points.len() != N {
        return cp_local_cross_tangent_cps(boundary, support, twin_pcurve, side);
    }
    let mut out = [Vec3::new(0.0, 0.0, 0.0); N];
    for (i, cp) in fitted.control_points.iter().enumerate() {
        out[i] = Vec3::new(cp.x, cp.y, cp.z);
    }
    out
}

fn score_cross_tangent_cps<const N: usize>(
    boundary: &[Point3; N],
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    side: BoundarySide,
    cross_cps: &[Vec3; N],
) -> f64 {
    let n_samples = 24usize;
    let mut worst = 0.0_f64;
    for i in 0..n_samples {
        let t = 0.025 + 0.95 * (i as f64) / (n_samples - 1) as f64;
        let point = bezier_point(boundary, t);
        let tangent = bezier_tangent(boundary, t);
        let cross = bezier_vec(cross_cps, t);
        let Some(support_normal) = support_normal_at_boundary_t(support, twin_pcurve, &point, t)
        else {
            return f64::INFINITY;
        };
        let cap_normal = match side {
            BoundarySide::UZero | BoundarySide::UOne => cross.cross(&tangent),
            BoundarySide::VZero | BoundarySide::VOne => tangent.cross(&cross),
        };
        let Some(cap_unit) = cap_normal.normalized() else {
            return f64::INFINITY;
        };
        let Some(support_unit) = support_normal.normalized() else {
            return f64::INFINITY;
        };
        worst = worst.max(1.0 - cap_unit.dot(&support_unit).abs());
    }
    worst
}

fn cp_local_cross_tangent_cps<const N: usize>(
    boundary: &[Point3; N],
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    side: BoundarySide,
) -> [Vec3; N] {
    let mut out = [Vec3::new(0.0, 0.0, 0.0); N];
    for (i, item) in out.iter_mut().enumerate() {
        let t = if N > 1 {
            i as f64 / (N - 1) as f64
        } else {
            0.0
        };
        if let Some(cross) =
            boundary_owned_cross_tangent(boundary, i, support, twin_pcurve, t, side)
        {
            *item = cross;
        }
    }
    out
}

fn bezier_vec<const N: usize>(cps: &[Vec3; N], t: f64) -> Vec3 {
    let degree = N - 1;
    let mut out = Vec3::new(0.0, 0.0, 0.0);
    for (i, cp) in cps.iter().enumerate() {
        out = out + *cp * bernstein(degree, i, t);
    }
    out
}

fn bezier_point<const N: usize>(cps: &[Point3; N], t: f64) -> Point3 {
    let degree = N - 1;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for (i, cp) in cps.iter().enumerate() {
        let b = bernstein(degree, i, t);
        x += b * cp.x;
        y += b * cp.y;
        z += b * cp.z;
    }
    Point3::new(x, y, z)
}

fn bezier_tangent<const N: usize>(cps: &[Point3; N], t: f64) -> Vec3 {
    if N < 2 {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    let degree = N - 1;
    let deriv_degree = degree - 1;
    let mut out = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..degree {
        let delta = point_delta(cps[i + 1], cps[i]) * degree as f64;
        out = out + delta * bernstein(deriv_degree, i, t);
    }
    out
}

fn bernstein(degree: usize, i: usize, t: f64) -> f64 {
    binomial(degree, i) * t.powi(i as i32) * (1.0 - t).powi((degree - i) as i32)
}

fn binomial(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    let mut out = 1.0;
    for i in 0..k {
        out *= (n - i) as f64 / (i + 1) as f64;
    }
    out
}

fn boundary_owned_cross_tangent<const N: usize>(
    boundary: &[Point3; N],
    index: usize,
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    t: f64,
    side: BoundarySide,
) -> Option<Vec3> {
    let tangent = boundary_cp_tangent(boundary, index)?;
    let normal = support_normal_at_boundary_t(support, twin_pcurve, &boundary[index], t)?;
    let cross = oriented_cross_tangent(tangent, normal, side);
    cross.normalized().map(|unit| unit * tangent.length())
}

fn oriented_cross_tangent(tangent: Vec3, normal: Vec3, side: BoundarySide) -> Vec3 {
    match side {
        BoundarySide::UZero | BoundarySide::UOne => tangent.cross(&normal),
        BoundarySide::VZero | BoundarySide::VOne => normal.cross(&tangent),
    }
}

fn boundary_cp_tangent<const N: usize>(boundary: &[Point3; N], index: usize) -> Option<Vec3> {
    if N < 2 || index >= N {
        return None;
    }
    let degree = (N - 1) as f64;
    let tangent = if index == 0 {
        point_delta(boundary[1], boundary[0]) * degree
    } else if index + 1 == N {
        point_delta(boundary[N - 1], boundary[N - 2]) * degree
    } else {
        point_delta(boundary[index + 1], boundary[index - 1]) * (0.5 * degree)
    };
    (tangent.length() > 1e-12).then_some(tangent)
}

pub(super) fn support_normal_at_boundary_t(
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    point: &Point3,
    boundary_t: f64,
) -> Option<Vec3> {
    let (t0, t1) = twin_pcurve.parameter_range();
    let forward_t = t0 + boundary_t * (t1 - t0);
    let reverse_t = t1 - boundary_t * (t1 - t0);
    let forward_uv = twin_pcurve.point_at(forward_t);
    let reverse_uv = twin_pcurve.point_at(reverse_t);
    let forward_p = support.point_at(forward_uv.x, forward_uv.y);
    let reverse_p = support.point_at(reverse_uv.x, reverse_uv.y);
    let uv = if point.distance_to(&forward_p) <= point.distance_to(&reverse_p) {
        forward_uv
    } else {
        reverse_uv
    };
    support.normal_at(uv.x, uv.y).normalized()
}

fn point_delta(a: Point3, b: Point3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

pub(super) fn offset_point(point: Point3, direction: Vec3, scale: f64) -> Point3 {
    Point3::new(
        point.x + scale * direction.x,
        point.y + scale * direction.y,
        point.z + scale * direction.z,
    )
}
