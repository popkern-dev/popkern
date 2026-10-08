//! Terminal section of a NON-cylindrical strip on an oblique end plane
//! (PK-740).
//!
//! A cylinder strip ending on a plane perpendicular to its axis ends on
//! an exact circular arc (`hole_rim_terminal_disc`). A torus strip — the
//! root fillet of a fin's cap edge on a cylinder housing — ends on the
//! fin's SIDE plane, which is offset from radial, so the plane does not
//! contain a meridian and the section `strip ∩ plane` is a spiric curve,
//! not a circle. This module builds that section between the strip's two
//! rail ends:
//!
//! - samples march across the strip's chart on the straight uv segment
//!   between the rail ends' charts, each Newton-corrected ONTO the plane
//!   along the chart direction perpendicular to the segment (the plane
//!   is transversal to the strip there — a near-meridian section);
//! - a cubic interpolates the samples, densified until every
//!   between-sample probe sits within `0.1 × linear` of BOTH the plane
//!   and the strip surface;
//! - any non-convergence, a probe outside the band at the densest
//!   level, or an end missing its rail end is a structured decline the
//!   caller turns into a transactional failure.
//!
//! **Tolerance:** the caller's linear tolerance; the fit gate is a tenth
//! of it so the authored edge lies well inside check 9's band.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::fitting::interpolate_curve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

/// Densification ladder (sample counts) for the section fit.
const SAMPLE_LADDER: [usize; 5] = [9, 17, 33, 65, 129];

/// The section `surface ∩ plane(origin, normal)` from `q1` to `q2`, as
/// a cubic NURBS within `0.1 × linear` of both.
pub(super) fn strip_plane_section(
    surface: &SurfaceGeometry,
    origin: Point3,
    normal: Vec3,
    q1: Point3,
    q2: Point3,
    linear: f64,
) -> Result<NurbsCurve, String> {
    let n = normal
        .normalized()
        .ok_or_else(|| "terminal section: degenerate end-plane normal".to_string())?;
    let uv1 = project_point_to_uv(surface, &q1)
        .ok_or_else(|| "terminal section: rail end 1 has no chart point".to_string())?;
    let raw2 = project_point_to_uv(surface, &q2)
        .ok_or_else(|| "terminal section: rail end 2 has no chart point".to_string())?;
    // Same period branch as the first end: the short way round.
    let (pu, pv) = surface.periodicity();
    let near = |x: f64, to: f64, p: Option<f64>| match p {
        Some(p) if p > 0.0 => x + ((to - x) / p).round() * p,
        _ => x,
    };
    let uv2 = (near(raw2.0, uv1.0, pu), near(raw2.1, uv1.1, pv));
    let gate = (linear * 0.1).max(1e-12);
    for (q, uv) in [(q1, uv1), (q2, uv2)] {
        let off = surface.point_at(uv.0, uv.1).distance_to(&q);
        if off > linear * 10.0 {
            return Err(format!(
                "terminal section: rail end ({:.4},{:.4},{:.4}) is {off:.3e} off the strip",
                q.x, q.y, q.z
            ));
        }
    }
    let d = (uv2.0 - uv1.0, uv2.1 - uv1.1);
    let len = (d.0 * d.0 + d.1 * d.1).sqrt();
    if len <= 0.0 {
        return Err("terminal section: rail ends share one chart point".to_string());
    }
    let perp = (-d.1 / len, d.0 / len);
    let signed = |p: Point3| (p - origin).dot(&n);
    // One on-plane point at chart fraction `s` (ends pinned exactly).
    let on_plane = |s: f64| -> Option<Point3> {
        if s <= 0.0 {
            return Some(q1);
        }
        if s >= 1.0 {
            return Some(q2);
        }
        let base = (uv1.0 + d.0 * s, uv1.1 + d.1 * s);
        let mut w = 0.0;
        for _ in 0..40 {
            let (u, v) = (base.0 + perp.0 * w, base.1 + perp.1 * w);
            let f = signed(surface.point_at(u, v));
            if f.abs() <= gate * 1e-3 {
                return Some(surface.point_at(u, v));
            }
            let df = n.dot(&(surface.partial_u(u, v) * perp.0 + surface.partial_v(u, v) * perp.1));
            if df.abs() < 1e-14 {
                return None;
            }
            w -= f / df;
        }
        None
    };
    let mut last_err = String::new();
    for &m in &SAMPLE_LADDER {
        let pts: Option<Vec<Point3>> = (0..m)
            .map(|k| on_plane(k as f64 / (m - 1) as f64))
            .collect();
        let Some(pts) = pts else {
            return Err("terminal section: on-plane Newton did not converge".to_string());
        };
        // Chord-length parameters on [0, 1].
        let mut params = vec![0.0];
        for w in pts.windows(2) {
            let last = params.last().copied().unwrap_or(0.0);
            params.push(last + w[0].distance_to(&w[1]));
        }
        let total = params.last().copied().unwrap_or(0.0);
        if total <= 0.0 {
            return Err("terminal section: degenerate section".to_string());
        }
        for p in &mut params {
            *p /= total;
        }
        let curve = interpolate_curve(&pts, &params, 3)
            .map_err(|e| format!("terminal section: interpolation failed: {e}"))?;
        // Probe between samples against the plane and the surface.
        let mut worst = 0.0_f64;
        for k in 0..m - 1 {
            let t = 0.5 * (params[k] + params[k + 1]);
            let p = curve.point_at(t);
            let plane_off = signed(p).abs();
            let surf_off = project_point_to_uv(surface, &p)
                .map(|(u, v)| surface.point_at(u, v).distance_to(&p))
                .unwrap_or(f64::INFINITY);
            worst = worst.max(plane_off).max(surf_off);
        }
        if worst <= gate {
            return Ok(curve);
        }
        last_err = format!("terminal section: fit {worst:.3e} off at {m} samples");
    }
    Err(last_err)
}
