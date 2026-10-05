//! Window clipping, seam splitting, and record emission for the
//! radial level-curve fast path (PK-461) — split from
//! `radial_level_curve.rs` under the 500-line file budget. The
//! mechanism rationale lives in the parent module's header.

use popkern_math::PortableMath;
use std::f64::consts::{PI, TAU};

use super::radial_level_curve_trace::{bisect_radial, branch_angle};
use super::{CylinderFrame, StationDir, StationSample, wall_uv};
use crate::geometry::fitting::fit_curve;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// Clip a trace to the cylinder's queried axial window, then split
/// it at exact chart-seam crossings (θ ≡ 0 mod 2π). Each returned
/// run carries `(sample, principal-branch u2)` with the seam
/// crossing emitted on BOTH branches — the tail sample exactly at
/// the period edge, the head sample exactly at 0 — so every record's
/// pcurve stays on one chart branch and terminates on the seam.
#[allow(clippy::too_many_arguments)]
pub(super) fn split_trace(
    wall: &SurfaceGeometry,
    dir: StationDir,
    trace: &[StationSample],
    r_range: (f64, f64),
    frame: &CylinderFrame,
    cyl_v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<Vec<(StationSample, f64)>> {
    // Axial window clip (the window is already padded by the caller;
    // pad again like the cross-drill arm so a boundary-hugging trace
    // is never dropped by sampling luck).
    let pad = (cyl_v_range.1 - cyl_v_range.0).abs() * 0.5 + tolerance.linear();
    let in_window = |s: &StationSample| {
        s.v_cyl > cyl_v_range.0.min(cyl_v_range.1) - pad
            && s.v_cyl < cyl_v_range.0.max(cyl_v_range.1) + pad
    };
    let mut segments: Vec<Vec<&StationSample>> = Vec::new();
    let mut current: Vec<&StationSample> = Vec::new();
    for s in trace {
        if in_window(s) {
            current.push(s);
        } else if !current.is_empty() {
            segments.push(std::mem::take(&mut current));
        }
    }
    if !current.is_empty() {
        segments.push(current);
    }

    let mut runs: Vec<Vec<(StationSample, f64)>> = Vec::new();
    for seg in segments {
        if seg.len() < 2 {
            continue;
        }
        let mut run: Vec<(StationSample, f64)> = Vec::new();
        let wrap_of = |theta: f64| (theta / TAU).floor();
        let mut wrap = wrap_of(seg[0].theta);
        run.push((seg[0].clone(), seg[0].theta - wrap * TAU));
        for w in seg.windows(2) {
            let (a, b) = (w[0], w[1]);
            let wrap_b = wrap_of(b.theta);
            if wrap_b != wrap {
                // Exactly one seam between a and b (|Δθ| < π/2). Find
                // the crossing station and emit it on both branches.
                let target = if b.theta > a.theta {
                    (wrap + 1.0) * TAU
                } else {
                    wrap * TAU
                };
                if let Some(cross) = solve_seam_crossing(wall, dir, a, b, target, r_range, frame) {
                    let tail_u2 = target - wrap * TAU; // exactly 0 or TAU
                    let head_u2 = target - wrap_b * TAU;
                    // Drop a neighbouring regular station that sits on
                    // top of the crossing (degenerate chord).
                    let min_gap = 0.05 * a.point.distance_to(&b.point);
                    if run
                        .last()
                        .is_some_and(|(s, _)| s.point.distance_to(&cross.point) < min_gap)
                    {
                        run.pop();
                    }
                    run.push((cross.clone(), tail_u2));
                    runs.push(std::mem::take(&mut run));
                    run.push((cross.clone(), head_u2));
                    if b.point.distance_to(&cross.point) >= min_gap {
                        run.push((b.clone(), b.theta - wrap_b * TAU));
                    }
                } else {
                    // Crossing solve failed: decline the seam cut and
                    // keep the pair honest by splitting at the raw
                    // samples (the fitter's chord filter drops blips;
                    // the honesty gate rejects anything off-band).
                    runs.push(std::mem::take(&mut run));
                    run.push((b.clone(), b.theta - wrap_b * TAU));
                }
                wrap = wrap_b;
            } else {
                run.push((b.clone(), b.theta - wrap * TAU));
            }
        }
        if !run.is_empty() {
            runs.push(run);
        }
    }
    runs
}

/// Solve the station where the trace's unwrapped angle equals
/// `target` between two bracketing samples, by bisection on the
/// station coordinate (each probe re-roots the radial crossing).
pub(super) fn solve_seam_crossing(
    wall: &SurfaceGeometry,
    dir: StationDir,
    a: &StationSample,
    b: &StationSample,
    target: f64,
    r_range: (f64, f64),
    frame: &CylinderFrame,
) -> Option<StationSample> {
    let (s_a, r_a) = match dir {
        StationDir::AlongV => (a.v_wall, a.u_wall),
        StationDir::AlongU => (a.u_wall, a.v_wall),
    };
    let (s_b, r_b) = match dir {
        StationDir::AlongV => (b.v_wall, b.u_wall),
        StationDir::AlongU => (b.u_wall, b.v_wall),
    };
    // A crossing landing exactly on a sample (a trace starting on
    // the seam) needs no solve — bisecting a zero-endpoint interval
    // walks to the wrong end and mislabels the crossing angle.
    for sample in [a, b] {
        if (sample.theta - target).abs() < 1e-9 {
            let mut exact = sample.clone();
            exact.theta = target;
            return Some(exact);
        }
    }
    let span = (r_range.1 - r_range.0).abs().max(1e-30);
    let probe = |s: f64| -> Option<(f64, StationSample)> {
        // Radial bracket: widen around the interpolated root until a
        // sign change appears (the root moves slowly between two
        // adjacent dense stations; a constant level coordinate needs
        // the widening, so the seed width never collapses to zero).
        let t = (s - s_a) / (s_b - s_a);
        let r_c = r_a + (r_b - r_a) * t;
        let mut width = (2.0 * (r_a - r_b).abs()).max(1e-4 * span);
        let eval = |r: f64| {
            let (u, v) = wall_uv(dir, s, r);
            frame.radial(&wall.point_at(u, v)) - frame.radius
        };
        let mut bracket = None;
        for _ in 0..16 {
            let lo = (r_c - width).max(r_range.0);
            let hi = (r_c + width).min(r_range.1);
            let sign = |d: f64| d > 0.0;
            if sign(eval(lo)) != sign(eval(hi)) {
                bracket = Some((lo, hi));
                break;
            }
            width *= 2.0;
        }
        let (lo, hi) = bracket?;
        let r = bisect_radial(wall, dir, s, lo, hi, frame);
        let (u, v) = wall_uv(dir, s, r);
        let p = wall.point_at(u, v);
        let (theta_p, v_cyl) = frame.theta_v(&p);
        // Local unwrap against a: |θ − θ_a| < π/2 on this interval.
        let mut d = theta_p - branch_angle(a.theta);
        while d > PI {
            d -= TAU;
        }
        while d < -PI {
            d += TAU;
        }
        let theta = a.theta + d;
        Some((
            theta - target,
            StationSample {
                u_wall: u,
                v_wall: v,
                point: p,
                theta,
                v_cyl,
            },
        ))
    };
    let (mut lo, mut hi) = (s_a, s_b);
    let (f_lo, _) = probe(lo)?;
    let mut best: Option<StationSample> = None;
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let (f_mid, sample) = probe(mid)?;
        if (f_mid > 0.0) == (f_lo > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
        best = Some(sample);
    }
    let mut sample = best?;
    sample.theta = target;
    Some(sample)
}

/// Build a `MarchSample` with the wall/cylinder charts on the right
/// operand slots. `theta_branch` is the principal-branch cylinder
/// angle for the sample's run.
pub(super) fn make_sample(wall_first: bool, s: &StationSample, theta_branch: f64) -> MarchSample {
    if wall_first {
        MarchSample {
            u1: s.u_wall,
            v1: s.v_wall,
            u2: theta_branch,
            v2: s.v_cyl,
            point: s.point,
        }
    } else {
        MarchSample {
            u1: theta_branch,
            v1: s.v_cyl,
            u2: s.u_wall,
            v2: s.v_wall,
            point: s.point,
        }
    }
}

/// Verify a fitted record against both surfaces through its pcurves.
pub(super) fn record_on_surfaces(
    record: &crate::geometry::curve_on_surface::CurveOnSurface,
    wall_first: bool,
    wall: &SurfaceGeometry,
    frame: &CylinderFrame,
    surface_band: f64,
    drift_band: f64,
) -> bool {
    use crate::geometry::curve::CurveEval;
    let (t0, t1) = record.curve_3d.parameter_range();
    const CHECKS: usize = 65;
    let mut worst_radial = 0.0_f64;
    let mut worst_wall_drift = 0.0_f64;
    let mut worst_cyl_drift = 0.0_f64;
    for i in 0..=CHECKS {
        let t = t0 + (t1 - t0) * i as f64 / CHECKS as f64;
        let p = record.curve_3d.point_at(t);
        // On-cylinder: exact analytic distance to the surface.
        worst_radial = worst_radial.max((frame.radial(&p) - frame.radius).abs());
        let (wall_uv_curve, cyl_uv_curve) = if wall_first {
            (&record.curve_uv1, &record.curve_uv2)
        } else {
            (&record.curve_uv2, &record.curve_uv1)
        };
        // Same-parameter pcurve consistency (parameterization drift).
        let wuv = wall_uv_curve.point_at(t);
        worst_wall_drift = worst_wall_drift.max(wall.point_at(wuv.x, wuv.y).distance_to(&p));
        let cuv = cyl_uv_curve.point_at(t);
        let cyl_point = frame.origin
            + frame.e1 * (frame.radius * cuv.x.pcos())
            + frame.e2 * (frame.radius * cuv.x.psin())
            + frame.axis_raw * cuv.y;
        worst_cyl_drift = worst_cyl_drift.max(cyl_point.distance_to(&p));
    }
    let ok = worst_radial <= surface_band
        && worst_wall_drift <= drift_band
        && worst_cyl_drift <= drift_band;
    if !ok && crate::trace::trace_enabled("radial-level") {
        eprintln!(
            "  [radial-level] honesty worst: radial {worst_radial:.3e} (band {surface_band:.3e}) \
             wall drift {worst_wall_drift:.3e} cyl drift {worst_cyl_drift:.3e} (band {drift_band:.3e})"
        );
    }
    ok
}

/// Refit both pcurves through EVERY exact sample (full
/// interpolation, the same chord-length parameterization the shared
/// fitter used). The default least-squares control ratio leaves a
/// ~1e-6 parametric truncation that a multi-turn wall chart
/// amplifies through its v-derivative into visible same-parameter
/// drift; interpolation removes the truncation while the exact
/// samples keep the pcurves on-chart.
pub(super) fn interpolate_pcurves(
    samples: &[MarchSample],
    record: &mut crate::geometry::curve_on_surface::CurveOnSurface,
) {
    let m = samples.len();
    if m < 2 {
        return;
    }
    let mut params = vec![0.0];
    let mut total = 0.0;
    for w in samples.windows(2) {
        total += w[0].point.distance_to(&w[1].point);
        params.push(total);
    }
    if total < 1e-30 {
        return;
    }
    for p in &mut params {
        *p /= total;
    }
    let degree = 3.min(m - 1);
    let uv1: Vec<crate::geometry::Point3> = samples
        .iter()
        .map(|s| crate::geometry::Point3::new(s.u1, s.v1, 0.0))
        .collect();
    let uv2: Vec<crate::geometry::Point3> = samples
        .iter()
        .map(|s| crate::geometry::Point3::new(s.u2, s.v2, 0.0))
        .collect();
    if let Ok(c) = fit_curve(&uv1, &params, degree, m) {
        record.curve_uv1 = c;
    }
    if let Ok(c) = fit_curve(&uv2, &params, degree, m) {
        record.curve_uv2 = c;
    }
}
