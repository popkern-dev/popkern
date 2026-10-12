//! Station scanning and dense level-curve tracing for the radial
//! level-curve fast path (PK-461) — split from
//! `radial_level_curve.rs` under the 500-line file budget. The
//! mechanism rationale lives in the parent module's header.

use std::f64::consts::{PI, TAU};

use super::{AxisymmetricFrame, StationDir, StationSample, trace_decline, wall_uv};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Scan one station line: bracket cells where `ρ − R` changes sign.
/// `None` when the graze gate trips (a near-band sample not adjacent
/// to a bracket — tangency/coincidence structure).
pub(super) fn station_scan(
    wall: &SurfaceGeometry,
    dir: StationDir,
    s: f64,
    r_range: (f64, f64),
    n_r: usize,
    frame: &AxisymmetricFrame,
    graze: f64,
) -> Option<Vec<usize>> {
    let mut rho = Vec::with_capacity(n_r);
    for i in 0..n_r {
        let r = r_range.0 + (r_range.1 - r_range.0) * i as f64 / (n_r - 1) as f64;
        let (u, v) = wall_uv(dir, s, r);
        rho.push(frame.level(&wall.point_at(u, v)));
    }
    // Sign classes fold an exact-zero sample into the negative side
    // so a crossing that lands exactly on a grid sample still yields
    // one bracket (a flat-at-radius stretch still declines through
    // the graze gate below).
    let sign = |d: f64| d > 0.0;
    let brackets: Vec<usize> = (0..n_r - 1)
        .filter(|&i| sign(rho[i]) != sign(rho[i + 1]))
        .collect();
    for (i, d) in rho.iter().enumerate() {
        if d.abs() <= graze {
            let excused = brackets.contains(&i) || (i > 0 && brackets.contains(&(i - 1)));
            if !excused {
                return None;
            }
        }
    }
    Some(brackets)
}

/// Bisect `ρ = R` inside a bracket cell. The cell has a strict sign
/// change, so bisection converges unconditionally.
pub(super) fn bisect_radial(
    wall: &SurfaceGeometry,
    dir: StationDir,
    s: f64,
    mut lo: f64,
    mut hi: f64,
    frame: &AxisymmetricFrame,
) -> f64 {
    let eval = |r: f64| {
        let (u, v) = wall_uv(dir, s, r);
        frame.level(&wall.point_at(u, v))
    };
    let f_lo = eval(lo);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if (eval(mid) > 0.0) == (f_lo > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Unwrapped chart excursion of a trace, in periods.
pub(super) fn trace_wraps(t: &[StationSample]) -> f64 {
    match (t.first(), t.last()) {
        (Some(a), Some(b)) => (b.theta - a.theta).abs() / TAU,
        _ => 0.0,
    }
}

pub(super) enum DensePass {
    Ok(Vec<Vec<StationSample>>),
    Aliased,
    Decline,
}

/// Root-find every crossing at `n` stations, verifying the crossing
/// structure at each and unwrapping the chart angle continuously.
#[allow(clippy::too_many_arguments)]
pub(super) fn dense_pass(
    wall: &SurfaceGeometry,
    dir: StationDir,
    s_range: (f64, f64),
    r_range: (f64, f64),
    n_r: usize,
    n_dense: usize,
    k: usize,
    frame: &AxisymmetricFrame,
    graze: f64,
) -> DensePass {
    let mut traces: Vec<Vec<StationSample>> =
        (0..k).map(|_| Vec::with_capacity(n_dense + 1)).collect();
    let mut prev_brackets: Option<Vec<usize>> = None;
    let cell = (r_range.1 - r_range.0) / (n_r - 1) as f64;
    for j in 0..=n_dense {
        let s = s_range.0 + (s_range.1 - s_range.0) * j as f64 / n_dense as f64;
        let Some(brackets) = station_scan(wall, dir, s, r_range, n_r, frame, graze) else {
            trace_decline("dense graze trip", k, 0.0);
            return DensePass::Decline;
        };
        if brackets.len() != k {
            trace_decline("dense crossing-count change", k, 0.0);
            return DensePass::Decline;
        }
        if let Some(prev) = &prev_brackets {
            for c in 0..k {
                if prev[c].abs_diff(brackets[c]) > 1 {
                    trace_decline("dense bracket drift", k, 0.0);
                    return DensePass::Decline;
                }
            }
        }
        for (c, &cell_idx) in brackets.iter().enumerate() {
            let lo = r_range.0 + cell * cell_idx as f64;
            let r = bisect_radial(wall, dir, s, lo, lo + cell, frame);
            let (u, v) = wall_uv(dir, s, r);
            let p = wall.point_at(u, v);
            let (theta_p, v_cyl) = frame.theta_v(&p);
            let theta = match traces[c].last() {
                None => theta_p,
                Some(prev) => {
                    let mut d = theta_p - branch_angle(prev.theta);
                    while d > PI {
                        d -= TAU;
                    }
                    while d < -PI {
                        d += TAU;
                    }
                    if d.abs() >= PI / 2.0 {
                        return DensePass::Aliased;
                    }
                    prev.theta + d
                }
            };
            traces[c].push(StationSample {
                u_wall: u,
                v_wall: v,
                point: p,
                theta,
                v_cyl,
            });
        }
        prev_brackets = Some(brackets);
    }
    DensePass::Ok(traces)
}

/// Principal-branch image of an unwrapped angle, in (-π, π].
pub(super) fn branch_angle(theta: f64) -> f64 {
    let mut t = theta.rem_euclid(TAU);
    if t > PI {
        t -= TAU;
    }
    t
}
