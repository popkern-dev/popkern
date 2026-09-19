//! Radial level-curve fast path: NURBS band × cylinder (PK-461).
//!
//! The thread-groove class — a helical-sweep flank wall crossing its
//! shank cylinder — marched at sag-quality steps (~2·√(tol·r)): an M6
//! cut cost 10–30 s and the 200k-step march budget bounded one cut to
//! ~35 turns. On such a pair the radial field
//! `ρ(u,v) = distance(S(u,v), cylinder axis)` classifies the whole
//! intersection: wherever the level set `ρ = R` is a transversal
//! *graph* over one chart direction, the complete answer is one curve
//! per crossing, traced by per-station 1D root-finds on the built
//! NURBS — no seeding, no marching, no phantom-basin rebranching, and
//! no step budget. Multi-wrap traces split exactly at the cylinder
//! seam per wrap (each record's pcurve stays on the principal chart
//! branch, terminating exactly on the seam — the same contract the
//! marched-trace splitter refines toward).
//!
//! Scope: the arm fires only when the trace wraps the cylinder chart
//! more than 1.5 periods. Sub-wrap pairs (blend strips, graze
//! classes, spring-grind sections) keep the long-tuned marching path;
//! this arm owns exactly the class marching cannot afford. Every gate
//! failure declines to marching — never a partial serve.
//! `POPKERN_ANALYTIC_HELICAL_SSI=0` pins the marching path for A/B
//! measurement and containment (PK-239 precedent; per-call read).
//!
//! Design: docs/decisions/pk461-radial-level-curve-ssi.md.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, fit_march_samples,
};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[path = "radial_level_curve_emit.rs"]
mod radial_level_curve_emit;
#[path = "radial_level_curve_trace.rs"]
mod radial_level_curve_trace;
use radial_level_curve_emit::{make_sample, record_on_surfaces, split_trace};
use radial_level_curve_trace::{DensePass, dense_pass, station_scan, trace_wraps};

/// Graze band as a multiple of linear tolerance: a non-bracket sample
/// this close to the level radius means tangency/coincidence
/// structure the transversal-graph model cannot claim (matches the
/// tier-B coincidence rejection threshold).
const GRAZE_BAND_FACTOR: f64 = 10.0;
/// Stations in the structure-detection grid (and radial scan floor).
const COARSE_STATIONS: usize = 33;
/// Dense stations per chart wrap for the emitted trace: cubic
/// interpolation error at this density is orders below the honesty
/// band on any radius the graze gate admits.
const SAMPLES_PER_WRAP: usize = 256;
/// Hard station cap; a trace needing more declines (traced).
const MAX_STATIONS: usize = 500_000;
/// Minimum unwrapped chart excursion, in periods, for the arm to
/// claim the pair (the multi-wrap scope).
const MIN_WRAPS: f64 = 1.5;
/// Serve containment (PK-494): the downstream boolean/mesh substrate
/// is verified through 40-wrap serves (M6 rod corpus probed per turn
/// count); beyond it boundary-first meshing fails Delaunay coverage
/// on the ultra-anisotropic wall charts and the served mesh TEARS —
/// worse than a decline. Declining here sends the pair to the
/// marching solver, whose step budget truncates and the merge's
/// strand gate declines honestly — exactly the pre-PK-461 behavior
/// for that range. PK-494 lifts this with the anisotropy-aware mesh
/// substrate.
const MAX_SERVED_WRAPS: f64 = 40.5;

/// Cylinder frame: unit axis with the same chart basis `point_at`
/// uses, plus the raw axis for the `v = rel·axis/|axis|²` inversion.
struct CylinderFrame {
    origin: Point3,
    axis_raw: Vec3,
    axis: Vec3,
    e1: Vec3,
    e2: Vec3,
    radius: f64,
}

impl CylinderFrame {
    fn radial(&self, p: &Point3) -> f64 {
        let rel = *p - self.origin;
        (rel - self.axis * rel.dot(&self.axis)).length()
    }
    /// Principal-branch angle in (-π, π] and axial chart coordinate.
    fn theta_v(&self, p: &Point3) -> (f64, f64) {
        let rel = *p - self.origin;
        let theta = rel.dot(&self.e2).atan2(rel.dot(&self.e1));
        let v = rel.dot(&self.axis_raw) / self.axis_raw.length_squared();
        (theta, v)
    }
}

/// Which wall chart direction the trace runs along (the other is the
/// radial direction the level set crosses).
#[derive(Clone, Copy, PartialEq, Eq)]
enum StationDir {
    /// Stations along v; radial scans along u (thread flanks: u is
    /// the profile, v the helical path).
    AlongV,
    /// Stations along u; radial scans along v.
    AlongU,
}

/// One converged trace point with both charts attached.
#[derive(Clone)]
struct StationSample {
    u_wall: f64,
    v_wall: f64,
    point: Point3,
    /// Continuously unwrapped chart angle on the cylinder.
    theta: f64,
    v_cyl: f64,
}

/// Try the radial level-curve path. `None` falls through to the
/// marching solver; `Some` is the complete answer for the pair.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_radial_level_curve(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (wall_first, wall, wall_u_range, wall_v_range, cyl, cyl_v_range) =
        match (surface1, surface2) {
            (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Cylinder { .. }) => {
                (true, surface1, u1_range, v1_range, surface2, v2_range)
            }
            (SurfaceGeometry::Cylinder { .. }, SurfaceGeometry::Nurbs(_)) => {
                (false, surface2, u2_range, v2_range, surface1, v1_range)
            }
            _ => return None,
        };
    let arm_enabled = std::env::var_os("POPKERN_ANALYTIC_HELICAL_SSI")
        .map(|v| v != "0")
        .unwrap_or(true);
    if !arm_enabled {
        return None;
    }
    // Geometrically closed wall charts (coiled spring walls, closed
    // loft sections) keep their tuned seam machinery (PK-443/444):
    // their level sets can be chart-wrapping loops, not graphs.
    {
        let (pu, pv) = crate::geometry::chart_periods::geometric_chart_periods(wall);
        if pu.is_some() || pv.is_some() {
            return None;
        }
    }
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = cyl
    else {
        return None;
    };
    let axis_unit = axis.normalized()?;
    let (e1, e2) = orthonormal_basis(axis_unit);
    let frame = CylinderFrame {
        origin: *origin,
        axis_raw: *axis,
        axis: axis_unit,
        e1,
        e2,
        radius: *radius,
    };
    let graze = GRAZE_BAND_FACTOR * tolerance.linear();

    // Control-hull prescreen: a rational surface lies in its control
    // hull, and distance-to-axis is convex, so the surface cannot
    // reach a radius beyond its farthest control point.
    let SurfaceGeometry::Nurbs(ns) = wall else {
        return None;
    };
    let hull_max = ns
        .control_points
        .iter()
        .map(|p| frame.radial(p))
        .fold(0.0_f64, f64::max);
    if frame.radius > hull_max + graze {
        return None;
    }

    // Clamp the queried wall window to the knot domain (never
    // extrapolate a NURBS chart).
    let (ku, kv) = wall.parameter_range();
    let clamp = |r: (f64, f64), k: (f64, f64)| -> Option<(f64, f64)> {
        let lo = r.0.max(k.0);
        let hi = r.1.min(k.1);
        (hi > lo).then_some((lo, hi))
    };
    let wall_u = clamp(wall_u_range, ku)?;
    let wall_v = clamp(wall_v_range, kv)?;

    for dir in [StationDir::AlongV, StationDir::AlongU] {
        if let Some(result) = trace_direction(
            wall_first,
            wall,
            wall_u,
            wall_v,
            dir,
            &frame,
            cyl_v_range,
            surface1,
            surface2,
            tolerance,
            graze,
        ) {
            return Some(result);
        }
    }
    None
}

/// Station/radial ranges for a direction.
fn dir_ranges(dir: StationDir, wall_u: (f64, f64), wall_v: (f64, f64)) -> ((f64, f64), (f64, f64)) {
    match dir {
        StationDir::AlongV => (wall_v, wall_u),
        StationDir::AlongU => (wall_u, wall_v),
    }
}

/// Map (station, radial) to the wall's (u, v).
fn wall_uv(dir: StationDir, s: f64, r: f64) -> (f64, f64) {
    match dir {
        StationDir::AlongV => (r, s),
        StationDir::AlongU => (s, r),
    }
}

/// Trace all level curves along one station direction. `None`
/// declines the direction (structure, graze, wrap scope, or honesty).
#[allow(clippy::too_many_arguments)]
fn trace_direction(
    wall_first: bool,
    wall: &SurfaceGeometry,
    wall_u: (f64, f64),
    wall_v: (f64, f64),
    dir: StationDir,
    frame: &CylinderFrame,
    cyl_v_range: (f64, f64),
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    tolerance: &Tolerance,
    graze: f64,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (s_range, r_range) = dir_ranges(dir, wall_u, wall_v);
    let SurfaceGeometry::Nurbs(ns) = wall else {
        return None;
    };
    // Radial scan density: resolve every knot span.
    let radial_knots = match dir {
        StationDir::AlongV => &ns.knots_u,
        StationDir::AlongU => &ns.knots_v,
    };
    let spans = radial_knots
        .windows(2)
        .filter(|w| w[1] > w[0] && w[0] >= r_range.0 - 1e-12 && w[1] <= r_range.1 + 1e-12)
        .count()
        .max(1);
    let n_r = (4 * spans + 1).clamp(COARSE_STATIONS, 129);

    // Coarse structure check: identical crossing count at every
    // station, bracket drift ≤ 1 cell between neighbours.
    let station_at = |j: usize, n: usize| -> f64 {
        s_range.0 + (s_range.1 - s_range.0) * j as f64 / (n - 1) as f64
    };
    let mut coarse: Vec<Vec<usize>> = Vec::with_capacity(COARSE_STATIONS);
    for j in 0..COARSE_STATIONS {
        let b = station_scan(
            wall,
            dir,
            station_at(j, COARSE_STATIONS),
            r_range,
            n_r,
            frame,
            graze,
        )?;
        coarse.push(b);
    }
    let k = coarse[0].len();
    if k == 0 {
        return None;
    }
    if coarse.iter().any(|b| b.len() != k) {
        trace_decline("coarse crossing-count mismatch", k, 0.0);
        return None;
    }
    for w in coarse.windows(2) {
        for (a, b) in w[0].iter().zip(&w[1]) {
            if a.abs_diff(*b) > 1 {
                trace_decline("coarse bracket drift", k, 0.0);
                return None;
            }
        }
    }

    // Dense pass, re-run at higher density until the chart angle is
    // unambiguous between neighbours (|Δθ| < π/2) — an unwrap through
    // aliased stations would silently lose whole turns.
    let mut n_dense = (8 * COARSE_STATIONS).max(SAMPLES_PER_WRAP / 4);
    let traces = loop {
        match dense_pass(wall, dir, s_range, r_range, n_r, n_dense, k, frame, graze) {
            DensePass::Ok(traces) => {
                let wraps = traces
                    .iter()
                    .map(|t| trace_wraps(t))
                    .fold(0.0_f64, f64::max);
                let want =
                    ((wraps * SAMPLES_PER_WRAP as f64).ceil() as usize).max(SAMPLES_PER_WRAP / 4);
                if want > MAX_STATIONS {
                    trace_decline("station cap", k, wraps);
                    return None;
                }
                if want as f64 > n_dense as f64 * 1.1 {
                    n_dense = want;
                    continue;
                }
                break traces;
            }
            DensePass::Aliased => {
                n_dense = n_dense.checked_mul(4)?;
                if n_dense > MAX_STATIONS {
                    trace_decline("aliasing cap", k, 0.0);
                    return None;
                }
            }
            DensePass::Decline => return None,
        }
    };

    // Multi-wrap scope: EVERY level curve must wrap beyond the gate,
    // or the whole pair stays on the marching path (mixed emission is
    // not a thing — the arm answers for the pair or not at all).
    for t in &traces {
        if trace_wraps(t) <= MIN_WRAPS {
            trace_decline("sub-wrap scope", k, trace_wraps(t));
            return None;
        }
        if trace_wraps(t) > MAX_SERVED_WRAPS {
            trace_decline("serve containment (PK-494)", k, trace_wraps(t));
            return None;
        }
    }

    // Emit: window-clip, seam-split, crease-split, fit, and verify
    // every record.
    let mut records = Vec::new();
    for t in &traces {
        let runs = split_trace(wall, dir, t, r_range, frame, cyl_v_range, tolerance);
        let runs: Vec<Vec<MarchSample>> = runs
            .into_iter()
            .map(|run| {
                run.iter()
                    .map(|(sample, theta_branch)| make_sample(wall_first, sample, *theta_branch))
                    .collect()
            })
            .collect();
        // C0 crease crossings cut the runs too (PK-442 contract —
        // a smooth fit across a faceted chart's corner rounds the
        // true tangent break). No-op for smooth walls.
        let runs =
            crate::geometry::intersection::surface_surface_crease_split::split_samples_at_c0_creases(
                runs, surface1, surface2, tolerance,
            );
        for samples in runs {
            let chord: f64 = samples
                .windows(2)
                .map(|w| w[0].point.distance_to(&w[1].point))
                .sum();
            if samples.len() < 2 || chord <= tolerance.linear() {
                continue;
            }
            match fit_march_samples(&samples, false) {
                Ok(curve) => records.push((samples, curve)),
                Err(_) => {
                    trace_decline("record fit failed", k, 0.0);
                    return None;
                }
            }
        }
    }
    if records.is_empty() {
        // Structure detected but everything clipped away: the level
        // set lies outside the queried cylinder window — an exact
        // empty answer.
        return Some(Vec::new());
    }
    // Honesty gate: the fitted records must lie on the cylinder
    // (analytic distance, serving-band scale) and their pcurves must
    // track the 3D curve within the same-parameter drift bound that
    // marched records live under (the 100×tol hull-gate scale — the
    // two fits are independent least squares, so their
    // parameterizations drift tangentially without either leaving
    // its surface). A multi-turn wall chart amplifies pcurve fit
    // error by its v-derivative (v ∈ [0,1] over the whole coil, so
    // |∂S/∂v| is the coil length): on a drift violation the pcurves
    // are refitted INTERPOLATING every exact sample, which removes
    // the least-squares truncation the amplification magnifies. A
    // record still off-band after that declines the whole arm.
    let surface_band = GRAZE_BAND_FACTOR * tolerance.linear();
    let drift_band = 100.0 * tolerance.linear();
    for (samples, record) in &mut records {
        if !record_on_surfaces(record, wall_first, wall, frame, surface_band, drift_band) {
            radial_level_curve_emit::interpolate_pcurves(samples, record);
            if !record_on_surfaces(record, wall_first, wall, frame, surface_band, drift_band) {
                trace_decline("honesty band", k, 0.0);
                return None;
            }
        }
    }
    if crate::trace::trace_enabled("radial-level") {
        let wraps = traces
            .iter()
            .map(|t| trace_wraps(t))
            .fold(0.0_f64, f64::max);
        eprintln!(
            "  [radial-level] served {} records ({k} crossings, {wraps:.2} wraps)",
            records.len()
        );
    }
    Some(
        records
            .into_iter()
            .map(|(_, record)| SurfaceSurfaceIntersection::Curve(record))
            .collect(),
    )
}

/// Gated decline diagnostics under the `radial-level` trace topic.
fn trace_decline(reason: &str, k: usize, wraps: f64) {
    if crate::trace::trace_enabled("radial-level") {
        eprintln!("  [radial-level] decline: {reason} (crossings {k}, wraps {wraps:.2})");
    }
}
