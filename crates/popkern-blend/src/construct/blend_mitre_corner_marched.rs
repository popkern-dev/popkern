//! Marched mitre geometry (PK-528): the mitre point and seams for strips
//! with no analytic common form. That covers mismatched-radius cylinders
//! (their axes miss one another) and curvature-continuous or
//! variable-radius NURBS strips.
//!
//! - **Mitre point.** Gauss–Newton on the strips' signed distances
//!   `(P − foot_i)·n_i`, each foot a hint-seeded projection, starting
//!   from the stale corner. It is accepted only when every strip lies
//!   within `pos_tol` of P, measured as the full point-to-foot
//!   distance so that a projection clamped at a patch border cannot
//!   pass. For N = 3 the solution is generically an isolated point.
//!   For N ≥ 4 it exists only at symmetric corners.
//! - **Seam.** Chord samples from the rail crossing to P, each
//!   refined onto both strips by alternating projection (the PK-160
//!   marcher). Every sample is gated on both strips, the 3-D curve is
//!   a degree-3 interpolant, and the fit is gated at every span
//!   midpoint. The sample count doubles (17 → 65) until the gate
//!   holds, or the seam declines. Pcurves project the same samples
//!   into each strip chart, seeded at the crossing's rail-neighbour
//!   chart point so that periodic charts stay on the corner's branch.

use popkern_core::geometry::fitting::interpolate_curve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, SurfaceEval, Vec3};

use crate::types::BlendSeamCornerDeclineReason as SeamReason;

use super::super::blend_seam_corner_march::authoring::{chart_pcurve, chord_params};
use super::super::blend_seam_corner_swept::{project_seeded, refine_onto_both, surface_distance};

/// Gauss–Newton cap for the mitre point. The solve converges
/// quadratically from the corner basin; the cap bounds a solve with no
/// common point.
const MAX_ITERATIONS: usize = 60;

/// Sample counts tried per seam, coarsest first.
const SEAM_SAMPLES: [usize; 3] = [17, 33, 65];

/// The point near `start` on every strip, or `None`.
pub(super) fn solve_general_mitre_point(
    strips: &[&SurfaceGeometry],
    start: Point3,
    reach: f64,
    pos_tol: f64,
) -> Option<Point3> {
    let mut hints: Vec<Option<(f64, f64)>> = vec![None; strips.len()];
    let mut p = start;
    for _ in 0..MAX_ITERATIONS {
        let mut jtj = [[0.0_f64; 3]; 3];
        let mut jtf = [0.0_f64; 3];
        let mut worst = 0.0_f64;
        for (geom, hint) in strips.iter().zip(hints.iter_mut()) {
            let (u, v) = project_seeded(geom, &p, *hint)?;
            *hint = Some((u, v));
            let foot = geom.point_at(u, v);
            let n = geom.normal_at(u, v).normalized()?;
            let f = (p - foot).dot(&n);
            worst = worst.max((p - foot).length());
            let row = [n.x, n.y, n.z];
            for i in 0..3 {
                jtf[i] += row[i] * f;
                for j in 0..3 {
                    jtj[i][j] += row[i] * row[j];
                }
            }
        }
        if worst <= pos_tol * 1e-3 {
            break;
        }
        let delta = super::point::solve3(jtj, [-jtf[0], -jtf[1], -jtf[2]])?;
        p = p + Vec3::new(delta[0], delta[1], delta[2]);
        if (p - start).length() > reach {
            return None;
        }
    }
    strips
        .iter()
        .zip(&hints)
        .all(|(geom, hint)| surface_distance(geom, &p, *hint) <= pos_tol)
        .then_some(p)
}

/// A marched seam: its samples from the rail crossing to the mitre
/// point, with their chord parameters.
pub(super) struct MarchedSeam {
    pub(super) samples: Vec<Point3>,
    pub(super) params: Vec<f64>,
    /// The seam's declared band (at least the crossing's recorded one).
    pub(super) band: f64,
}

/// Cap on a derived seam band: under the kernel's 1e-4 Precision
/// ceiling, as for every declared-accuracy record in the blend
/// pipeline (`final_declared_accuracy`, `edge_reconcile_fit`).
const DECLARED_ACCURACY_CAP: f64 = 9.9e-5;

/// March the seam between two strips from `crossing` to `mitre`.
///
/// `recorded` is the largest tolerance band recorded at the crossing
/// (at least `pos_tol`). Fitted strips (G2 or variable-radius NURBS)
/// meet their rails only within their fit accuracy: measured 6.9e-5
/// off a G2 strip at an oblique corner, whose rail recorded 6.8e-5.
/// A crossing on such geometry, whose recorded band exceeds `pos_tol`,
/// is accepted when its residual is sub-precision. The seam then
/// declares 1.5× that residual, capped under the Precision ceiling,
/// as its band. On exact geometry the crossing must meet `pos_tol`.
/// Interior samples converge onto the surfaces themselves and are
/// gated at `pos_tol`.
pub(super) fn march_mitre_seam(
    geom_a: &SurfaceGeometry,
    geom_b: &SurfaceGeometry,
    crossing: Point3,
    mitre: Point3,
    pos_tol: f64,
    recorded: f64,
) -> Result<MarchedSeam, SeamReason> {
    let residual =
        surface_distance(geom_a, &crossing, None).max(surface_distance(geom_b, &crossing, None));
    let band = if residual <= recorded {
        recorded
    } else if recorded > pos_tol && residual <= DECLARED_ACCURACY_CAP {
        (1.5 * residual).min(DECLARED_ACCURACY_CAP).max(recorded)
    } else {
        return Err(SeamReason::AnchorsOffSeam);
    };
    for count in SEAM_SAMPLES {
        let samples = seam_samples(geom_a, geom_b, crossing, mitre, count, pos_tol)?;
        let params = chord_params(&samples);
        let curve =
            interpolate_curve(&samples, &params, 3).map_err(|_| SeamReason::MarchUnresolved)?;
        let fits = params.windows(2).all(|w| {
            let q = curve.point_at(0.5 * (w[0] + w[1]));
            surface_distance(geom_a, &q, None) <= band && surface_distance(geom_b, &q, None) <= band
        });
        if fits {
            return Ok(MarchedSeam {
                samples,
                params,
                band,
            });
        }
    }
    Err(SeamReason::MarchUnresolved)
}

fn seam_samples(
    geom_a: &SurfaceGeometry,
    geom_b: &SurfaceGeometry,
    from: Point3,
    to: Point3,
    count: usize,
    pos_tol: f64,
) -> Result<Vec<Point3>, SeamReason> {
    let chord = (to - from)
        .normalized()
        .ok_or(SeamReason::MarchUnresolved)?;
    let mut hint_a = project_seeded(geom_a, &from, None);
    let mut hint_b = project_seeded(geom_b, &from, None);
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let f = i as f64 / (count - 1) as f64;
        let q = from + (to - from) * f;
        let p = if i > 0 && i < count - 1 {
            let start = out.last().copied().map_or(q, |prev: Point3| {
                // Start from the previous sample carried onto this
                // slice: the curve is continuous, so it lies closer
                // than the chord point where the seam bows.
                prev + chord * (q - prev).dot(&chord)
            });
            let p = slice_newton(
                geom_a,
                geom_b,
                start,
                q,
                chord,
                &mut hint_a,
                &mut hint_b,
                pos_tol,
            )
            .ok_or(SeamReason::MarchUnresolved)?;
            if surface_distance(geom_a, &p, hint_a) > pos_tol
                || surface_distance(geom_b, &p, hint_b) > pos_tol
            {
                return Err(SeamReason::MarchUnresolved);
            }
            p
        } else {
            q
        };
        out.push(p);
    }
    Ok(out)
}

/// Newton iteration cap for one seam sample.
const SLICE_ITERATIONS: usize = 40;

/// The point on both strips in the slice plane through `q` normal to
/// `chord`. Newton on `(sd_a, sd_b, (p − q)·chord)` converges
/// quadratically however small the angle between the strips. Near a
/// rail crossing both strips are tangent to the same support plane,
/// and alternating projection, whose rate is cos²θ, crawls there
/// (measured: stalled 1.2e-5..3.5e-5 off after 400 rounds). The
/// alternating marcher is the fallback when the Newton system is
/// singular.
#[allow(clippy::too_many_arguments)]
fn slice_newton(
    geom_a: &SurfaceGeometry,
    geom_b: &SurfaceGeometry,
    start: Point3,
    q: Point3,
    chord: Vec3,
    hint_a: &mut Option<(f64, f64)>,
    hint_b: &mut Option<(f64, f64)>,
    pos_tol: f64,
) -> Option<Point3> {
    let mut p = start;
    for _ in 0..SLICE_ITERATIONS {
        let (ua, va) = project_seeded(geom_a, &p, *hint_a)?;
        let (ub, vb) = project_seeded(geom_b, &p, *hint_b)?;
        *hint_a = Some((ua, va));
        *hint_b = Some((ub, vb));
        let n_a = geom_a.normal_at(ua, va).normalized()?;
        let n_b = geom_b.normal_at(ub, vb).normalized()?;
        let f = [
            (p - geom_a.point_at(ua, va)).dot(&n_a),
            (p - geom_b.point_at(ub, vb)).dot(&n_b),
            (p - q).dot(&chord),
        ];
        if f[0].abs().max(f[1].abs()) <= pos_tol * 1e-3 && f[2].abs() <= pos_tol {
            return Some(p);
        }
        let jac = [
            [n_a.x, n_a.y, n_a.z],
            [n_b.x, n_b.y, n_b.z],
            [chord.x, chord.y, chord.z],
        ];
        let Some(delta) = super::point::solve3(jac, [-f[0], -f[1], -f[2]]) else {
            return refine_onto_both(geom_a, geom_b, p, hint_a, hint_b, pos_tol);
        };
        p = p + Vec3::new(delta[0], delta[1], delta[2]);
    }
    refine_onto_both(geom_a, geom_b, p, hint_a, hint_b, pos_tol)
}

/// The seam's 3-D curve, running from `first` to the other end.
pub(super) fn seam_curve(seam: &MarchedSeam, reversed: bool) -> Result<CurveOnEdge, SeamReason> {
    let (pts, params) = oriented(seam, reversed);
    let curve = interpolate_curve(&pts, &params, 3).map_err(|_| SeamReason::MarchUnresolved)?;
    let (t0, t1) = curve.parameter_range();
    Ok(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve),
        t_start: t0,
        t_end: t1,
    })
}

/// The seam's pcurve in one strip chart. `crossing_uv` is the rail
/// neighbour's chart point at the crossing, which seeds the chart
/// branch.
pub(super) fn seam_pcurve(
    geom: &SurfaceGeometry,
    seam: &MarchedSeam,
    reversed: bool,
    crossing_uv: Option<(f64, f64)>,
) -> Result<Pcurve, SeamReason> {
    // Project from the crossing end so that the seed applies first.
    // Each projection is unwrapped onto its predecessor's period
    // branch (the seed's for the first). A projection may land on the
    // other side of a periodic chart's seam, as the crossing does when
    // it sits on the u = 0 / 2π rail: measured 3.7e-3 pcurve deviation
    // from a 2π jump inside one seam.
    let (period_u, period_v) = geom.parameter_periodicity();
    let unwrap = |value: f64, near: f64, period: Option<f64>| match period {
        Some(t) if t > 0.0 => value - ((value - near) / t).round() * t,
        _ => value,
    };
    let mut uv: Vec<(f64, f64)> = Vec::with_capacity(seam.samples.len());
    let mut hint = crossing_uv;
    for p in &seam.samples {
        let (mut u, mut v) =
            project_seeded(geom, p, hint).ok_or(SeamReason::PcurveClosureFailed)?;
        if let Some((hu, hv)) = hint {
            u = unwrap(u, hu, period_u);
            v = unwrap(v, hv, period_v);
        }
        hint = Some((u, v));
        uv.push((u, v));
    }
    let (uv, params) = if reversed {
        (
            uv.into_iter().rev().collect::<Vec<_>>(),
            seam.params
                .iter()
                .rev()
                .map(|t| 1.0 - t)
                .collect::<Vec<_>>(),
        )
    } else {
        (uv, seam.params.clone())
    };
    chart_pcurve(&uv, &params).map_err(|_| SeamReason::PcurveClosureFailed)
}

fn oriented(seam: &MarchedSeam, reversed: bool) -> (Vec<Point3>, Vec<f64>) {
    if reversed {
        (
            seam.samples.iter().rev().copied().collect(),
            seam.params.iter().rev().map(|t| 1.0 - t).collect(),
        )
    } else {
        (seam.samples.clone(), seam.params.clone())
    }
}
