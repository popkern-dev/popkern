//! Curve fitting for marched surface-surface intersection samples.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::fitting::fit_curve;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::surface::SurfaceGeometry;

/// Fit NURBS curves (3D + UV1 + UV2) to a sequence of march samples.
///
/// The fit is a pure function of the samples. The MARCHING caller
/// welds its closed loops with [`weld_closed_loop`] before fitting
/// (its tip gap is a marching artifact); analytic arms (spiric
/// sections, cross-drill branches) sample their loops exactly and
/// must NOT be welded — the perturbed tail parameterization trips
/// their sub-tolerance fitted-curve gates (measured: the spiric
/// dispatch pin declined to marching and emitted fragment chains).
pub(crate) fn fit_march_samples(
    samples: &[MarchSample],
    is_closed: bool,
) -> Result<CurveOnSurface, CoreError> {
    if samples.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "Too few samples to fit".into(),
        ));
    }

    // Chord-length parameterization
    let mut params = vec![0.0];
    let mut total_length = 0.0;
    for i in 1..samples.len() {
        total_length += samples[i].point.distance_to(&samples[i - 1].point);
        params.push(total_length);
    }
    if total_length < 1e-30 {
        return Err(CoreError::PreconditionFailed(
            "Degenerate curve (zero length)".into(),
        ));
    }
    // Normalize to [0, 1]
    for p in &mut params {
        *p /= total_length;
    }

    let degree = 3.min(samples.len() - 1);
    let m = samples.len();

    // 3D curve, residual-adaptive: a fixed 1-in-5 control ratio (and
    // a hard 200-point cap) under-resolved tight curvature — a blade
    // girdle rounding a 6e-3-radius edge fitted 3.7e-4 off the true
    // intersection and the imprinted edge failed the curve-surface
    // band (PK-440). Double the control count until the worst sample
    // residual is far inside any serving band (1e-7 of the curve's
    // own scale), up to full interpolation.
    let pts_3d: Vec<Point3> = samples.iter().map(|s| s.point).collect();
    let target = 1.0e-7 * total_length.max(1.0);
    let mut n_cp = (m / 5).clamp(degree + 1, m);
    let curve_3d = loop {
        use crate::geometry::curve::CurveEval as _;
        let c = fit_curve(&pts_3d, &params, degree, n_cp)?;
        let worst = pts_3d
            .iter()
            .zip(&params)
            .map(|(p, &t)| p.distance_to(&c.point_at(t)))
            .fold(0.0_f64, f64::max);
        if worst <= target || n_cp >= m {
            break c;
        }
        n_cp = (n_cp * 2).min(m);
    };

    // UV1 curve (using Point3(u1, v1, 0.0) convention)
    let pts_uv1: Vec<Point3> = samples
        .iter()
        .map(|s| Point3::new(s.u1, s.v1, 0.0))
        .collect();
    let curve_uv1 = fit_curve(&pts_uv1, &params, degree, n_cp)?;

    // UV2 curve
    let pts_uv2: Vec<Point3> = samples
        .iter()
        .map(|s| Point3::new(s.u2, s.v2, 0.0))
        .collect();
    let curve_uv2 = fit_curve(&pts_uv2, &params, degree, n_cp)?;

    Ok(CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed,
    })
}

/// Weld a marched closed loop: drop trailing samples within half an
/// average step of the start, then append the closing sample — the
/// exact start point, with the start's UVs when the loop also closes
/// in UV space (per surface, judged against the UV step scale) and
/// the tail's UVs otherwise (seam-winding loop: 3D closes, UV does
/// not).
pub(crate) fn weld_closed_loop(samples: &[MarchSample]) -> Vec<MarchSample> {
    let mut out = samples.to_vec();
    let first = out[0].clone();
    let avg_step = {
        let total: f64 = out
            .windows(2)
            .map(|w| w[0].point.distance_to(&w[1].point))
            .sum();
        (total / (out.len() - 1) as f64).max(1e-30)
    };
    while out.len() > 4 {
        let last = &out[out.len() - 1];
        if last.point.distance_to(&first.point) < avg_step * 0.5 {
            out.pop();
        } else {
            break;
        }
    }
    let tail = out[out.len() - 1].clone();
    let uv_gap = |a: (f64, f64), b: (f64, f64)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
    let uv_step = |pick: fn(&MarchSample) -> (f64, f64)| {
        let total: f64 = out
            .windows(2)
            .map(|w| uv_gap(pick(&w[0]), pick(&w[1])))
            .sum();
        (total / (out.len() - 1) as f64).max(1e-30)
    };
    let uv1_closes =
        uv_gap((tail.u1, tail.v1), (first.u1, first.v1)) < uv_step(|s| (s.u1, s.v1)) * 4.0;
    let uv2_closes =
        uv_gap((tail.u2, tail.v2), (first.u2, first.v2)) < uv_step(|s| (s.u2, s.v2)) * 4.0;
    out.push(MarchSample {
        point: first.point,
        u1: if uv1_closes { first.u1 } else { tail.u1 },
        v1: if uv1_closes { first.v1 } else { tail.v1 },
        u2: if uv2_closes { first.u2 } else { tail.u2 },
        v2: if uv2_closes { first.v2 } else { tail.v2 },
    });
    out
}

/// Split a marched sample run wherever a periodic UV coordinate jumps
/// by more than half its period between neighbours (PK-304): the
/// Newton steps normalize angles into `[0, period)`, so such a jump is
/// a seam crossing, and each emitted segment stays continuous on the
/// principal branch. Downstream fragment-junction welding (the split
/// stage's island/open-chain pre-passes) reassembles the pieces.
pub(crate) fn split_samples_at_period_jumps(
    samples: Vec<MarchSample>,
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    tolerance: &crate::tolerance::Tolerance,
) -> Vec<Vec<MarchSample>> {
    let (p1u, p1v) = surface1.periodicity();
    let (p2u, p2v) = surface2.periodicity();
    if p1u.is_none() && p1v.is_none() && p2u.is_none() && p2v.is_none() {
        return vec![samples];
    }
    let periods = [p1u, p1v, p2u, p2v];
    let coords = |s: &MarchSample| [s.u1, s.v1, s.u2, s.v2];
    let jumps = |a: &MarchSample, b: &MarchSample| -> bool {
        coords(a)
            .iter()
            .zip(coords(b))
            .zip(periods)
            .any(|((&x, y), period)| match period {
                Some(p) if p.is_finite() && p > 0.0 => (x - y).abs() > p * 0.5,
                _ => false,
            })
    };
    // Refine the exact seam crossing between the two straddling
    // samples and emit it as BOTH the tail's last sample and the
    // head's first sample (each on its own branch). Cutting at the
    // raw samples left every segment endpoint up to a marching step
    // (~1e-3) short of the seam — far beyond snap tolerance — so the
    // merge graph minted a dangling record node NEXT TO the exact
    // edge-pierce node, the fragment walk dropped the record as an
    // open chain, and the face served UNSPLIT with its disc fragment
    // missing (measured: torus × cylinder faces adjacent to the u=0
    // seam; 15 stranded half-edges in the assembled union).
    let refine_crossing =
        |a: &MarchSample, b: &MarchSample| -> Option<(MarchSample, MarchSample)> {
            let ca = coords(a);
            let cb = coords(b);
            // Unwrap b to a's branch and find the FIRST seam crossing.
            let mut best: Option<(f64, usize)> = None;
            let mut cb_unwrapped = cb;
            for i in 0..4 {
                let Some(p) = periods[i] else { continue };
                if !(p.is_finite() && p > 0.0) || (ca[i] - cb[i]).abs() <= p * 0.5 {
                    continue;
                }
                let shift = ((ca[i] - cb[i]) / p).round() * p;
                let bu = cb[i] + shift;
                cb_unwrapped[i] = bu;
                // Seam multiples of p strictly between ca[i] and bu.
                let lo = ca[i].min(bu);
                let hi = ca[i].max(bu);
                let k = (lo / p).floor() + 1.0;
                let seam = k * p;
                if seam <= lo || seam >= hi {
                    continue;
                }
                let t = (seam - ca[i]) / (bu - ca[i]);
                if best.is_none_or(|(bt, _)| t < bt) {
                    best = Some((t, i));
                }
            }
            let (t, _) = best?;
            let lerp = |x: f64, y: f64| x + t * (y - x);
            let guess_point = crate::geometry::Point3::new(
                lerp(a.point.x, b.point.x),
                lerp(a.point.y, b.point.y),
                lerp(a.point.z, b.point.z),
            );
            let g: Vec<f64> = (0..4).map(|i| lerp(ca[i], cb_unwrapped[i])).collect();
            let dir = (b.point - a.point).normalized()?;
            let refined = crate::geometry::intersection::marching::marching_newton::march_project(
                surface1,
                g[0],
                g[1],
                surface2,
                g[2],
                g[3],
                &guess_point,
                &dir,
                tolerance,
            )?;
            // Both copies carry the SAME 3D point; each side's UVs are
            // rebranched nearest its neighbour sample.
            let rebranch = |val: f64, near: f64, period: Option<f64>| match period {
                Some(p) if p.is_finite() && p > 0.0 => val - ((val - near) / p).round() * p,
                _ => val,
            };
            let for_side = |near: &MarchSample| MarchSample {
                u1: rebranch(refined.u1, near.u1, p1u),
                v1: rebranch(refined.v1, near.v1, p1v),
                u2: rebranch(refined.u2, near.u2, p2u),
                v2: rebranch(refined.v2, near.v2, p2v),
                point: refined.point,
            };
            Some((for_side(a), for_side(b)))
        };
    let mut segments: Vec<Vec<MarchSample>> = Vec::new();
    let mut current: Vec<MarchSample> = Vec::new();
    for sample in samples {
        if let Some(last) = current.last()
            && jumps(last, &sample)
        {
            let mut head_extra: Option<MarchSample> = None;
            if let Some((tail_x, head_x)) = refine_crossing(last, &sample) {
                current.push(tail_x);
                head_extra = Some(head_x);
            }
            segments.push(std::mem::take(&mut current));
            if let Some(hx) = head_extra {
                current.push(hx);
            }
        }
        current.push(sample);
    }
    if !current.is_empty() {
        segments.push(current);
    }
    segments
}

/// Segment a marched trace at UV period jumps, weld a closed loop's
/// marching tip gap, fit each piece, and push the results (PK-304).
///
/// Segmentation: the Newton steps normalize each sample's angles into
/// `[0, period)`, so a seam-crossing trace jumps by a full period
/// between neighbours; fitting ONE curve through the jump made the UV
/// track sweep the entire chart at the wrap parameter — downstream
/// face-window clipping saw phantom traversals through every patch
/// (measured: a torus × cylinder loop shed 0.002-wide sliver windows
/// into a patch on the OPPOSITE side of the torus, and the boolean
/// died in split). Principal-branch fragments are the same
/// chain-of-fragments shape the split stage's island machinery
/// already reassembles.
///
/// Weld: a closed loop's tip gap is marching slop (~1e-3), far above
/// the linear tolerance — without an exact weld, downstream
/// closed-curve detection treated the loop as open, param-sorted its
/// boundary crossings, and clipped a phantom exterior span. Applied
/// here, never inside [`fit_march_samples`]: analytic arms share the
/// fitter with exactly-closed loops that must stay untouched.
/// Rejoin the first and last segments of a seam-split CLOSED trace
/// across the seed/closure cut (see the call site in
/// `fit_marched_trace`). Merges only when the junction is
/// branch-continuous on every uv axis — a genuine seam jump there
/// (seed exactly on a seam) is left split.
fn rejoin_closed_trace_ends(mut segments: Vec<Vec<MarchSample>>) -> Vec<Vec<MarchSample>> {
    let Some(first) = segments.first().and_then(|s| s.first()) else {
        return segments;
    };
    let Some(last) = segments.last().and_then(|s| s.last()) else {
        return segments;
    };
    const BRANCH_SLACK: f64 = std::f64::consts::FRAC_PI_2;
    let continuous = (first.u1 - last.u1).abs() < BRANCH_SLACK
        && (first.v1 - last.v1).abs() < BRANCH_SLACK
        && (first.u2 - last.u2).abs() < BRANCH_SLACK
        && (first.v2 - last.v2).abs() < BRANCH_SLACK;
    if !continuous {
        return segments;
    }
    let head = segments.remove(0);
    match segments.last_mut() {
        Some(tail) => tail.extend(head),
        None => return vec![head],
    }
    segments.rotate_right(1);
    segments
}

pub(crate) fn fit_marched_trace(
    samples: &[MarchSample],
    is_closed: bool,
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    tolerance: &crate::tolerance::Tolerance,
    results: &mut Vec<crate::geometry::intersection::SurfaceSurfaceIntersection>,
) {
    // Segmentation is scoped to pairs where BOTH surfaces carry a
    // periodic axis (the measured defect class: torus × cylinder,
    // torus × torus). Mixed pairs (plane × torus, plane × cylinder)
    // keep their long-tuned wrapped-record path — the A/B pin
    // `filleted_operand_boolean_serves_identically_under_analytic_arms`
    // measures that segmenting them regresses serving (open loop +
    // radius-scale curve deviations). The chart-cell clip design
    // (docs/architecture/periodic-chart-cell-clipping.md) subsumes
    // this scoping.
    let both_periodic = {
        let (p1u, p1v) = surface1.periodicity();
        let (p2u, p2v) = surface2.periodicity();
        (p1u.is_some() || p1v.is_some()) && (p2u.is_some() || p2v.is_some())
    };
    let segments = if both_periodic {
        split_samples_at_period_jumps(samples.to_vec(), surface1, surface2, tolerance)
    } else {
        vec![samples.to_vec()]
    };
    // A CLOSED trace that seam-splits into multiple segments leaves
    // its first and last segments as two pieces of ONE arc, cut at
    // the arbitrary seed/closure point instead of a seam. The cut
    // point depends only on where seeding happened to land, and when
    // the seed falls near (but not on) a seam the tail piece is a
    // micro-segment that downstream record handling drops — leaving
    // a gap in the intersection loop and a half-claimed strand
    // decline (measured: torus × cylinder union under the exact
    // partial-window torus NURBS, whose seeds land at different
    // patch midpoints than the sampled fallback's; the loop lost a
    // 0.016-wide connector at the cylinder seam). Rejoin the ends
    // when their junction is branch-continuous; a seed exactly ON a
    // seam yields a near-period jump there and correctly skips.
    let segments = if is_closed && segments.len() > 1 {
        rejoin_closed_trace_ends(segments)
    } else {
        segments
    };
    // Drop degenerate blips: a jump adjacent to a sample already AT
    // the seam leaves a segment whose whole chord is one point (the
    // refined crossing plus its own near-duplicate). Fitted, it became
    // a zero-length shared record self-loop that the fragment walk
    // chewed into a mixed 10-he loop (measured: torus × cylinder face
    // adjacent to the u = 0 seam at the tangent junction).
    let min_chord = tolerance.linear();
    let segments: Vec<Vec<MarchSample>> = segments
        .into_iter()
        .filter(|seg| {
            let chord: f64 = seg
                .windows(2)
                .map(|w| w[0].point.distance_to(&w[1].point))
                .sum();
            chord > min_chord
        })
        .collect();
    let single_segment = segments.len() == 1;
    let segments: Vec<Vec<MarchSample>> = if is_closed && single_segment {
        segments
            .into_iter()
            .map(|seg| {
                if seg.len() > 4 {
                    weld_closed_loop(&seg)
                } else {
                    seg
                }
            })
            .collect()
    } else {
        segments
    };
    if crate::trace::trace_enabled("surface-intersection-segments") {
        for (i, seg) in segments.iter().enumerate() {
            if let (Some(f), Some(l)) = (seg.first(), seg.last()) {
                let m = &seg[seg.len() / 2];
                eprintln!(
                    "  [ssi-seg] {i}/{}: n={} uv1 ({:.3},{:.3})->({:.3},{:.3}) mid ({:.3},{:.3}) p ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                    segments.len(),
                    seg.len(),
                    f.u1,
                    f.v1,
                    l.u1,
                    l.v1,
                    m.u1,
                    m.v1,
                    f.point.x,
                    f.point.y,
                    f.point.z,
                    l.point.x,
                    l.point.y,
                    l.point.z
                );
            }
        }
    }
    for segment in segments {
        if segment.len() < 2 {
            continue;
        }
        if let Ok(curve_on_surface) = fit_march_samples(&segment, is_closed && single_segment) {
            results.push(
                crate::geometry::intersection::SurfaceSurfaceIntersection::Curve(curve_on_surface),
            );
        }
    }
}
