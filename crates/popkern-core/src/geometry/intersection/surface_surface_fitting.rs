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
    let n_cp = (samples.len() / 5).clamp(degree + 1, 200);

    // 3D curve
    let pts_3d: Vec<Point3> = samples.iter().map(|s| s.point).collect();
    let curve_3d = fit_curve(&pts_3d, &params, degree, n_cp)?;

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
) -> Vec<Vec<MarchSample>> {
    let (p1u, p1v) = surface1.periodicity();
    let (p2u, p2v) = surface2.periodicity();
    if p1u.is_none() && p1v.is_none() && p2u.is_none() && p2v.is_none() {
        return vec![samples];
    }
    let jumps = |a: &MarchSample, b: &MarchSample| -> bool {
        let jump = |x: f64, y: f64, period: Option<f64>| match period {
            Some(p) if p.is_finite() && p > 0.0 => (x - y).abs() > p * 0.5,
            _ => false,
        };
        jump(a.u1, b.u1, p1u)
            || jump(a.v1, b.v1, p1v)
            || jump(a.u2, b.u2, p2u)
            || jump(a.v2, b.v2, p2v)
    };
    let mut segments: Vec<Vec<MarchSample>> = Vec::new();
    let mut current: Vec<MarchSample> = Vec::new();
    for sample in samples {
        if let Some(last) = current.last()
            && jumps(last, &sample)
        {
            segments.push(std::mem::take(&mut current));
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
pub(crate) fn fit_marched_trace(
    samples: &[MarchSample],
    is_closed: bool,
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
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
        split_samples_at_period_jumps(samples.to_vec(), surface1, surface2)
    } else {
        vec![samples.to_vec()]
    };
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
    if std::env::var_os("PK_SSI_SEG_TRACE").is_some() {
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
