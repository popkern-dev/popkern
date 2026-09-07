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
fn rejoin_closed_trace_ends(
    mut segments: Vec<Vec<MarchSample>>,
    periods: [Option<f64>; 4],
) -> Vec<Vec<MarchSample>> {
    let Some(first) = segments.first().and_then(|s| s.first()) else {
        return segments;
    };
    let Some(last) = segments.last().and_then(|s| s.last()) else {
        return segments;
    };
    // Branch slack is PERIOD-AWARE: the fixed π/2 slack is meaningless
    // on a period-1 NURBS chart (a full-period junction jump of ~1.0
    // passed it, gluing a wrong-branch seed sample onto the tail —
    // PK-454). On a periodic axis anything approaching half a period
    // is a genuine seam jump, never junction slop.
    const BRANCH_SLACK: f64 = std::f64::consts::FRAC_PI_2;
    let slack = |p: Option<f64>| match p {
        Some(p) if p.is_finite() && p > 0.0 => BRANCH_SLACK.min(0.45 * p),
        _ => BRANCH_SLACK,
    };
    let continuous = (first.u1 - last.u1).abs() < slack(periods[0])
        && (first.v1 - last.v1).abs() < slack(periods[1])
        && (first.u2 - last.u2).abs() < slack(periods[2])
        && (first.v2 - last.v2).abs() < slack(periods[3]);
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
    //
    // A geometrically-closed NURBS chart is different from that
    // analytic scoping: its records CANNOT wrap (the pcurve would
    // leave the chart), and the march wraps its seam
    // (`geometry::chart_periods`), so any trace on such a surface must
    // split at the seam jumps — the splitter emits the refined seam
    // crossing on both branches, terminating each record exactly on
    // the seam (PK-443: a spring-grind section arc crossed the wall
    // seam mid-way; unsplit, the far side stayed on the wrong branch).
    let both_periodic = {
        let (p1u, p1v) = surface1.periodicity();
        let (p2u, p2v) = surface2.periodicity();
        (p1u.is_some() || p1v.is_some()) && (p2u.is_some() || p2v.is_some())
    };
    let chart_closed_nurbs = {
        let closed = |s: &SurfaceGeometry| {
            matches!(s, SurfaceGeometry::Nurbs(_)) && {
                let (u, v) = crate::geometry::chart_periods::geometric_chart_periods(s);
                u.is_some() || v.is_some()
            }
        };
        closed(surface1) || closed(surface2)
    };
    // MULTI-WRAP records segment even on mixed pairs: a trace that
    // WINDS a periodic chart repeatedly — a thread-groove flank's
    // helical intersection with its shank advances one full period
    // per turn — fits a pcurve through every seam jump and the
    // split/merge strands the seam pieces between the crossings
    // (measured: 2-turn groove cut, 3 lethal strands incl. the seam
    // edge span between its two helix crossings, PK-460). The
    // discriminator is NET WINDING (the unwrapped span exceeding
    // 1.5 periods), not jump count: the tuned wrapped-record classes
    // include CLOSED loops that cross the seam twice with zero net
    // winding (a plane×cylinder ellipse; the grazing strip×cap
    // operand-order pin regressed under a jump-count gate), while a
    // genuine multi-wrap trace accumulates a full period per turn.
    let multi_wrap = {
        let winds = |axis: usize, period: Option<f64>| -> bool {
            let Some(p) = period.filter(|p| p.is_finite() && *p > 0.0) else {
                return false;
            };
            let get = |s: &MarchSample| match axis {
                0 => s.u1,
                1 => s.v1,
                2 => s.u2,
                _ => s.v2,
            };
            let mut unwrapped = 0.0_f64;
            let mut lo = 0.0_f64;
            let mut hi = 0.0_f64;
            for w in samples.windows(2) {
                let mut d = get(&w[1]) - get(&w[0]);
                d -= p * (d / p).round();
                unwrapped += d;
                lo = lo.min(unwrapped);
                hi = hi.max(unwrapped);
            }
            hi - lo > 1.5 * p
        };
        let (p1u, p1v) = surface1.periodicity();
        let (p2u, p2v) = surface2.periodicity();
        winds(0, p1u) || winds(1, p1v) || winds(2, p2u) || winds(3, p2v)
    };
    let segments = if both_periodic || chart_closed_nurbs || multi_wrap {
        crate::geometry::intersection::surface_surface_seam_split::split_samples_at_period_jumps(
            samples.to_vec(),
            surface1,
            surface2,
            tolerance,
        )
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
    if crate::trace::trace_enabled("seam-split") {
        for (i, seg) in segments.iter().enumerate() {
            let (f, l) = (&seg[0], &seg[seg.len() - 1]);
            eprintln!(
                "  [fit-seg] pre-rejoin {i}: {} pts uv1 ({:.4},{:.4})->({:.4},{:.4}) 3d ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) is_closed={is_closed}",
                seg.len(),
                f.u1,
                f.v1,
                l.u1,
                l.v1,
                f.point.x,
                f.point.y,
                f.point.z,
                l.point.x,
                l.point.y,
                l.point.z
            );
        }
    }
    let segments = if is_closed && segments.len() > 1 {
        let combined = |s: &SurfaceGeometry| -> (Option<f64>, Option<f64>) {
            let (au, av) = s.periodicity();
            let (cu, cv) = crate::geometry::chart_periods::geometric_chart_periods(s);
            (au.or(cu), av.or(cv))
        };
        let (p1u, p1v) = combined(surface1);
        let (p2u, p2v) = combined(surface2);
        rejoin_closed_trace_ends(segments, [p1u, p1v, p2u, p2v])
    } else {
        segments
    };
    if crate::trace::trace_enabled("seam-split") {
        for (i, seg) in segments.iter().enumerate() {
            let (f, l) = (&seg[0], &seg[seg.len() - 1]);
            eprintln!(
                "  [fit-seg] post-rejoin {i}: {} pts uv1 ({:.4},{:.4})->({:.4},{:.4})",
                seg.len(),
                f.u1,
                f.v1,
                l.u1,
                l.v1
            );
        }
    }
    // C0 crease crossings cut the runs too: a smooth fit across a
    // faceted chart's corner rounds the true tangent break (PK-442 —
    // see `surface_surface_crease_split`).
    let segments =
        crate::geometry::intersection::surface_surface_crease_split::split_samples_at_c0_creases(
            segments, surface1, surface2, tolerance,
        );
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
    // Honesty gate: an emitted segment's points must lie ON the
    // surfaces. Seam-adjacent wrap jitter at a C0 seam CORNER (a
    // faceted prism wall's crease doubling as the chart seam —
    // PK-442) left micro segments evaluated on the EXTRAPOLATED facet
    // beyond the chart — provably off-surface, since a NURBS surface
    // lies inside its control-point hull. The hull test is the
    // discriminator: a tangential GRAZE record's uv registration is
    // legitimately sloppy (the strip × cap contact carries 6e-2 uv
    // slop on the tangent chart) but its points stay in-hull and the
    // record must survive. Drops are recorded under the
    // surface-intersection-segments trace topic.
    let hull_gate = tolerance.linear() * 100.0;
    let hull_excludes = |s: &SurfaceGeometry, p: &crate::geometry::Point3| -> bool {
        let SurfaceGeometry::Nurbs(ns) = s else {
            return false;
        };
        let mut lo = crate::geometry::Point3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut hi = crate::geometry::Point3::new(f64::MIN, f64::MIN, f64::MIN);
        for c in &ns.control_points {
            lo = crate::geometry::Point3::new(lo.x.min(c.x), lo.y.min(c.y), lo.z.min(c.z));
            hi = crate::geometry::Point3::new(hi.x.max(c.x), hi.y.max(c.y), hi.z.max(c.z));
        }
        p.x < lo.x - hull_gate
            || p.x > hi.x + hull_gate
            || p.y < lo.y - hull_gate
            || p.y > hi.y + hull_gate
            || p.z < lo.z - hull_gate
            || p.z > hi.z + hull_gate
    };
    let segments: Vec<Vec<MarchSample>> = segments
        .into_iter()
        .filter(|seg| {
            let off = [0, seg.len() / 2, seg.len() - 1].into_iter().any(|i| {
                hull_excludes(surface1, &seg[i].point) || hull_excludes(surface2, &seg[i].point)
            });
            if off && crate::trace::trace_enabled("surface-intersection-segments") {
                eprintln!(
                    "  [ssi-seg] DROP off-hull segment: n={} p ({:.4},{:.4},{:.4})",
                    seg.len(),
                    seg[seg.len() / 2].point.x,
                    seg[seg.len() / 2].point.y,
                    seg[seg.len() / 2].point.z
                );
            }
            !off
        })
        .collect();
    let single_segment = segments.len() == 1;
    let segments: Vec<Vec<MarchSample>> = if is_closed && single_segment {
        segments
            .into_iter()
            .map(|seg| {
                // The weld closes a marching TIP GAP (~one step). A
                // seam-REJOINED closed trace already terminates on the
                // pinned crossing at BOTH ends: its 3D ends coincide
                // exactly while its UVs deliberately differ by the
                // period branch. Welding it pops the exact crossing
                // sample ("trailing samples within half a step of the
                // start" — distance 0 qualifies) and re-appends the
                // closure with the PREVIOUS sample's uv, so the fitted
                // pcurve terminates short of the seam (measured,
                // PK-454: a 0.0027 chart gap at the z = 0.5 band
                // grind's seam junction tore the whole wall lattice —
                // 1023 open mesh edges on a validator-clean serve).
                let already_exact = seg
                    .first()
                    .zip(seg.last())
                    .is_some_and(|(f, l)| f.point.distance_to(&l.point) <= tolerance.linear());
                if seg.len() > 4 && !already_exact {
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
        eprintln!(
            "  [ssi-seg] is_closed={} both_periodic={} n_seg={}",
            is_closed,
            both_periodic,
            segments.len()
        );
        for (i, seg) in segments.iter().enumerate() {
            if let (Some(f), Some(l)) = (seg.first(), seg.last()) {
                let m = &seg[seg.len() / 2];
                eprintln!(
                    "  [ssi-seg] {i}/{}: n={} uv1 ({:.3},{:.3})->({:.3},{:.3}) mid ({:.3},{:.3}) p ({:.9},{:.9},{:.9})->({:.9},{:.9},{:.9})",
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
