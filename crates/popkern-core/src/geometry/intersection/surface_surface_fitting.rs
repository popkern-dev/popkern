//! Curve fitting for marched surface-surface intersection samples.

use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::surface::SurfaceGeometry;

#[path = "surface_surface_fit_samples.rs"]
mod surface_surface_fit_samples;
pub(crate) use surface_surface_fit_samples::{fit_march_samples, weld_closed_loop};

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
///
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
    // MULTI-WRAP records segment even on mixed pairs: the tuned
    // wrapped-record path handles at most ONE period jump (the
    // plane×cylinder ellipse class it was measured on); a trace that
    // wraps a periodic chart repeatedly — a thread-groove flank's
    // helical intersection with its shank crosses the cylinder seam
    // once per turn — fits a pcurve through every jump and the
    // split/merge strands the seam pieces between the crossings
    // (measured: 2-turn groove cut, 3 lethal strands incl. the seam
    // edge span between its two helix crossings, PK-460). Two or
    // more jumps on any periodic axis is the discriminator: the
    // single-jump tuned class is untouched.
    let multi_wrap = {
        let get = |s: &MarchSample, axis: usize| match axis {
            0 => s.u1,
            1 => s.v1,
            2 => s.u2,
            _ => s.v2,
        };
        // A trace is multi-wrap on an axis only when BOTH hold: it
        // jumps the period at least twice, AND its minimal-image
        // unwrapped excursion exceeds 1.5 periods. The jump count
        // alone misfires on GRAZING traces: a near-tangent contact's
        // chart registration is legitimately sloppy (±6e-2 measured on
        // the strip×cap hull-gate class) and jitters across the seam
        // repeatedly without ever traversing — segmenting it
        // shattered the contact curve and the operand-order pin lost
        // its analytic-first crossing. True multi-wrap traversal (a
        // thread helix crossing the shank seam once per turn)
        // accumulates whole periods of unwrapped excursion; jitter
        // accumulates none.
        let wraps = |axis: usize, period: Option<f64>| -> bool {
            let Some(p) = period.filter(|p| p.is_finite() && *p > 0.0) else {
                return false;
            };
            let jump_count = samples
                .windows(2)
                .filter(|w| (get(&w[1], axis) - get(&w[0], axis)).abs() > 0.5 * p)
                .count();
            if jump_count < 2 {
                return false;
            }
            let (mut acc, mut lo, mut hi) = (0.0_f64, 0.0_f64, 0.0_f64);
            for w in samples.windows(2) {
                let du = get(&w[1], axis) - get(&w[0], axis);
                acc += du - p * (du / p).round();
                lo = lo.min(acc);
                hi = hi.max(acc);
            }
            hi - lo > 1.5 * p
        };
        let (p1u, p1v) = surface1.periodicity();
        let (p2u, p2v) = surface2.periodicity();
        wraps(0, p1u) || wraps(1, p1v) || wraps(2, p2u) || wraps(3, p2v)
    };
    // OPEN seam-crossing traces on mixed pairs segment too: the march
    // normalizes each sample into the period, so an open trace that
    // walks across the seam (a cone-seam spiral on a drafted bore)
    // carries one genuine period jump the wrapped-record path cannot
    // absorb on a NURBS chart (PK-460: the wrapped cone-seam pcurve
    // left the chart and the record stranded).
    let open_seam_crossing = !is_closed && {
        // Same jitter guard as the multi-wrap gate: a jump only
        // counts when the trace's minimal-image unwrapped excursion
        // shows a REAL traversal (> half a period) — tangent-chart
        // registration jitters across the seam without traversing
        // (the grazing strip×cap operand-order pin regressed under a
        // bare jump-count gate, again).
        let crosses = |axis: usize, period: Option<f64>| -> bool {
            let Some(p) = period.filter(|p| p.is_finite() && *p > 0.0) else {
                return false;
            };
            let get = |s: &MarchSample| match axis {
                0 => s.u1,
                1 => s.v1,
                2 => s.u2,
                _ => s.v2,
            };
            let jumps = samples
                .windows(2)
                .filter(|w| (get(&w[1]) - get(&w[0])).abs() > 0.5 * p)
                .count();
            if jumps == 0 {
                return false;
            }
            let (mut acc, mut lo, mut hi) = (0.0_f64, 0.0_f64, 0.0_f64);
            for w in samples.windows(2) {
                let mut d = get(&w[1]) - get(&w[0]);
                d -= p * (d / p).round();
                acc += d;
                lo = lo.min(acc);
                hi = hi.max(acc);
            }
            hi - lo > 0.5 * p
        };
        let (p1u, p1v) = surface1.periodicity();
        let (p2u, p2v) = surface2.periodicity();
        crosses(0, p1u) || crosses(1, p1v) || crosses(2, p2u) || crosses(3, p2v)
    };
    let segments = if both_periodic || chart_closed_nurbs || multi_wrap || open_seam_crossing {
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
