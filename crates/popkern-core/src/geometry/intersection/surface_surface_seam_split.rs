//! Seam-jump splitting of marched sample runs (extracted from
//! `surface_surface_fitting.rs` under the 500-line file budget).

use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::surface::SurfaceGeometry;

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
    // Periods include geometrically-closed NURBS charts (the march
    // wraps their seam; see `geometry::chart_periods`). Kernel NURBS
    // charts are clamped to [0, span], so seam positions at whole
    // multiples of the period hold for them exactly as for analytic
    // angle charts.
    let (p1u, p1v) = crate::geometry::chart_periods::geometric_chart_periods(surface1);
    let (p2u, p2v) = crate::geometry::chart_periods::geometric_chart_periods(surface2);
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
    let refine_crossing = |a: &MarchSample,
                           b: &MarchSample|
     -> Option<(MarchSample, MarchSample)> {
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
        let (t, axis) = best?;
        let lerp = |x: f64, y: f64| x + t * (y - x);
        let guess_point = crate::geometry::Point3::new(
            lerp(a.point.x, b.point.x),
            lerp(a.point.y, b.point.y),
            lerp(a.point.z, b.point.z),
        );
        let g: Vec<f64> = (0..4).map(|i| lerp(ca[i], cb_unwrapped[i])).collect();
        // Pin the jumping coordinate AT the seam: the free
        // projection slides along the intersection curve and left
        // the cut ~2e-4 off the seam, so neither piece's endpoint
        // matched the seam-edge hit node and micro pieces fit to
        // garbage (measured: the spring-grind seam sliver,
        // PK-443). Fall back to the free projection when the
        // pinned solve declines (curve tangent to the seam).
        let seam_value = {
            // Reconstruct the seam multiple the crossing search chose.
            let p = periods[axis]?;
            let lo = ca[axis].min(cb_unwrapped[axis]);
            (((lo / p).floor()) + 1.0) * p
        };
        let pin = crate::geometry::intersection::marching::BoundaryPin {
            on_surface_one: axis < 2,
            axis_u: axis % 2 == 0,
            value: seam_value,
        };
        let pinned =
            crate::geometry::intersection::marching::marching_newton::newton_boundary_crossing(
                surface1, g[0], g[1], surface2, g[2], g[3], pin, tolerance,
            );
        let trace = crate::trace::trace_enabled("seam-split");
        let (refined, via) = match pinned {
            Some(s)
                if s.point.distance_to(&guess_point)
                    <= a.point.distance_to(&b.point) * 4.0 + tolerance.linear() =>
            {
                (s, "pinned")
            }
            other => {
                if trace {
                    eprintln!(
                        "  [seam-split] pinned {} (axis {axis} seam {seam_value:.6})",
                        match &other {
                            Some(s) => format!(
                                "REJECTED gate: d={:.3e} vs gap={:.3e}",
                                s.point.distance_to(&guess_point),
                                a.point.distance_to(&b.point)
                            ),
                            None => "DECLINED".to_string(),
                        }
                    );
                }
                let dir = (b.point - a.point).normalized()?;
                let s = crate::geometry::intersection::marching::marching_newton::march_project(
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
                (s, "free")
            }
        };
        if trace {
            use crate::geometry::surface::SurfaceEval;
            let s1p = surface1.point_at(refined.u1, refined.v1);
            eprintln!(
                "  [seam-split] {via} crossing uv1 ({:.6},{:.6}) uv2 ({:.6},{:.6}) p ({:.6},{:.6},{:.6}) resid1 {:.2e}",
                refined.u1,
                refined.v1,
                refined.u2,
                refined.v2,
                refined.point.x,
                refined.point.y,
                refined.point.z,
                s1p.distance_to(&refined.point)
            );
        }
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
