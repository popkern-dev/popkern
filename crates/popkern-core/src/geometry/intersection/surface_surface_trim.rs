//! Window trimming for marched surface-surface intersection traces.
//!
//! The march runs in windows padded by an entry band (see
//! `intersect_surfaces_with_handler`) so an out-of-window seed can walk
//! into the caller's windows — but the EMITTED curve must not overhang
//! them: callers key curve endpoints against their own trim windows.
//!
//! Dropping the outside samples alone leaves each trimmed end up to one
//! marching step SHORT of the window boundary. For analytic pairs the
//! windows carry 50% padding and the downstream UV clip computes exact
//! face crossings, so the slack is invisible. For a NURBS chart the
//! query window IS the chart, and a closed intersection loop cut at a
//! geometrically-closed chart's seam terminated ~5e-4 shy of the seam
//! on BOTH ends (measured: barrel wall × cube faces, PK-443) — the
//! merge graph minted two dangling endpoint nodes next to the exact
//! edge-pierce node and the interior-graze filter discarded the whole
//! loop. Refining each trimmed end onto the exact window boundary with
//! the pinned Newton solve makes the record terminate ON the seam
//! point, where endpoint node adoption finds the existing junction.

use crate::geometry::intersection::marching::marching_newton::newton_boundary_crossing;
use crate::geometry::intersection::marching::{BoundaryPin, MarchSample};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// Trim leading/trailing samples outside the caller's windows in a
/// non-periodic direction, refining each trimmed end onto the exact
/// window boundary when the pinned Newton solve converges. On solve
/// failure (curve tangent to the boundary, no convergence) the end
/// keeps the unrefined inside sample — the pre-existing behavior.
#[allow(clippy::too_many_arguments)]
pub(crate) fn trim_to_windows(
    samples: &[MarchSample],
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<MarchSample> {
    // Chart-aware periodicity: a geometrically-closed NURBS chart's
    // seam direction wraps in the march, so it is never trimmed here —
    // the period-jump splitter cuts those records at the seam with the
    // refined crossing instead.
    let (periodic_u1, periodic_v1) = crate::geometry::chart_periods::marching_periodicity(surface1);
    let (periodic_u2, periodic_v2) = crate::geometry::chart_periods::marching_periodicity(surface2);
    let outside = |s: &MarchSample| {
        let out =
            |x: f64, range: (f64, f64), periodic: bool| !periodic && (x < range.0 || x > range.1);
        out(s.u1, u1_range, periodic_u1)
            || out(s.v1, v1_range, periodic_v1)
            || out(s.u2, u2_range, periodic_u2)
            || out(s.v2, v2_range, periodic_v2)
    };
    let Some(start) = samples.iter().position(|s| !outside(s)) else {
        return Vec::new();
    };
    let end = samples.len() - samples.iter().rev().position(|s| !outside(s)).unwrap_or(0);

    let refine = |outside_sample: &MarchSample, inside_sample: &MarchSample| {
        refine_window_crossing(
            outside_sample,
            inside_sample,
            surface1,
            u1_range,
            v1_range,
            surface2,
            u2_range,
            v2_range,
            (periodic_u1, periodic_v1, periodic_u2, periodic_v2),
            tolerance,
        )
    };

    let trace = crate::trace::trace_enabled("trim-refine");
    let mut trimmed: Vec<MarchSample> = Vec::with_capacity(end - start + 2);
    if start > 0 {
        let refined = refine(&samples[start - 1], &samples[start]);
        if trace {
            eprintln!(
                "  [trim-refine] lead: cut {} sample(s), refine {}",
                start,
                if refined.is_some() { "ok" } else { "DECLINED" }
            );
        }
        if let Some(refined) = refined {
            trimmed.push(refined);
        }
    } else if trace {
        eprintln!(
            "  [trim-refine] lead: no outside samples (first uv2 ({:.6},{:.6}))",
            samples[0].u2, samples[0].v2
        );
    }
    trimmed.extend_from_slice(&samples[start..end]);
    if end < samples.len() {
        let refined = refine(&samples[end], &samples[end - 1]);
        if trace {
            eprintln!(
                "  [trim-refine] tail: cut {} sample(s), refine {}",
                samples.len() - end,
                if refined.is_some() { "ok" } else { "DECLINED" }
            );
        }
        if let Some(refined) = refined {
            trimmed.push(refined);
        }
    } else if trace {
        eprintln!(
            "  [trim-refine] tail: no outside samples (last uv2 ({:.6},{:.6}))",
            samples[samples.len() - 1].u2,
            samples[samples.len() - 1].v2
        );
    }
    trimmed
}

/// Refine the exact intersection∩window-boundary point between a
/// straddling (outside, inside) sample pair. Pins the first violated
/// non-periodic chart coordinate at its violated bound (a corner exit
/// violates two; refinement failure on the first falls back to the
/// unrefined trim). Returns `None` when the solve declines, when the
/// refined point wanders beyond the straddling gap (a different
/// intersection branch), or when it duplicates the inside sample.
#[allow(clippy::too_many_arguments)]
fn refine_window_crossing(
    outside_sample: &MarchSample,
    inside_sample: &MarchSample,
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    periodic: (bool, bool, bool, bool),
    tolerance: &Tolerance,
) -> Option<MarchSample> {
    let coords = [
        (outside_sample.u1, u1_range, periodic.0, true, true),
        (outside_sample.v1, v1_range, periodic.1, true, false),
        (outside_sample.u2, u2_range, periodic.2, false, true),
        (outside_sample.v2, v2_range, periodic.3, false, false),
    ];
    let pin = coords
        .iter()
        .find_map(|&(x, range, is_periodic, on_surface_one, axis_u)| {
            if is_periodic {
                return None;
            }
            let value = if x < range.0 {
                range.0
            } else if x > range.1 {
                range.1
            } else {
                return None;
            };
            Some(BoundaryPin {
                on_surface_one,
                axis_u,
                value,
            })
        })?;
    let refined = newton_boundary_crossing(
        surface1,
        inside_sample.u1,
        inside_sample.v1,
        surface2,
        inside_sample.u2,
        inside_sample.v2,
        pin,
        tolerance,
    )?;
    // The refined crossing must stay within the straddling gap's
    // neighbourhood — a solve that converged onto a DIFFERENT branch
    // of the intersection is worse than the unrefined trim.
    let gap = outside_sample.point.distance_to(&inside_sample.point);
    let dist = refined.point.distance_to(&inside_sample.point);
    if !dist.is_finite() || dist > gap * 4.0 + tolerance.linear() {
        return None;
    }
    // A duplicate of the inside sample would produce a zero chord in
    // the chord-length parameterization downstream.
    if dist < 1e-12 {
        return None;
    }
    // Wrap chart-periodic coordinates back into the window: the pin
    // constrains ONE coordinate, and the Newton solve may express a
    // free periodic coordinate on the continuation branch just past
    // the seam (measured: a corner crossing at a helical wall's rim ×
    // seam refined to u = −2e-4; without the wrap the sample run
    // carried no seam jump, the period splitter never cut the record
    // at the seam, and the overhanging sliver stranded the assembly —
    // PK-443). Wrapping restores the jump the splitter keys on.
    let mut refined = refined;
    let wrap = |x: &mut f64, range: (f64, f64), is_periodic: bool| {
        let span = range.1 - range.0;
        if is_periodic && span > 0.0 && (*x < range.0 || *x > range.1) {
            *x -= ((*x - range.0) / span).floor() * span;
        }
    };
    wrap(&mut refined.u1, u1_range, periodic.0);
    wrap(&mut refined.v1, v1_range, periodic.1);
    wrap(&mut refined.u2, u2_range, periodic.2);
    wrap(&mut refined.v2, v2_range, periodic.3);
    Some(refined)
}
