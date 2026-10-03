//! Honest multi-window construction for the UV clip — split from
//! `clip_windows.rs` (500-line cap). The flush-rim and flush-span
//! arms (PK-393 and the PK-427 windowed sibling) live with it.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_uv::UvCrossing;

/// Honest multi-window construction from ≥3 boundary crossings: the
/// sorted crossings partition the curve's parameter line into
/// candidate intervals; each interval is classified by THREE interior
/// containment samples and kept only on a unanimous verdict.
/// Returns `None` when any interval classifies ambiguously
/// (boundary-hugging curves — containment verdicts are coin flips
/// there; the caller falls back to the tuned crossing-hull logic) so
/// this path never overrides a measured heuristic with noise.
///
/// Closed curves classify circularly (the wrap interval runs from the
/// last crossing to the first crossing plus one period); a window
/// through the seam keeps its exit beyond the range end, matching the
/// [`clip_curve_wrapped`] extraction contract. `Some(vec![])` means
/// every interval is decisively outside.
pub(super) fn decisive_windows(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    crossings: &[UvCrossing],
    is_closed: bool,
) -> Option<Vec<(f64, f64)>> {
    let (t_min, t_max) = curve_uv.parameter_range();
    let span = t_max - t_min;
    if !span.is_finite() || span <= 0.0 {
        return None;
    }
    // Tangential grazes arrive as near-equal crossing pairs; dedup so
    // they cannot fabricate zero-width intervals with flipped parity.
    let mut ts: Vec<f64> = crossings.iter().map(|c| c.t_curve).collect();
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    ts.dedup_by(|b, a| (*b - *a).abs() < span * 1e-7);
    if ts.len() < 2 {
        return None;
    }
    let mut intervals: Vec<(f64, f64)> = Vec::new();
    if !is_closed {
        intervals.push((t_min, ts[0]));
    }
    for w in ts.windows(2) {
        intervals.push((w[0], w[1]));
    }
    if is_closed {
        intervals.push((ts[ts.len() - 1], ts[0] + span));
    } else {
        intervals.push((ts[ts.len() - 1], t_max));
    }

    // Same containment convention as every tuned arm of
    // `uv_clip_one_face`: the record's own period branch, judged by
    // `point_inside_face_uv`.
    let inside_at = |t: f64| -> bool {
        let tt = if t > t_max { t - span } else { t };
        crate::split_inner_loop::point_inside_face_uv(
            brep,
            face,
            curve_uv.point_at(tt),
            &brep.tolerance,
        )
    };
    let min_len = span * 1e-9;
    // Boundary-hugging guard for OPEN records (PK-114): the flush
    // boolean's boundary-running traces classify some spans
    // decisively OUTSIDE while lying within a hair of the trim
    // boundary — dropping those spans starves the downstream
    // gap-substitute completion, which pairs full-span hull records
    // against boundary-registration gaps (measured: 46 open mesh
    // edges + 32% volume error on the flush gusset union). An
    // out-of-face interval only counts as decisive when its midpoint
    // is genuinely CLEAR of the face boundary (2% of the trim
    // scale); anything closer is ambiguous and falls back to the
    // tuned crossing hull. Closed records keep their PK-111 contract
    // untouched.
    let (period_u, period_v) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.periodicity())
        .unwrap_or((None, None));
    let boundary_clearance = super::clip_windows_clearance::boundary_clearance_samples(brep, face);
    let near_boundary_at = |t: f64| -> bool {
        let Some((boundary_samples, clearance)) = boundary_clearance.as_ref() else {
            return false;
        };
        let tt = if t > t_max { t - span } else { t };
        let p = curve_uv.point_at(tt);
        super::clip_windows_clearance::near_boundary(
            boundary_samples,
            *clearance,
            (p.x, p.y),
            (period_u, period_v),
        )
    };

    // TRANSVERSAL crossings (PK-632): the hugging guard below exists
    // for boundary-RUNNING traces — records that meet the trim at a
    // graze and then lie along it. An out-of-face interval bounded by
    // crossings where the record cuts the boundary at an angle is a
    // genuine gap between two in-face spans however close the rim is
    // (measured: a toothed prism's cap × the mating gear's flank line
    // crosses tooth flanks at ~40°; its tooth-gap interval sits within
    // the 2 % clearance of a flank, the guard fell back to the crossing
    // hull, the hull bridged two teeth through the gap, and the
    // intersection served open or declined). Transversality is the
    // sine of the uv angle between the record and the boundary pcurve
    // at the crossing; an unknown boundary (an inner-loop half-edge)
    // counts as grazing, keeping the guard.
    //
    // A gap is crossed from one boundary EDGE to ANOTHER (flank to
    // flank across a notch). An out-of-face interval that leaves and
    // rejoins the SAME edge is the record running alongside that edge
    // — the hugging class itself, however steep its crossings
    // (measured: the centrifugal impeller's vane-wall spirals leave
    // the disc's rim circle at sin 0.22–0.29 and run outside it within
    // 0.66 of the clearance band before re-crossing the same circle;
    // all 70 such intervals were same-edge, all 272 toothed-prism gaps
    // were edge-to-edge). Same-edge intervals keep the guard.
    let boundary = crate::clip::face_uv_boundary(brep, face);
    let transversal = |c: &UvCrossing| -> bool {
        use popkern_core::geometry::curve::CurveEval as _;
        let Some(pcurve) = boundary
            .as_ref()
            .and_then(|b| b.iter().find(|(h, _)| *h == c.boundary_he))
            .map(|(_, pc)| pc.to_nurbs())
        else {
            return false;
        };
        let tb = pcurve.tangent_at(c.t_boundary);
        let tr = curve_uv.tangent_at(c.t_curve);
        let (lb, lr) = (
            (tb.x * tb.x + tb.y * tb.y).sqrt(),
            (tr.x * tr.x + tr.y * tr.y).sqrt(),
        );
        if lb < 1e-12 || lr < 1e-12 {
            return false;
        }
        (tb.x * tr.y - tb.y * tr.x).abs() / (lb * lr) >= 0.17
    };
    let crossing_transversal_at = |t: f64| -> bool {
        let near: Vec<&UvCrossing> = crossings
            .iter()
            .filter(|c| (c.t_curve - t).abs() < span * 1e-7)
            .collect();
        !near.is_empty() && near.iter().all(|c| transversal(c))
    };

    let rejoins_same_edge = |a: f64, b: f64| -> bool {
        let at = |t: f64| {
            crossings
                .iter()
                .filter(move |c| (c.t_curve - t).abs() < span * 1e-7)
        };
        at(a).any(|ca| at(b).any(|cb| cb.boundary_he == ca.boundary_he))
    };

    let mut verdicts: Vec<Option<bool>> = Vec::with_capacity(intervals.len());
    for &(a, b) in &intervals {
        if b - a < min_len {
            // Endpoint-graze sliver: no room for an honest verdict;
            // it inherits its neighbor below.
            verdicts.push(None);
            continue;
        }
        let flags = [
            inside_at(a + (b - a) * 0.25),
            inside_at(a + (b - a) * 0.5),
            inside_at(a + (b - a) * 0.75),
        ];
        if flags[0] != flags[1] || flags[1] != flags[2] {
            return None;
        }
        // The interval's ends that are crossings (an open record's
        // first and last intervals end at the record's own end).
        let bounded_by_transversals = (a <= t_min || crossing_transversal_at(a))
            && (b >= t_max || crossing_transversal_at(b))
            && !rejoins_same_edge(a, b);
        if !is_closed
            && !flags[0]
            && near_boundary_at(a + (b - a) * 0.5)
            && !bounded_by_transversals
        {
            return None;
        }
        verdicts.push(Some(flags[0]));
    }
    // Slivers inherit the preceding decisive verdict (leading slivers
    // the following one) — at 1e-9 of the span the attribution is
    // immaterial, but runs must be contiguous.
    let first_decisive = verdicts.iter().flatten().next().copied()?;
    let mut resolved: Vec<bool> = Vec::with_capacity(verdicts.len());
    let mut prev = first_decisive;
    for v in &verdicts {
        let r = v.unwrap_or(prev);
        resolved.push(r);
        prev = r;
    }

    let n = intervals.len();
    let mut windows: Vec<(f64, f64)> = Vec::new();
    if is_closed {
        // FLUSH SPAN: an interval that classifies outside while its
        // midpoint hugs the trim boundary is a span running ALONG the
        // face's rim — the imprint curve of a tool face flush with
        // this one. PK-393 measured the whole-curve case (the z=0.4
        // circle of a flush-top-only cut classified 3 crossings →
        // all-outside on the tool lateral and was dropped; the host
        // top never got the hole and the wall self-intersected).
        // The WINDOWED case is a flush cap meeting this face over
        // only part of the record (measured: a nesting lip's corner
        // quarter-arc on a quarter-cylinder face — the closed circle
        // hugs the face's top edge for one quarter and is decisively
        // clear-outside elsewhere; dropping it left the containing
        // face's imprint ring an open chain and the union declined
        // with an unmatched trim vertex). On-boundary containment
        // verdicts are coin flips, so a hugging midpoint is treated
        // as on-face per interval; genuinely-outside spans have
        // clear midpoints and still drop.
        // The per-interval keep demands a TIGHT band — 100× linear
        // tolerance, not the 2 % heuristic clearance: a genuinely
        // flush span lies ON the boundary to roundoff, while tuned
        // holder/socket fixtures carry near-boundary spans well
        // inside the 2 % band that the downstream machinery expects
        // DROPPED (measured 2026-09-03: keeping them stranded two
        // fixtures in merge). The PK-393 whole-curve arm keeps its
        // original 2 % contract via `near_boundary_at` below.
        // And only when NO interval is decisively inside: a flush
        // record that also carries genuine in-face spans is the tuned
        // gap-substitute class (the fanned gusset's boundary-running
        // traces) — keeping its hugging spans as windows starves the
        // completion machinery exactly as PK-114 measured for open
        // records. The nesting-lip class has no inside spans at all.
        let any_inside = resolved.iter().any(|&r| r);
        let tight = brep.tolerance.linear() * 100.0;
        let tight_hug_at = |t: f64| -> bool {
            if any_inside {
                return false;
            }
            let Some((boundary_samples, _)) = boundary_clearance.as_ref() else {
                return false;
            };
            let tt = if t > t_max { t - span } else { t };
            let p = curve_uv.point_at(tt);
            super::clip_windows_clearance::near_boundary(
                boundary_samples,
                tight,
                (p.x, p.y),
                (period_u, period_v),
            )
        };
        let effective: Vec<bool> = intervals
            .iter()
            .zip(&resolved)
            .map(|(&(a, b), &r)| r || tight_hug_at(a + (b - a) * 0.5))
            .collect();
        // PK-393 whole-curve flush rim: every interval outside while
        // every midpoint hugs within the tuned 2 % clearance.
        if resolved.iter().all(|&r| !r)
            && intervals
                .iter()
                .all(|&(a, b)| near_boundary_at(a + (b - a) * 0.5))
        {
            return Some(vec![(t_min, t_max)]);
        }
        if effective.iter().all(|&r| r) {
            return Some(vec![(t_min, t_max)]);
        }
        let resolved = effective;
        let f = resolved.iter().position(|&r| !r)?;
        // Circular sweep from a decisively-outside interval; indices
        // past the array end continue one period up, where the wrap
        // interval's exit already points.
        let mut run: Option<(f64, f64)> = None;
        for i in f..f + n {
            let j = i % n;
            let lift = if i >= n { span } else { 0.0 };
            let (a, b) = (intervals[j].0 + lift, intervals[j].1 + lift);
            if resolved[j] {
                run = Some(run.map_or((a, b), |(ra, _)| (ra, b)));
            } else if let Some(w) = run.take() {
                windows.push(w);
            }
        }
        if let Some(w) = run {
            windows.push(w);
        }
        // Re-anchor each window's enter into the curve's own range;
        // a seam-running window keeps exit > t_max by design.
        for w in &mut windows {
            let k = ((w.0 - t_min) / span).floor();
            w.0 -= k * span;
            w.1 -= k * span;
        }
        windows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    } else {
        let mut run: Option<(f64, f64)> = None;
        for (i, &(a, b)) in intervals.iter().enumerate() {
            if resolved[i] {
                run = Some(run.map_or((a, b), |(ra, _)| (ra, b)));
            } else if let Some(w) = run.take() {
                windows.push(w);
            }
        }
        if let Some(w) = run {
            windows.push(w);
        }
    }
    windows.retain(|w| w.1 - w.0 > min_len);
    Some(windows)
}
