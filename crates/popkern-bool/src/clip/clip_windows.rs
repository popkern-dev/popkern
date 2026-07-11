//! Window-SET clipping of intersection records to face boundaries.
//!
//! A record curve can enter and leave a face's trimmed domain several
//! times — an analytic sphere×plane circle crosses a cube face's
//! bounds four times, a line through a reflex outline twice per bite.
//! The per-face verdict is therefore a set of in-face parameter
//! windows ([`UvClipResult::Windows`]), and the joint verdict is the
//! pairwise intersection of the two faces' window sets (mod the
//! parameter period for closed records). Every surviving window
//! becomes its own clipped `CurveOnSurface` — collapsing them to one
//! span (the pre-PK-111 behavior) manufactured records that ran
//! THROUGH regions outside the face and over foreign junctions.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_uv::{
    UvClipResult, UvCrossing, clip_curve_to_crossings, clip_curve_wrapped, uv_clip_one_face,
};

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
        if resolved.iter().all(|&r| r) {
            return Some(vec![(t_min, t_max)]);
        }
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

/// Intersect two in-face window sets on the record's parameter line —
/// on its parameter CIRCLE for closed records, where either side's
/// windows may run through the seam (exit beyond the range end) and a
/// single window pair can overlap in two separate arcs.
fn intersect_window_sets(
    wa: &[(f64, f64)],
    wb: &[(f64, f64)],
    t_min: f64,
    t_max: f64,
    is_closed: bool,
) -> Vec<(f64, f64)> {
    let period = t_max - t_min;
    let full = |w: &(f64, f64)| w.1 - w.0 >= period - 1e-12;
    // A full-range window is no constraint at all: the other set
    // passes through unchanged (in particular a seam-running window
    // survives as ONE wrapped window, not two split arcs).
    if let [w] = wa
        && full(w)
    {
        return wb.to_vec();
    }
    if let [w] = wb
        && full(w)
    {
        return wa.to_vec();
    }
    let shifts: &[f64] = if is_closed { &[-1.0, 0.0, 1.0] } else { &[0.0] };
    let mut out: Vec<(f64, f64)> = Vec::new();
    for a in wa {
        for b in wb {
            for &s in shifts {
                let e = a.0.max(b.0 + s * period);
                let x = a.1.min(b.1 + s * period);
                if x - e <= 1e-15 {
                    continue;
                }
                let (mut e, mut x) = (e, x);
                if is_closed {
                    let k = ((e - t_min) / period).floor();
                    e -= k * period;
                    x -= k * period;
                }
                out.push((e, x));
            }
        }
    }
    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let eps = period.abs().max(1e-12) * 1e-9;
    out.dedup_by(|b, a| (b.0 - a.0).abs() < eps && (b.1 - a.1).abs() < eps);
    out
}

/// Clip an intersection record to both faces' trimmed UV domains.
/// Returns EVERY in-face span as its own `CurveOnSurface` (empty =
/// the record misses at least one face entirely). A record judged
/// in-face over its whole closed range survives as a single closed
/// clone.
pub(crate) fn clip_curve_uv_space(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
) -> Vec<CurveOnSurface> {
    let trace = std::env::var_os("PK_CLIP_TRACE").is_some();
    if trace {
        trace_record_geometry(brep_b, face_a, face_b, curve);
    }
    let clip_a = uv_clip_one_face(brep_a, face_a, &curve.curve_uv1, curve.is_closed);
    let clip_b = uv_clip_one_face(brep_b, face_b, &curve.curve_uv2, curve.is_closed);

    let name = |r: &UvClipResult| match r {
        UvClipResult::NoPcurves => "NoPcurves".to_string(),
        UvClipResult::Outside => "Outside".to_string(),
        UvClipResult::Windows(ws) => format!(
            "Windows[{}]",
            ws.iter()
                .map(|(e, x)| format!("({e:.3},{x:.3})"))
                .collect::<Vec<_>>()
                .join(",")
        ),
    };
    if matches!(clip_a, UvClipResult::Outside) || matches!(clip_b, UvClipResult::Outside) {
        if trace {
            eprintln!(
                "  [uvclip] dropped: A={} B={}",
                name(&clip_a),
                name(&clip_b)
            );
        }
        return Vec::new();
    }
    if matches!(clip_a, UvClipResult::NoPcurves) && matches!(clip_b, UvClipResult::NoPcurves) {
        return vec![curve.clone()];
    }
    if trace {
        eprintln!(
            "  [uvclip] {face_a:?}x{face_b:?}: A={} B={}",
            name(&clip_a),
            name(&clip_b)
        );
    }

    let (t_min, t_max) = curve.curve_3d.parameter_range();
    let windows_of = |r: UvClipResult| -> Vec<(f64, f64)> {
        match r {
            UvClipResult::Windows(ws) => ws,
            // NoPcurves on one side: that face imposes no constraint.
            _ => vec![(t_min, t_max)],
        }
    };
    let wa = windows_of(clip_a);
    let wb = windows_of(clip_b);
    let joint = intersect_window_sets(&wa, &wb, t_min, t_max, curve.is_closed);
    if trace && joint.len() > 1 {
        eprintln!(
            "  [uvclip] {} in-face windows survive the joint clip",
            joint.len()
        );
    }

    let mut out = Vec::new();
    for (t_enter, t_exit) in joint {
        if t_exit - t_enter < 1e-15 {
            continue;
        }
        if curve.is_closed && (t_enter - t_min).abs() < 1e-12 && (t_exit - t_max).abs() < 1e-12 {
            out.push(curve.clone());
            continue;
        }
        let wrapped = t_exit > t_max + 1e-12;
        let extract = |c: &NurbsCurve| -> Option<NurbsCurve> {
            if wrapped {
                clip_curve_wrapped(c, t_enter, t_exit)
            } else {
                clip_curve_to_crossings(c, t_enter, t_exit)
            }
        };
        let (Some(clipped_3d), Some(clipped_uv1), Some(clipped_uv2)) = (
            extract(&curve.curve_3d),
            extract(&curve.curve_uv1),
            extract(&curve.curve_uv2),
        ) else {
            continue;
        };
        out.push(CurveOnSurface {
            curve_3d: clipped_3d,
            curve_uv1: clipped_uv1,
            curve_uv2: clipped_uv2,
            is_closed: false,
        });
    }
    out
}

/// PK_CLIP_TRACE preamble: the record's 3D/UV extents and the B
/// face's boundary bbox.
fn trace_record_geometry(
    brep_b: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
) {
    let (a3, b3) = curve.curve_3d.parameter_range();
    let (a1, b1) = curve.curve_uv1.parameter_range();
    let p3a = curve.curve_3d.point_at(a3);
    let p3b = curve.curve_3d.point_at(b3);
    let p1a = curve.curve_uv1.point_at(a1);
    let p1b = curve.curve_uv1.point_at(b1);
    let (a2, b2) = curve.curve_uv2.parameter_range();
    let p2a = curve.curve_uv2.point_at(a2);
    let p2b = curve.curve_uv2.point_at(b2);
    eprintln!(
        "    [uvspace] {face_a:?}x{face_b:?} 3d t[{a3:.1},{b3:.1}] ({:.1},{:.1},{:.1})->({:.1},{:.1},{:.1}) | uv1 t[{a1:.1},{b1:.1}] ({:.1},{:.1})->({:.1},{:.1}) | uv2 ({:.3},{:.1})->({:.3},{:.1})",
        p3a.x,
        p3a.y,
        p3a.z,
        p3b.x,
        p3b.y,
        p3b.z,
        p1a.x,
        p1a.y,
        p1b.x,
        p1b.y,
        p2a.x,
        p2a.y,
        p2b.x,
        p2b.y
    );
    let ring = crate::clip::face_uv_boundary(brep_b, face_b);
    let mut lo = (f64::MAX, f64::MAX);
    let mut hi = (f64::MIN, f64::MIN);
    for (_, pc) in ring.iter().flatten() {
        let (a, b) = pc.parameter_range();
        for i in 0..4 {
            let q = pc.point_at(a + (b - a) * i as f64 / 3.0);
            lo = (lo.0.min(q.x), lo.1.min(q.y));
            hi = (hi.0.max(q.x), hi.1.max(q.y));
        }
    }
    eprintln!(
        "    [uvspace] B-poly bbox ({:.3},{:.1})..({:.3},{:.1})",
        lo.0, lo.1, hi.0, hi.1
    );
}
