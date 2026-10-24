//! Chart-cell containment for faces on periodic surfaces (PK-304/305).
//!
//! A face of a periodic analytic surface occupies one CELL of the
//! (u, v) chart lattice; whether a curve span lies on the face is a
//! statement *modulo the lattice* in both axes. This module computes
//! the cell (continuous boundary ring + bounding box + periods) and
//! answers folded containment queries — replacing the ±π seam dance
//! and u-mean branch shifts for cell-resolvable faces. Design:
//! `docs/architecture/periodic-chart-cell-clipping.md`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// One boundary pcurve's chained sample run.
pub(crate) struct RingRun {
    he: HalfEdgeHandle,
    t0: f64,
    t1: f64,
    samples: Vec<(f64, f64)>,
}

/// One face's chart cell: its continuous boundary ring in a single
/// branch, the ring's bounding box, and the surface periods.
pub(crate) struct ChartCell {
    /// Boundary sample runs, chained to continuity (one run per
    /// pcurve; consecutive runs meet end-to-start), each tagged with
    /// its half-edge and the pcurve parameter range.
    ring: Vec<RingRun>,
    /// Bounding box of the ring.
    pub u_range: (f64, f64),
    /// Bounding box of the ring (v).
    pub v_range: (f64, f64),
    u_period: Option<f64>,
    v_period: Option<f64>,
    /// True when the ring closes in the plane (no net period
    /// winding) — the precondition for even-odd containment.
    closes: bool,
    /// True when the ring runs counter-clockwise (positive area).
    counter_clockwise: bool,
    /// Coincidence band for boundary riding (the model's linear
    /// tolerance, read in uv): see `chart_cell_riding`.
    band: f64,
}

/// Compute the chart cell for `face`, or `None` when the face is not
/// cell-resolvable (aperiodic surface, no pcurves, or a ring that
/// winds the lattice — e.g. a full-period wall — which keeps the
/// existing tuned machinery).
pub(crate) fn face_chart_cell(
    brep: &BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
) -> Option<ChartCell> {
    let sof = brep.faces.get(face)?.surface.as_ref()?;
    let (u_period, v_period) = sof.geometry.periodicity();
    // DOUBLY-periodic surfaces only (torus): singly-periodic charts
    // (cylinder, sphere, cone) keep the long-tuned crossing machinery
    // — routing them through cells regressed the flush-gusset union's
    // watertight/volume oracle (measured 2026-08-21), while the
    // measured defect class the cell fixes (44% chart uncoverage,
    // double-claimed v-branch windows) is exclusive to two periodic
    // axes.
    if u_period.is_none() || v_period.is_none() {
        return None;
    }
    let boundary = crate::clip::face_uv_boundary(brep, face)?;
    if boundary.is_empty() {
        return None;
    }
    // Chain the ring to continuity with EXACT whole-period rebranching
    // only: a seam edge's second use sits an exact multiple of the
    // period from the walk (shared pcurve), while the arbitrary jumps
    // at a cone apex or sphere pole are NOT exact multiples — blind
    // rounding folds such contours and mis-classifies everything
    // (the measured cone/sphere lesson in `point_inside_face_uv`).
    // A ring with a non-exact jump is not cell-resolvable here.
    let exact_multiple = |delta: f64, period: Option<f64>| -> Option<f64> {
        match period {
            _ if delta.abs() < 1.0e-6 => Some(0.0),
            Some(p) if p.is_finite() && p > 0.0 => {
                let k = (delta / p).round();
                (k != 0.0 && (delta - k * p).abs() < 1.0e-6).then_some(k * p)
            }
            _ => None,
        }
    };
    let samples_per_edge = 32;
    let mut ring: Vec<RingRun> = Vec::new();
    let (mut u_lo, mut u_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut v_lo, mut v_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut prev_end: Option<(f64, f64)> = None;
    for (he, pcurve) in &boundary {
        let (t0, t1) = pcurve.parameter_range();
        if (t1 - t0).abs() < tolerance.linear() {
            continue;
        }
        let start = pcurve.point_at(t0);
        let (mut du, mut dv) = (0.0, 0.0);
        if let Some((pu, pv)) = prev_end {
            du = exact_multiple(pu - start.x, u_period)?;
            dv = exact_multiple(pv - start.y, v_period)?;
        }
        let mut run = Vec::with_capacity(samples_per_edge + 1);
        let mut prev_sample: Option<(f64, f64)> = None;
        for i in 0..=samples_per_edge {
            let t = t0 + (t1 - t0) * (i as f64) / (samples_per_edge as f64);
            let q = pcurve.point_at(t);
            // Intra-run unwrap: a stored pcurve may itself cross the
            // seam mid-run (the primitive torus's seam-straddling
            // quadrant faces do), wrapping its samples by a whole
            // period. Nearest-branch each sample to its predecessor —
            // measured: without this, the seam-band quadrants' rings
            // spanned 3/4 of the WRONG side of the chart and 44% of
            // the torus chart was covered by NO face cell, dropping
            // every seam-hugging intersection arc.
            let unwrap = |x: f64, prev: f64, period: Option<f64>| match period {
                Some(p) if p.is_finite() && p > 0.0 => x - ((x - prev) / p).round() * p,
                _ => x,
            };
            let (mut x, mut y) = (q.x + du, q.y + dv);
            if let Some((px, py)) = prev_sample {
                x = unwrap(x, px, u_period);
                y = unwrap(y, py, v_period);
            }
            prev_sample = Some((x, y));
            u_lo = u_lo.min(x);
            u_hi = u_hi.max(x);
            v_lo = v_lo.min(y);
            v_hi = v_hi.max(y);
            run.push((x, y));
        }
        prev_end = run.last().copied();
        ring.push(RingRun {
            he: *he,
            t0,
            t1,
            samples: run,
        });
    }
    if ring.is_empty() || !u_lo.is_finite() || !v_lo.is_finite() {
        return None;
    }
    // Closure: the chained ring's tail must meet its head without a
    // whole-period jump (net lattice winding zero). A winding ring
    // (full-period wall) is not a planar polygon — not cell-resolvable.
    let head = *ring.first().and_then(|r| r.samples.first())?;
    let tail = *ring.last().and_then(|r| r.samples.last())?;
    let gap_gate = |period: Option<f64>, span: f64| {
        period
            .map(|p| p * 0.25)
            .unwrap_or(span.abs().max(1.0) * 0.5)
    };
    let closes = (head.0 - tail.0).abs() < gap_gate(u_period, u_hi - u_lo)
        && (head.1 - tail.1).abs() < gap_gate(v_period, v_hi - v_lo);
    // A ring spanning (almost) the full period in an axis is treated
    // as winding for that axis — even-odd off the folded query is
    // then unreliable at the wrap.
    let winds = |period: Option<f64>, lo: f64, hi: f64| {
        period.map(|p| hi - lo >= p * 0.999).unwrap_or(false)
    };
    if !closes || winds(u_period, u_lo, u_hi) || winds(v_period, v_lo, v_hi) {
        return None;
    }
    // Enclosed-area guard: a full primitive wall's boundary is its
    // seam (twice) plus degenerate pole rims — the chained "ring"
    // collapses to a zero-area line (measured: the sphere wall's cell
    // came out u[6.283, 6.283] and every record dropped Outside). A
    // cell must enclose real chart area; degenerate rings keep the
    // existing machinery, which understands full-chart walls.
    let mut area2 = 0.0;
    let mut prev: Option<(f64, f64)> = None;
    let mut first: Option<(f64, f64)> = None;
    for run in &ring {
        for &p in &run.samples {
            if let Some(q) = prev {
                area2 += q.0 * p.1 - p.0 * q.1;
            }
            if first.is_none() {
                first = Some(p);
            }
            prev = Some(p);
        }
    }
    if let (Some(last), Some(head)) = (prev, first) {
        area2 += last.0 * head.1 - head.0 * last.1;
    }
    let min_area = 1.0e-6
        * u_period
            .or(v_period)
            .map(|p| p * p)
            .unwrap_or(1.0)
            .max(((u_hi - u_lo) * (v_hi - v_lo)).abs());
    if area2.abs() * 0.5 < min_area.max(1.0e-9) {
        return None;
    }
    Some(ChartCell {
        ring,
        u_range: (u_lo, u_hi),
        v_range: (v_lo, v_hi),
        u_period,
        v_period,
        closes: true,
        counter_clockwise: area2 > 0.0,
        band: tolerance.linear(),
    })
}

impl ChartCell {
    /// Fold a UV point into this cell's lattice branch (nearest
    /// whole-period image of each coordinate to the cell centre).
    pub(crate) fn fold(&self, u: f64, v: f64) -> (f64, f64) {
        let fold_axis = |x: f64, (lo, hi): (f64, f64), period: Option<f64>| match period {
            Some(p) if p.is_finite() && p > 0.0 => {
                let centre = 0.5 * (lo + hi);
                x - ((x - centre) / p).round() * p
            }
            _ => x,
        };
        (
            fold_axis(u, self.u_range, self.u_period),
            fold_axis(v, self.v_range, self.v_period),
        )
    }

    /// Folded containment: fold the point into the cell branch, then
    /// even-odd ray casting against the continuous ring.
    pub(crate) fn contains(&self, u: f64, v: f64) -> bool {
        debug_assert!(self.closes);
        let (fu, fv) = self.fold(u, v);
        // Cheap bbox reject with a tolerance-free strict box (the
        // even-odd pass decides marginal points).
        if fu < self.u_range.0 - 1e-9
            || fu > self.u_range.1 + 1e-9
            || fv < self.v_range.0 - 1e-9
            || fv > self.v_range.1 + 1e-9
        {
            return false;
        }
        // Even-odd over the ring segments (runs are continuous and
        // consecutive runs meet, so per-run windows cover the whole
        // ring; the head-tail closure segment is added explicitly).
        let mut crossings = 0usize;
        let mut test_segment = |a: (f64, f64), b: (f64, f64)| {
            if (a.1 > fv) != (b.1 > fv) {
                let f = (fv - a.1) / (b.1 - a.1);
                let x = a.0 + f * (b.0 - a.0);
                if x > fu {
                    crossings += 1;
                }
            }
        };
        let mut prev: Option<(f64, f64)> = None;
        let mut first: Option<(f64, f64)> = None;
        for run in &self.ring {
            for &p in &run.samples {
                if let Some(q) = prev {
                    test_segment(q, p);
                }
                if first.is_none() {
                    first = Some(p);
                }
                prev = Some(p);
            }
        }
        if let (Some(last), Some(head)) = (prev, first) {
            test_segment(last, head);
        }
        crossings % 2 == 1
    }
}

/// Windows of `curve_uv`'s parameter range that lie inside the cell,
/// by dense sampling + bisection refinement at the in/out transitions.
/// A stretch riding the cell boundary takes one half-open verdict
/// (`chart_cell_riding`). Empty when the curve never enters the cell.
pub(crate) fn cell_windows(
    cell: &ChartCell,
    curve_uv: &popkern_core::geometry::nurbs_curve::NurbsCurve,
) -> Vec<(f64, f64)> {
    use popkern_core::geometry::curve::CurveEval as _;
    let (t0, t1) = curve_uv.parameter_range();
    let n = (curve_uv.control_points.len() * 4).clamp(96, 1024);
    let ts: Vec<f64> = (0..=n)
        .map(|i| t0 + (t1 - t0) * (i as f64) / (n as f64))
        .collect();
    let nudges = chart_cell_riding::riding_nudges(cell, curve_uv, &ts);
    let at = |t: f64, nudge: Option<(f64, f64)>| {
        let q = curve_uv.point_at(t);
        let (du, dv) = nudge.unwrap_or((0.0, 0.0));
        cell.contains(q.x + du, q.y + dv)
    };
    let refine = |mut lo: f64, mut hi: f64, lo_in: bool, nudge: Option<(f64, f64)>| -> f64 {
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if at(mid, nudge) == lo_in {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    let mut windows: Vec<(f64, f64)> = Vec::new();
    let mut prev_in = at(ts[0], nudges[0]);
    let mut inside_since = prev_in.then_some(t0);
    for i in 1..=n {
        let now_in = at(ts[i], nudges[i]);
        if now_in != prev_in {
            // A transition at a riding stretch's end refines on the
            // stretch's own (nudged) verdict.
            let cross = refine(ts[i - 1], ts[i], prev_in, nudges[i - 1].or(nudges[i]));
            if prev_in {
                if let Some(start) = inside_since.take() {
                    windows.push((start, cross));
                }
            } else {
                inside_since = Some(cross);
            }
        }
        prev_in = now_in;
    }
    if let Some(start) = inside_since {
        windows.push((start, t1));
    }
    windows
}

/// The pair-clip entry: cell-fold windows for a cell-resolvable face,
/// `None` to keep the caller's tuned machinery. On a CLOSED record a
/// window through the record's parameter seam is ONE wrapped window
/// (exit beyond the range end, the joint clip's convention): the cell
/// of a face straddling the chart seam holds the arc of a u-circle
/// whose parameter starts at u = 0 as two windows, and the pieces met
/// at a junction no edge carries (measured, PK-942: a torus halved on
/// its equator under rotation_z(0.7) served 8 curveless edges there).
pub(crate) fn cell_clip(
    brep: &BRep,
    face: FaceHandle,
    uv: &popkern_core::geometry::nurbs_curve::NurbsCurve,
    is_closed: bool,
) -> Option<super::clip_uv::UvClipResult> {
    use popkern_core::geometry::curve::CurveEval as _;
    let cell = face_chart_cell(brep, face, &brep.tolerance)?;
    let mut windows = cell_windows(&cell, uv);
    let (t0, t1) = uv.parameter_range();
    if is_closed
        && windows.len() >= 2
        && let (Some(&(first_in, first_out)), Some(&(last_in, last_out))) =
            (windows.first(), windows.last())
        && first_in == t0
        && last_out == t1
    {
        windows.pop();
        windows[0] = (last_in, first_out + (t1 - t0));
    }
    if popkern_core::trace::trace_enabled("clip") {
        eprintln!(
            "  [cell-clip] {face:?} cell u[{:.3},{:.3}] v[{:.3},{:.3}] -> {} window(s)",
            cell.u_range.0,
            cell.u_range.1,
            cell.v_range.0,
            cell.v_range.1,
            windows.len()
        );
    }
    Some(if windows.is_empty() {
        super::clip_uv::UvClipResult::Outside
    } else {
        super::clip_uv::UvClipResult::Windows(windows)
    })
}

#[path = "chart_cell_crossings.rs"]
mod chart_cell_crossings;
pub(crate) use chart_cell_crossings::cell_crossings;

#[path = "chart_cell_riding.rs"]
mod chart_cell_riding;

#[path = "chart_cell_dissolve.rs"]
mod chart_cell_dissolve;
pub(crate) use chart_cell_dissolve::edge_dissolve_keeps_cell_resolvable;

#[cfg(test)]
mod tests {
    use super::*;

    /// Ground-truth dump + coverage check of the primitive torus's
    /// face cells: every chart point must be covered by some
    /// cell-resolvable face (ignored dump helper).
    #[test]
    fn primitive_torus_cells_cover_the_chart() {
        let tol = Tolerance::default();
        let torus = popkern_core::topology::BRep::torus(1.0, 0.35, tol).expect("torus");
        let mut cells = Vec::new();
        for (fh, _) in &torus.faces {
            match face_chart_cell(&torus, fh, &tol) {
                Some(c) => {
                    eprintln!(
                        "{fh:?}: u[{:.3},{:.3}] v[{:.3},{:.3}]",
                        c.u_range.0, c.u_range.1, c.v_range.0, c.v_range.1
                    );
                    cells.push((fh, c));
                }
                None => {
                    eprintln!("{fh:?}: NOT cell-resolvable");
                    if let Some(bnd) = crate::clip::face_uv_boundary(&torus, fh) {
                        for (he, pc) in &bnd {
                            let (t0, t1) = pc.parameter_range();
                            let p0 = pc.point_at(t0);
                            let p1 = pc.point_at(t1);
                            eprintln!(
                                "    he {he:?} ({:.4},{:.4})->({:.4},{:.4})",
                                p0.x, p0.y, p1.x, p1.y
                            );
                        }
                    }
                }
            }
        }
        // Chart coverage probe: every (u, v) grid point must be inside
        // at least one cell (fold-containment), else records there are
        // droppable by every face.
        let mut uncovered = 0;
        for i in 0..48 {
            for j in 0..48 {
                let u = (i as f64 + 0.5) / 48.0 * std::f64::consts::TAU;
                let v = (j as f64 + 0.5) / 48.0 * std::f64::consts::TAU;
                if !cells.iter().any(|(_, c)| c.contains(u, v)) {
                    uncovered += 1;
                    if uncovered <= 5 {
                        eprintln!("uncovered: u={u:.3} v={v:.3}");
                    }
                }
            }
        }
        assert_eq!(
            uncovered, 0,
            "{uncovered}/2304 chart points uncovered by face cells"
        );
    }
}
