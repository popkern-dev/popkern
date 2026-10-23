//! Exact containment for loops on a DOUBLY periodic chart (tori)
//! (PK-685).
//!
//! On a chart periodic in both directions no ray escapes the region:
//! a non-contractible loop (a band rim around the tube or the ring)
//! crosses any full-period ray an even number of times, so the
//! open-chart ray parity of `periodic_chart_parity` has no analogue.
//! Nor can the loops alone say which side is the region: a region and
//! its complement share their boundary. What decides it is the walk's
//! own handedness: the chart permutation walk emits every fragment
//! loop as a material-on-LEFT orbit of the σ-oriented chart
//! (`merge::walk::chart_walk`, σ self-calibrated from the face's stored
//! trim ring), so a point a short step to the σ-left of a boundary
//! segment, clear of every other segment, lies in the region by
//! construction.
//!
//! From that reference, membership of any point is exact. The region's
//! boundary (outer plus inner loops) is a mod-2 boundary, so a path's
//! crossing count with it has a parity that is independent of the
//! path: a closed path meets a bounding cycle an even number of times.
//! A point is in the region exactly when the straight lifted path from
//! the reference crosses the loops an even number of times. Loops whose
//! windings do not sum to zero mod 2 in both directions bound nothing
//! (an incomplete trim set), and the decision declines rather than
//! guess.
//!
//! A loop is a sequence of RUNS (one per half-edge). Inside a run,
//! segments are taken as the producer lifted them: a two-sample
//! endpoint run along a seam may legitimately span half a period or
//! more, and its short representative would reverse it (measured: a
//! half-torus band's seam pieces from v = π to 2π read as walking
//! −π, the windings came out odd and the region undecidable). Between
//! runs, where branch-glue jumps live, the junction segment is the
//! SHORT representative (`|Δu| ≤ period_u / 2`, `|Δv| ≤ period_v / 2`),
//! so a run authored on another period branch contributes its aligned
//! copy and a branch jump contributes nothing.

/// Chart periods `(period_u, period_v)` of a doubly periodic chart.
pub(crate) type Periods = (f64, f64);

/// Shift `d` by whole multiples of `period` into `[-period/2, period/2]`.
fn wrap(d: f64, period: f64) -> f64 {
    d - period * (d / period).round()
}

/// A closed chart loop: its runs' samples concatenated, with the index
/// of each run's first sample (see the module docs for how segments
/// inside and between runs are read).
#[derive(Debug, Clone, Default)]
pub(crate) struct ChartLoop {
    points: Vec<(f64, f64)>,
    run_start: Vec<bool>,
}

impl ChartLoop {
    /// Swap the loop's chart axes (a walk run on transposed axes,
    /// PK-922).
    pub(crate) fn transpose(&mut self) {
        for p in &mut self.points {
            std::mem::swap(&mut p.0, &mut p.1);
        }
    }

    /// The loop of `runs`, in walk order.
    pub(crate) fn from_runs<'a>(runs: impl IntoIterator<Item = &'a [(f64, f64)]>) -> Self {
        let mut points = Vec::new();
        let mut run_start = Vec::new();
        for run in runs {
            for (k, &p) in run.iter().enumerate() {
                points.push(p);
                run_start.push(k == 0);
            }
        }
        Self { points, run_start }
    }

    /// The loop's chart samples, in walk order.
    pub(crate) fn points(&self) -> &[(f64, f64)] {
        &self.points
    }
}

/// The loops' segments, each starting at its sample as stored: as-is
/// inside a run, the short representative at a run junction.
fn segments(loops: &[ChartLoop], (pu, pv): Periods) -> Vec<((f64, f64), (f64, f64))> {
    let mut out = Vec::new();
    for lp in loops {
        let n = lp.points.len();
        for i in 0..n {
            let j = (i + 1) % n;
            let (a, b) = (lp.points[i], lp.points[j]);
            let (du, dv) = if lp.run_start[j] {
                (wrap(b.0 - a.0, pu), wrap(b.1 - a.1, pv))
            } else {
                (b.0 - a.0, b.1 - a.1)
            };
            if du == 0.0 && dv == 0.0 {
                continue;
            }
            out.push((a, (a.0 + du, a.1 + dv)));
        }
    }
    out
}

/// Whether the loops bound a region on the torus: the net winding of
/// all loops together is even in both directions.
pub(crate) fn bounds_a_region(loops: &[ChartLoop], periods: Periods) -> bool {
    let (su, sv) = segments(loops, periods)
        .iter()
        .fold((0.0, 0.0), |(su, sv), (a, b)| {
            (su + b.0 - a.0, sv + b.1 - a.1)
        });
    let (wu, wv) = (
        (su / periods.0).round() as i64,
        (sv / periods.1).round() as i64,
    );
    wu % 2 == 0 && wv % 2 == 0
}

/// Squared distance from `p` to the segment `a`–`b`.
fn distance_squared_to_segment(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 {
        (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let (x, y) = (a.0 + t * dx - p.0, a.1 + t * dy - p.1);
    x * x + y * y
}

/// The lift of segment `s` whose midpoint lies within half a period of
/// `near` in both directions.
fn lift_near(
    s: ((f64, f64), (f64, f64)),
    near: (f64, f64),
    (pu, pv): Periods,
) -> ((f64, f64), (f64, f64)) {
    let mid = (0.5 * (s.0.0 + s.1.0), 0.5 * (s.0.1 + s.1.1));
    let (ku, kv) = (
        pu * ((near.0 - mid.0) / pu).round(),
        pv * ((near.1 - mid.1) / pv).round(),
    );
    ((s.0.0 + ku, s.0.1 + kv), (s.1.0 + ku, s.1.1 + kv))
}

/// A point inside the region: a step to the σ-LEFT of the longest
/// segment's midpoint (`sigma` the walk's chart handedness), half the
/// clearance to every other segment lift and at most a quarter of the
/// segment. Tries segments longest first; `None` when every candidate
/// is crowded below `1e-9` of its own length.
pub(crate) fn interior_reference(
    loops: &[ChartLoop],
    sigma: f64,
    periods: Periods,
) -> Option<(f64, f64)> {
    let segments = segments(loops, periods);
    let length =
        |s: &((f64, f64), (f64, f64))| ((s.1.0 - s.0.0).powi(2) + (s.1.1 - s.0.1).powi(2)).sqrt();
    let mut order: Vec<usize> = (0..segments.len()).collect();
    order.sort_by(|&i, &j| length(&segments[j]).total_cmp(&length(&segments[i])));
    for &i in &order {
        let (a, b) = segments[i];
        let len = length(&segments[i]);
        let mid = (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
        let clearance_sq = segments
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .flat_map(|(_, &s)| {
                let (p, q) = lift_near(s, mid, periods);
                // The nearest lift and its eight neighbours: a short
                // representative spans at most half a period per axis.
                (-1..=1).flat_map(move |du| {
                    (-1..=1).map(move |dv| {
                        let (ou, ov) = (du as f64 * periods.0, dv as f64 * periods.1);
                        distance_squared_to_segment(mid, (p.0 + ou, p.1 + ov), (q.0 + ou, q.1 + ov))
                    })
                })
            })
            .fold(f64::INFINITY, f64::min);
        let step = (0.25 * len).min(0.5 * clearance_sq.sqrt());
        if step <= 1e-9 * len {
            continue;
        }
        let left = (-(b.1 - a.1) / len, (b.0 - a.0) / len);
        return Some((mid.0 + sigma * step * left.0, mid.1 + sigma * step * left.1));
    }
    None
}

/// A point inside the region next to a SLIT: a boundary curve the
/// region's loops traverse in both directions (a seam edge, or a cut
/// walked in and back out), which therefore lies in the region's
/// interior with region on BOTH sides. `slit` is one traversal's chart
/// polyline; the point steps off its middle segment, to either side,
/// by half the clearance to every boundary segment other than the
/// slit's own two copies (at most a quarter of the segment). No walk
/// handedness is involved. `None` when the slit is degenerate or
/// crowded.
pub(crate) fn slit_reference(
    loops: &[ChartLoop],
    slit: &[(f64, f64)],
    periods: Periods,
) -> Option<(f64, f64)> {
    if slit.len() < 2 {
        return None;
    }
    let k = (slit.len() - 1) / 2;
    let (a, b) = (slit[k], slit[k + 1]);
    let len = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
    if len <= 0.0 {
        return None;
    }
    let mid = (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
    // The slit's own copies lie on the slit's chart line (any lift):
    // skip segments both of whose ends sit on it.
    let on_slit = |p: (f64, f64)| {
        slit.windows(2).any(|w| {
            (-1..=1).any(|du| {
                (-1..=1).any(|dv| {
                    let (ou, ov) = (du as f64 * periods.0, dv as f64 * periods.1);
                    distance_squared_to_segment(
                        p,
                        (w[0].0 + ou, w[0].1 + ov),
                        (w[1].0 + ou, w[1].1 + ov),
                    ) <= (1e-9 * len).powi(2)
                })
            })
        })
    };
    let clearance_sq = segments(loops, periods)
        .into_iter()
        .filter(|&(p, q)| !(on_slit(p) && on_slit(q)))
        .flat_map(|s| {
            let (p, q) = lift_near(s, mid, periods);
            (-1..=1).flat_map(move |du| {
                (-1..=1).map(move |dv| {
                    let (ou, ov) = (du as f64 * periods.0, dv as f64 * periods.1);
                    distance_squared_to_segment(mid, (p.0 + ou, p.1 + ov), (q.0 + ou, q.1 + ov))
                })
            })
        })
        .fold(f64::INFINITY, f64::min);
    let step = (0.25 * len).min(0.5 * clearance_sq.sqrt());
    if step <= 1e-9 * len {
        return None;
    }
    Some((
        mid.0 - step * (b.1 - a.1) / len,
        mid.1 + step * (b.0 - a.0) / len,
    ))
}

/// Signed area of the triangle `o`, `a`, `b` (twice).
fn orient(o: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
}

/// Whether `x` lies on the same side of the loops as `reference`: the
/// straight path from `reference` to the lift of `x` nearest it crosses
/// the loops' lifts an even number of times. Loop segments use the
/// half-open convention across the path's line, so a path through a
/// sample vertex counts once.
pub(crate) fn same_side(
    loops: &[ChartLoop],
    reference: (f64, f64),
    x: (f64, f64),
    periods: Periods,
) -> bool {
    let (pu, pv) = periods;
    let target = (
        reference.0 + wrap(x.0 - reference.0, pu),
        reference.1 + wrap(x.1 - reference.1, pv),
    );
    let centre = (
        0.5 * (reference.0 + target.0),
        0.5 * (reference.1 + target.1),
    );
    let mut crossings = 0usize;
    for s in segments(loops, periods) {
        let (p, q) = lift_near(s, centre, periods);
        for du in -1..=1 {
            for dv in -1..=1 {
                let (ou, ov) = (du as f64 * pu, dv as f64 * pv);
                let (a, b) = ((p.0 + ou, p.1 + ov), (q.0 + ou, q.1 + ov));
                let (sa, sb) = (orient(reference, target, a), orient(reference, target, b));
                if (sa > 0.0) == (sb > 0.0) {
                    continue;
                }
                let (sr, st) = (orient(a, b, reference), orient(a, b, target));
                if (sr > 0.0) != (st > 0.0) {
                    crossings += 1;
                }
            }
        }
    }
    crossings.is_multiple_of(2)
}

/// Where the loops cross the row `v`, folded into `[0, period_u)` and
/// sorted, coincident pairs cancelled (seam copies of one edge, or a
/// tangential touch, toggle twice). Half-open in `v`, so a row through
/// a sample vertex counts the crossing once.
pub(crate) fn row_crossings(loops: &[ChartLoop], v: f64, periods: Periods) -> Vec<f64> {
    let (pu, pv) = periods;
    let mut xs = Vec::new();
    for (a, b) in segments(loops, periods) {
        if a.1 == b.1 {
            continue;
        }
        let base = pv * ((v - 0.5 * (a.1 + b.1)) / pv).round();
        for k in -1..=1 {
            let shift = base + k as f64 * pv;
            let (va, vb) = (a.1 + shift, b.1 + shift);
            if (va > v) != (vb > v) {
                let u = a.0 + (v - va) / (vb - va) * (b.0 - a.0);
                xs.push(u - pu * (u / pu).floor());
            }
        }
    }
    crate::periodic_chart_parity::cancel_coincident_pairs(xs, pu * 1e-9)
}

/// The widest material arc of the row `v`: arcs between consecutive
/// row crossings (cyclic in `u`) are uniformly in or out, and each is
/// decided by [`same_side`] from `reference`. Returns the arc midpoint
/// `u`, or `None` when no arc of the row is material.
pub(crate) fn widest_material_arc(
    loops: &[ChartLoop],
    reference: (f64, f64),
    v: f64,
    periods: Periods,
) -> Option<f64> {
    let pu = periods.0;
    let xs = row_crossings(loops, v, periods);
    let arcs: Vec<(f64, f64)> = if xs.is_empty() {
        vec![(0.0, pu)]
    } else {
        (0..xs.len())
            .map(|i| (xs[i], xs.get(i + 1).copied().unwrap_or(xs[0] + pu)))
            .collect()
    };
    arcs.into_iter()
        .filter(|&(a, b)| same_side(loops, reference, (0.5 * (a + b), v), periods))
        .max_by(|x, y| (x.1 - x.0).total_cmp(&(y.1 - y.0)))
        .map(|(a, b)| 0.5 * (a + b))
}

#[cfg(test)]
#[path = "periodic_chart_parity_torus_tests.rs"]
mod tests;
