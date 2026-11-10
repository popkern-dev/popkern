//! A trim polygon whose loop runs along a chart pole row, its runs in
//! the sense the face flag fixes.
//!
//! At a pole (a sphere's, a cone's apex) every u is one point, so a
//! loop arriving on one branch and leaving on another runs along the
//! pole row for free, and how far is the loop's choice: the short way,
//! or a period more or less. The face flag decides it. The kernel's
//! loops keep the face on their left when the flag is set (counter-
//! clockwise in the chart), so along the top row (the face below it)
//! a run goes toward −u and along the bottom row toward +u, each within
//! one period; cleared, the senses swap. Folding each sample to its
//! predecessor's branch took the short way at every run instead: a
//! sphere face whose loop visits a pole three times (Viviani's cut with
//! the ball's pole at the double point) assembled a polygon winding its
//! lobes' tongues twice, the material-sense integral read the cut
//! −117.9 for about +2983 and its general glue declined (PK-1145). The
//! mesher charts the same runs the same way
//! (`boundary_first_pole_crossing::ruled_pole_runs`).

use std::f64::consts::TAU;

use super::self_intersection_uv_containment::{RawLoopSample, UvPoint};

/// Azimuth gap, as a fraction of the period, under which a run's two
/// ends sit on one meridian (a seam pair over the pole, or a cusp): a
/// slit or a whole period, decided by closure.
const ZERO_GAP: f64 = 1e-6;
/// Most zero-gap runs whose choices are enumerated.
const MAX_FREE_RUNS: u32 = 6;

/// The loop polygon of `raw` with every pole-row run (consecutive
/// samples on one degenerate row, `on_row`, including the closing pair)
/// taken in the sense the flag fixes, and every other step folded onto
/// its predecessor's branch. `None` when the loop runs along no pole
/// row, or no choice of the zero-gap runs closes the polygon with
/// signed area of the flag's sign (`ccw`).
pub(super) fn ruled_pole_polygon(
    raw: &[RawLoopSample],
    on_row: &dyn Fn(f64) -> bool,
    ccw: bool,
) -> Option<Vec<UvPoint>> {
    let n = raw.len();
    if n < 3 {
        return None;
    }
    let row_pair = |i: usize, j: usize| {
        let (a, b) = (raw[i].uv, raw[j].uv);
        on_row(a.v) && on_row(b.v) && (a.v - b.v).abs() <= 1e-6
    };
    // Each step i (from sample i − 1 to i, the closing step n from the
    // last sample to the first): a pole run's sense, or None.
    let sense = |i: usize| -> Option<f64> {
        let (from, to) = (i - 1, i % n);
        if !row_pair(from, to) {
            return None;
        }
        let row = raw[to].uv.v;
        // The face's side of the row: the first sample off it after the
        // run, else before it.
        let off = (1..n)
            .map(|k| raw[(to + k) % n].uv.v)
            .find(|&v| (v - row).abs() > 1e-6)?;
        let top = off < row;
        Some(if top == ccw { -1.0 } else { 1.0 })
    };
    let steps: Vec<Option<f64>> = (1..=n).map(sense).collect();
    if steps.iter().all(Option::is_none) {
        return None;
    }
    let wrap = |d: f64| d - TAU * (d / TAU).round();
    let free: Vec<usize> = (0..n)
        .filter(|&s| {
            steps[s].is_some() && wrap(raw[(s + 1) % n].uv.u - raw[s].uv.u).abs() <= ZERO_GAP * TAU
        })
        .collect();
    if free.len() as u32 > MAX_FREE_RUNS {
        return None;
    }
    (0..1u32 << free.len()).find_map(|code| {
        let mut polygon = vec![raw[0].uv];
        let mut u = raw[0].uv.u;
        for (s, step_sense) in steps.iter().enumerate() {
            let next = raw[(s + 1) % n].uv;
            let short = wrap(next.u - raw[s].uv.u);
            let du = match step_sense {
                None => short,
                Some(sign) => match free.iter().position(|&f| f == s) {
                    Some(bit) if code >> bit & 1 == 1 => sign * TAU,
                    Some(_) => 0.0,
                    None if short * sign > 0.0 => short,
                    None => short + sign * TAU,
                },
            };
            u += du;
            if s + 1 < n {
                polygon.push(UvPoint { u, v: next.v });
            }
        }
        // Closed: the walk returns to its start's u, round a region (a
        // slit at every zero-gap run has none).
        let closed = (u - raw[0].uv.u).abs() <= 1e-6 * TAU;
        let area = signed_area(&polygon);
        (closed && area.abs() > 1e-12 && (area > 0.0) == ccw).then_some(polygon)
    })
}

fn signed_area(polygon: &[UvPoint]) -> f64 {
    let n = polygon.len();
    0.5 * (0..n)
        .map(|i| {
            let (a, b) = (polygon[i], polygon[(i + 1) % n]);
            a.u * b.v - b.u * a.v
        })
        .sum::<f64>()
}

#[cfg(test)]
#[path = "self_intersection_uv_pole_runs_tests.rs"]
mod tests;
