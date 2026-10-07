//! Pole crossings: a sphere boundary curve that passes THROUGH a
//! chart pole (PK-646).
//!
//! A great-circle arc over the pole runs up one meridian and down the
//! opposite one; its samples near the pole sit half a period apart in
//! azimuth, and the chart segment between them is ambiguous — the
//! consecutive unwrap took ±π by rounding, and either way a straight
//! chart chord is not the curve (measured: jack_in_the_box body 29's
//! sphere corner, bounded by an arc of a tangent cylinder over the
//! pole, meshed 0.034 of OpenCASCADE's 0.261 and tore its three
//! neighbours open).
//!
//! In the chart the curve runs to the pole row, along it, and away:
//! every pole-row point is the pole, so the run costs nothing in 3-D.
//! The run's direction is the face's choice: going the short or the
//! long way round changes the ring's winding by one period, and the
//! two windings bound the two regions the curve splits the sphere
//! into. The kernel flag says the loop runs counter-clockwise in the
//! chart, so a winding-0 ring must run that way; otherwise the face
//! is the region a winding-±1 ring bounds, which the pole cap closes
//! through the other pole by the same flag.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::BoundaryRing;

/// Make every pole crossing of a sphere ring an explicit pole-row run
/// in the direction the face flag selects. The run's samples anchor
/// the later unwrap (`degenerate` is false for them), so it keeps the
/// choice. Returns whether the ring changed.
pub(in crate::boundary_first) fn resolve_pole_crossings(
    ring: &mut BoundaryRing,
    degenerate: &mut Vec<bool>,
    surface: &SurfaceGeometry,
    orientation: bool,
) -> bool {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};
    let SurfaceGeometry::Sphere { .. } = surface else {
        return false;
    };
    let n = ring.uv.len();
    if n < 3 || degenerate.len() != n {
        return false;
    }
    let regular: Vec<usize> = (0..n).filter(|&i| !degenerate[i]).collect();
    if regular.len() < 2 {
        return false;
    }
    let wrap = |d: f64| d - TAU * (d / TAU).round();
    // Consecutive regular pairs (cyclic), each with its short step and
    // the pole it crosses (None: not a crossing).
    let pairs: Vec<(usize, usize, f64, Option<f64>)> = (0..regular.len())
        .map(|k| {
            let (i, j) = (regular[k], regular[(k + 1) % regular.len()]);
            let d = wrap(ring.uv[j].0 - ring.uv[i].0);
            let (vi, vj) = (ring.uv[i].1, ring.uv[j].1);
            // A pole VERTEX between the pair (two rims meeting at the
            // pole: a corner patch's meridians) always runs the pole
            // row, whatever the azimuth gap: its chart corner IS that
            // row. Gating it on the over-the-pole test below decided
            // a quarter-sphere corner (gap exactly pi/2) by 1e-8 of
            // inversion noise, and the losing side charted the pole as
            // one point that cut the corner (a 25.8° in-face fold on
            // the boss plate; PK-626).
            if let Some(m) = degenerate_between(degenerate, i, j) {
                let pole = if ring.uv[m].1 < FRAC_PI_2 { 0.0 } else { PI };
                return (i, j, d, (d.abs() > POLE_VERTEX_GAP).then_some(pole));
            }
            let pole = if vi < FRAC_PI_4 && vj < FRAC_PI_4 {
                Some(0.0)
            } else if vi > PI - FRAC_PI_4 && vj > PI - FRAC_PI_4 {
                Some(PI)
            } else {
                None
            };
            (i, j, d, pole.filter(|_| d.abs() > FRAC_PI_2))
        })
        .collect();
    let crossings: Vec<usize> = (0..pairs.len()).filter(|&k| pairs[k].3.is_some()).collect();
    if crossings.is_empty() || crossings.len() > 3 {
        return false;
    }
    // Every choice of short/long way at each crossing; keep the one
    // the flag names.
    let mut chosen: Option<(BoundaryRing, Vec<bool>)> = None;
    let mut fallback: Option<(BoundaryRing, Vec<bool>)> = None;
    for mask in 0..(1u32 << crossings.len()) {
        let (candidate, degen, total) =
            assemble(ring, degenerate, surface, &pairs, &crossings, mask);
        let winding = (total / TAU).round() as i64;
        if winding == 0 {
            if (signed_area(&candidate.uv) > 0.0) == orientation && chosen.is_none() {
                chosen = Some((candidate, degen));
            }
        } else if winding.abs() == 1 && fallback.is_none() {
            fallback = Some((candidate, degen));
        }
    }
    match chosen.or(fallback) {
        Some((resolved, degen)) => {
            *ring = resolved;
            *degenerate = degen;
            true
        }
        None => false,
    }
}

/// The ring walked from its first regular sample with each pair's
/// short step, except at crossings, which run along the pole row the
/// short way (mask bit clear) or the long way (set). Returns the ring,
/// its degeneracy flags, and the total azimuth travelled.
fn assemble(
    ring: &BoundaryRing,
    degenerate: &[bool],
    surface: &SurfaceGeometry,
    pairs: &[(usize, usize, f64, Option<f64>)],
    crossings: &[usize],
    mask: u32,
) -> (BoundaryRing, Vec<bool>, f64) {
    use std::f64::consts::TAU;
    let n = ring.uv.len();
    let mut out = BoundaryRing {
        uv: Vec::new(),
        pos: Vec::new(),
        sources: Vec::new(),
    };
    let mut degen = Vec::new();
    let mut u = ring.uv[pairs[0].0].0;
    let mut total = 0.0;
    for (k, &(i, j, d, pole)) in pairs.iter().enumerate() {
        out.uv.push((u, ring.uv[i].1));
        out.pos.push(ring.pos[i]);
        degen.push(false);
        let step = match (pole, crossings.iter().position(|&c| c == k)) {
            (Some(v_pole), Some(bit)) => {
                let long = mask & (1 << bit) != 0;
                let step = if long { d - d.signum() * TAU } else { d };
                // Segment i → pole, the pole row, pole → j.
                let (into_pole, row, out_of_pole) = pole_run_sources(ring, degenerate, i, j);
                out.sources.push(into_pole);
                let point: Point3 = surface.point_at(u, v_pole);
                let stations = (step.abs() / (TAU / 128.0)).ceil().max(1.0) as usize;
                for s in 0..=stations {
                    out.uv.push((u + step * s as f64 / stations as f64, v_pole));
                    out.pos.push(point);
                    degen.push(false);
                    out.sources
                        .push(if s == stations { out_of_pole } else { row });
                }
                step
            }
            _ => {
                out.sources.push(ring.sources[i]);
                // Degenerate samples between i and j keep their place
                // (the unwrap rewrites their azimuth).
                let mut m = (i + 1) % n;
                while m != j {
                    out.uv.push(ring.uv[m]);
                    out.pos.push(ring.pos[m]);
                    out.sources.push(ring.sources[m]);
                    degen.push(degenerate[m]);
                    m = (m + 1) % n;
                }
                d
            }
        };
        u += step;
        total += step;
    }
    (out, degen, total)
}

/// Azimuth gap (radians) below which a pole vertex's two rims count
/// as one meridian: no pole row is needed between them.
const POLE_VERTEX_GAP: f64 = 1e-6;

/// The last chart-degenerate ring sample strictly between regular
/// samples `i` and `j` (cyclic), if any: a pole vertex.
fn degenerate_between(degenerate: &[bool], i: usize, j: usize) -> Option<usize> {
    let n = degenerate.len();
    let mut last = None;
    let mut m = (i + 1) % n;
    while m != j {
        if degenerate[m] {
            last = Some(m);
        }
        m = (m + 1) % n;
    }
    last
}

/// Canonical spans of a pole crossing's three chart pieces — regular
/// sample `i` to the pole, the pole row, the pole to regular sample
/// `j` — as `(into_pole, row, out_of_pole)`. The row's span is empty
/// (every row point is the pole), so it never requests a split.
///
/// - A curve passing OVER the pole (no ring sample at it) carries one
///   span `i → j`, shared at its midpoint parameter.
/// - A pole VERTEX (chart-degenerate samples between `i` and `j`,
///   where two boundary edges meet at the pole, e.g. the two meridian
///   rims of a sphere corner patch) has two spans, usually of two
///   different edges: `i`'s own span runs into the pole, and the last
///   degenerate sample's span runs out of it. Sharing `i`'s span on
///   both sides labelled the outgoing segment with the INCOMING edge,
///   so its long-segment splits landed on the incoming edge at the
///   pole and halved it every encroachment round, never refining the
///   outgoing one (PK-626: 4x-graded meridian samples to 3e-5 of the
///   edge, 1e-13-area slivers at the boss plate's corner patches).
fn pole_run_sources(
    ring: &BoundaryRing,
    degenerate: &[bool],
    i: usize,
    j: usize,
) -> (Span, Span, Span) {
    match degenerate_between(degenerate, i, j) {
        Some(p) => {
            let out_of_pole = ring.sources[p];
            let at_pole = (out_of_pole.0, out_of_pole.1, out_of_pole.1);
            (ring.sources[i], at_pole, out_of_pole)
        }
        None => {
            let (edge, t_i, t_j) = ring.sources[i];
            let t_mid = 0.5 * (t_i + t_j);
            ((edge, t_i, t_mid), (edge, t_mid, t_mid), (edge, t_mid, t_j))
        }
    }
}

/// One ring segment's canonical span (see `BoundaryRing::sources`).
type Span = (popkern_core::topology::EdgeHandle, f64, f64);

fn signed_area(uv: &[(f64, f64)]) -> f64 {
    let n = uv.len();
    (0..n)
        .map(|i| {
            let (a, b) = (uv[i], uv[(i + 1) % n]);
            a.0 * b.1 - b.0 * a.1
        })
        .sum::<f64>()
        * 0.5
}
