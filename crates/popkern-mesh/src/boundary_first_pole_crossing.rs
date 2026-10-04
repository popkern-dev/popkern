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
                // Segment i → pole, the pole row, pole → j; the pair's
                // span is shared at its midpoint parameter.
                let (edge, t_i, t_j) = ring.sources[i];
                let t_mid = 0.5 * (t_i + t_j);
                out.sources.push((edge, t_i, t_mid));
                let point: Point3 = surface.point_at(u, v_pole);
                let stations = (step.abs() / (TAU / 128.0)).ceil().max(1.0) as usize;
                for s in 0..=stations {
                    out.uv.push((u + step * s as f64 / stations as f64, v_pole));
                    out.pos.push(point);
                    degen.push(false);
                    out.sources.push(if s == stations {
                        (edge, t_mid, t_j)
                    } else {
                        (edge, t_mid, t_mid)
                    });
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
