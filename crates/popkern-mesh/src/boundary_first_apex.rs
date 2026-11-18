//! Apexes on the face chart: a cone's point, a geometric singularity
//! the atlas leaves to the face chart (chart-service.md §6a: the atlas
//! re-charts chart poles; an apex is a topological event, never
//! smoothed).
//!
//! A cone face bounded by one circle is its bounded side, the apex
//! (PK-636): an imported drill-point cone is one closed circle and
//! nothing else, its uv ring a full-period wrap no window charts. The
//! chart closes the ring through the apex: the unwrapped ring runs
//! `u0 → u0 ± p`, and an apex run returns along `v = 0` at the ring's
//! own u stations, every station the same 3-D point (the triangles on
//! that run are zero-area in 3-D and weld away). A ring that wraps a
//! full period through an apex vertex closes along the apex row
//! (PK-931).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::BoundaryRing;

/// A cone ring winding once round the cone, closed through its apex, or
/// `None` on another surface (or when the ring does not wrap exactly
/// once).
pub(super) fn close_round_apex(
    ring: &BoundaryRing,
    period: f64,
    surface: &SurfaceGeometry,
    chord_height: f64,
) -> Option<BoundaryRing> {
    let n = ring.uv.len();
    if n < 3 {
        return None;
    }
    let du = ring.uv[n - 1].0 - ring.uv[0].0;
    let winding = (du / period).round() as i64;
    if winding.abs() != 1 {
        return None;
    }
    // A cone face bounded by one circle is the BOUNDED side — the
    // apex; the other side runs to infinity. This holds whatever the
    // orientation flag says (measured: team.stp's drill-point cone
    // carries a flag the side rule reads backwards).
    let SurfaceGeometry::Cone { .. } = surface else {
        return None;
    };
    let v_apex = 0.0;
    let apex = surface.point_at(ring.uv[0].0, v_apex);
    // The ring does not repeat its first sample: close it one period
    // on (its last segment is the canonical span `sources[n − 1]`),
    // then return along the apex row over the ring's own u stations,
    // so the chart's apex row is as fine as its rim.
    let mut out = BoundaryRing {
        uv: ring.uv.clone(),
        pos: ring.pos.clone(),
        sources: ring.sources.clone(),
    };
    let closing = ring.sources[n - 1];
    let apex_span = (closing.0, closing.2, closing.2);
    let u_end = ring.uv[0].0 + winding as f64 * period;
    let v_rim = ring.uv[0].1;
    // The chart cut runs along one meridian, met from both sides: the
    // two walls carry the SAME v stations (a seam ladder, as in the
    // band slit) chosen so each wall chord stays within the chord
    // height of the surface — a single rim-to-apex chord flattened
    // the cap (measured: 8 % of a dome's volume).
    // A wall chord is also no longer than the ring's own longest
    // sample step: the wall's span is a zero-length apex run that the
    // mesher's segment refinement cannot split, and a straight cone
    // meridian has no sag to bisect, so one rim-to-apex segment's
    // protecting disk evacuated the cap's interior lattice and
    // Delaunay fanned chords through the solid (measured: jack body
    // 28's drill-point cone meshed 33.2 of 42.1, watertight).
    let longest_step = ring
        .pos
        .iter()
        .zip(ring.pos.iter().cycle().skip(1))
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max);
    let ladder = wall_stations(
        surface,
        ring.uv[0].0,
        (v_rim, v_apex),
        chord_height,
        longest_step,
    );
    out.uv.push((u_end, v_rim));
    out.pos.push(ring.pos[0]);
    out.sources.push(apex_span);
    for &v in &ladder {
        out.uv.push((u_end, v));
        out.pos.push(surface.point_at(u_end, v));
        out.sources.push(apex_span);
    }
    out.uv.push((u_end, v_apex));
    out.pos.push(apex);
    out.sources.push(apex_span);
    // The return run walks from u_end back to u0 strictly
    // monotonically: a ring holding runs of near-equal u that repeated
    // them left zero-length apex-row segments the coverage walk could
    // not cover (PK-646).
    let toward = (ring.uv[0].0 - u_end).signum();
    let step_floor = 1e-6 * period;
    let mut last_u = u_end;
    for k in (1..n).rev() {
        let u = ring.uv[k].0;
        if (u - last_u) * toward < step_floor || (ring.uv[0].0 - u) * toward < step_floor {
            continue;
        }
        out.uv.push((u, v_apex));
        out.pos.push(apex);
        out.sources.push(apex_span);
        last_u = u;
    }
    out.uv.push((ring.uv[0].0, v_apex));
    out.pos.push(apex);
    out.sources.push(apex_span);
    for &v in ladder.iter().rev() {
        out.uv.push((ring.uv[0].0, v));
        out.pos.push(surface.point_at(ring.uv[0].0, v));
        out.sources.push(apex_span);
    }
    Some(out)
}

/// Interior v stations of the meridian `u` from `v_rim` toward
/// `v_apex` (exclusive of both), bisected until every chord's
/// midpoint sits within `chord_height` of the surface and no chord is
/// longer than `longest_step`.
fn wall_stations(
    surface: &SurfaceGeometry,
    u: f64,
    (v_rim, v_apex): (f64, f64),
    chord_height: f64,
    longest_step: f64,
) -> Vec<f64> {
    let mut vs = vec![v_rim, v_apex];
    for _ in 0..12 {
        let mut next = vec![vs[0]];
        let mut split = false;
        for w in vs.windows(2) {
            let (a, b) = (surface.point_at(u, w[0]), surface.point_at(u, w[1]));
            let m = surface.point_at(u, 0.5 * (w[0] + w[1]));
            let mid = popkern_core::geometry::Point3::new(
                0.5 * (a.x + b.x),
                0.5 * (a.y + b.y),
                0.5 * (a.z + b.z),
            );
            if m.distance_to(&mid) > chord_height || a.distance_to(&b) > longest_step {
                next.push(0.5 * (w[0] + w[1]));
                split = true;
            }
            next.push(w[1]);
        }
        vs = next;
        if !split {
            break;
        }
    }
    vs[1..vs.len() - 1].to_vec()
}

/// A cone ring that wraps a full period through its APEX vertex and
/// that the structured revolution declined (PK-931): close it along
/// the apex row into a winding-0 ring the face flag orients, ready for
/// the general path. A cone minus a box whose corner edge runs down
/// the cone's axis leaves the cone wall rim -> seam up -> apex -> seam
/// down -> notch -> seam down: the notch makes the two seam sides
/// unequal, the revolution declines, and the band path (two wrapped
/// rings) has only one, so the wall fell to the uniform grid and meshed
/// 0.18 of its 2.05 volume share. `Err` hands the ring back unchanged.
pub(in crate::boundary_first) fn close_through_apex(
    mut ring: BoundaryRing,
    degenerate: &mut Vec<bool>,
    surface: &SurfaceGeometry,
    orientation: bool,
    period: f64,
) -> Result<Vec<BoundaryRing>, BoundaryRing> {
    use std::f64::consts::TAU;
    let SurfaceGeometry::Cone { .. } = surface else {
        return Err(ring);
    };
    let n = ring.uv.len();
    if n < 3 || degenerate.len() != n || (period - TAU).abs() > 1e-12 {
        return Err(ring);
    }
    let regular: Vec<usize> = (0..n).filter(|&i| !degenerate[i]).collect();
    if regular.len() < 2 {
        return Err(ring);
    }
    let wrap = |d: f64| d - TAU * (d / TAU).round();
    let pairs: Vec<(usize, usize, f64, Option<f64>)> = (0..regular.len())
        .map(|k| {
            let (i, j) = (regular[k], regular[(k + 1) % regular.len()]);
            let d = wrap(ring.uv[j].0 - ring.uv[i].0);
            let apex = degenerate_between(degenerate, i, j).map(|m| ring.uv[m].1);
            (i, j, d, apex)
        })
        .collect();
    let crossings: Vec<usize> = (0..pairs.len()).filter(|&k| pairs[k].3.is_some()).collect();
    if crossings.is_empty() || crossings.len() > 3 {
        return Err(ring);
    }
    // Each apex crossing runs the row by its short step or a full period
    // either way (a seam pair meets the apex with no azimuth gap).
    let mut chosen: Option<(BoundaryRing, Vec<bool>)> = None;
    for code in 0..3u32.pow(crossings.len() as u32) {
        let mut turns = vec![0i32; pairs.len()];
        let mut c = code;
        for &k in &crossings {
            turns[k] = (c % 3) as i32 - 1;
            c /= 3;
        }
        let (candidate, degen, total) = assemble_turns(&ring, degenerate, surface, &pairs, &turns);
        if (total / TAU).round() as i64 == 0 && (signed_area(&candidate.uv) > 0.0) == orientation {
            chosen = Some((candidate, degen));
            break;
        }
    }
    let Some((resolved, degen)) = chosen else {
        return Err(ring);
    };
    ring = resolved;
    *degenerate = degen;
    Ok(super::boundary_first_slit::dissolve_pinches(ring)
        .into_iter()
        .filter(|lobe| lobe.uv.len() >= 3)
        .collect())
}

/// [`assemble`] with each crossing's step its short step plus
/// `turns[k]` whole periods.
fn assemble_turns(
    ring: &BoundaryRing,
    degenerate: &[bool],
    surface: &SurfaceGeometry,
    pairs: &[(usize, usize, f64, Option<f64>)],
    turns: &[i32],
) -> (BoundaryRing, Vec<bool>, f64) {
    use std::f64::consts::TAU;
    let crossings: Vec<usize> = (0..pairs.len()).filter(|&k| pairs[k].3.is_some()).collect();
    let shifted: Vec<(usize, usize, f64, Option<f64>)> = pairs
        .iter()
        .zip(turns)
        .map(|(&(i, j, d, apex), &t)| (i, j, d + TAU * f64::from(t), apex))
        .collect();
    // `assemble` takes a crossing's step as given when its mask bit is
    // clear, so the shifted steps pass through unchanged.
    assemble(ring, degenerate, surface, &shifted, &crossings, 0)
}

/// The ring walked from its first regular sample with each pair's
/// short step, except at crossings, which run along the apex row the
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
    for (k, &(i, j, d, apex)) in pairs.iter().enumerate() {
        out.uv.push((u, ring.uv[i].1));
        out.pos.push(ring.pos[i]);
        degen.push(false);
        let step = match (apex, crossings.iter().position(|&c| c == k)) {
            (Some(v_apex), Some(bit)) => {
                let long = mask & (1 << bit) != 0;
                let step = if long { d - d.signum() * TAU } else { d };
                // Segment i → apex, the apex row, apex → j.
                let (into_apex, row, out_of_apex) = apex_run_sources(ring, degenerate, i, j);
                out.sources.push(into_apex);
                let point: Point3 = surface.point_at(u, v_apex);
                let stations = (step.abs() / (TAU / 128.0)).ceil().max(1.0) as usize;
                for s in 0..=stations {
                    out.uv.push((u + step * s as f64 / stations as f64, v_apex));
                    out.pos.push(point);
                    degen.push(false);
                    out.sources
                        .push(if s == stations { out_of_apex } else { row });
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

/// The last chart-degenerate ring sample strictly between regular
/// samples `i` and `j` (cyclic), if any: an apex vertex.
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

/// Canonical spans of an apex crossing's three chart pieces — regular
/// sample `i` to the apex, the apex row, the apex to regular sample
/// `j` — as `(into_apex, row, out_of_apex)`. The row's span is empty
/// (every row point is the apex), so it never requests a split.
///
/// - A curve passing OVER the apex (no ring sample at it) carries one
///   span `i → j`, shared at its midpoint parameter.
/// - A apex VERTEX (chart-degenerate samples between `i` and `j`,
///   where two boundary edges meet at the apex, e.g. the two meridian
///   rims of a sphere corner patch) has two spans, usually of two
///   different edges: `i`'s own span runs into the apex, and the last
///   degenerate sample's span runs out of it. Sharing `i`'s span on
///   both sides labelled the outgoing segment with the INCOMING edge,
///   so its long-segment splits landed on the incoming edge at the
///   apex and halved it every encroachment round, never refining the
///   outgoing one (PK-626: 4x-graded meridian samples to 3e-5 of the
///   edge, 1e-13-area slivers at the boss plate's corner patches).
fn apex_run_sources(
    ring: &BoundaryRing,
    degenerate: &[bool],
    i: usize,
    j: usize,
) -> (Span, Span, Span) {
    match degenerate_between(degenerate, i, j) {
        Some(p) => {
            let out_of_apex = ring.sources[p];
            let at_apex = (out_of_apex.0, out_of_apex.1, out_of_apex.1);
            (ring.sources[i], at_apex, out_of_apex)
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
