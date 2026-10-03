//! Pole caps: a face whose single boundary ring winds once around a
//! chart POLE (a sphere's colatitude 0 or π, a cone's apex) with the
//! pole inside the face (PK-636).
//!
//! Native kernel caps carry the pole as topology (a seam edge, an apex
//! vertex) and reach the structured revolution path. An imported cap
//! is one closed circle and nothing else — a drill-point cone, a
//! domed knob — and its uv ring is a full-period wrap that no window
//! charts, so it declined to the fallback grid and meshed EMPTY
//! (measured: jack_in_the_box.stp, ten knob/ball bodies open by 32 or
//! 64 edges each).
//!
//! The chart closes the ring through the pole: the unwrapped ring runs
//! `u0 → u0 ± p`, and a pole run returns along `v = v_pole` at the
//! ring's own u stations, every station the same 3-D point (the
//! triangles on that run are zero-area in 3-D and weld away). Which
//! pole a sphere cap contains follows from the kernel's flag contract
//! (the loop runs counter-clockwise in the chart when the flag is
//! set, so the face lies left of the walk; a walk toward +u has +v on
//! its left). A cone cap is always the apex side (the only bounded
//! one).

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::BoundaryRing;

#[path = "boundary_first_pole_crossing.rs"]
mod pole_crossing;
#[path = "boundary_first_pole_runs.rs"]
mod pole_runs;
pub(super) use pole_crossing::resolve_pole_crossings;
pub(super) use pole_runs::{complement_through_poles, enclose_both_poles};

/// The ring closed through the pole the face contains, or `None` when
/// the surface has no pole on that side (or the ring does not wrap
/// exactly once).
pub(super) fn close_through_pole(
    ring: &BoundaryRing,
    period: f64,
    surface: &SurfaceGeometry,
    orientation: bool,
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
    // The kernel flag says the loop runs counter-clockwise in the
    // chart (the winding-honesty contract), so the face lies LEFT of
    // the walk when the flag is set; left of a walk toward +u is +v.
    let face_toward_larger_v = orientation == (winding > 0);
    let v_pole = match surface {
        SurfaceGeometry::Sphere { .. } => {
            if face_toward_larger_v {
                std::f64::consts::PI
            } else {
                0.0
            }
        }
        // A cone face bounded by one circle is the BOUNDED side — the
        // apex; the other side runs to infinity. This holds whatever
        // the orientation flag says (measured: team.stp's drill-point
        // cone carries a flag the side rule reads backwards).
        SurfaceGeometry::Cone { .. } => 0.0,
        _ => return None,
    };
    let pole = surface.point_at(ring.uv[0].0, v_pole);
    // The ring does not repeat its first sample: close it one period
    // on (its last segment is the canonical span `sources[n − 1]`),
    // then return along the pole row over the ring's own u stations,
    // so the chart's pole row is as fine as its rim.
    let mut out = BoundaryRing {
        uv: ring.uv.clone(),
        pos: ring.pos.clone(),
        sources: ring.sources.clone(),
    };
    let closing = ring.sources[n - 1];
    let pole_span = (closing.0, closing.2, closing.2);
    let u_end = ring.uv[0].0 + winding as f64 * period;
    let v_rim = ring.uv[0].1;
    // The chart cut runs along one meridian, met from both sides: the
    // two walls carry the SAME v stations (a seam ladder, as in the
    // band slit) chosen so each wall chord stays within the chord
    // height of the surface — a single rim-to-pole chord flattened
    // the cap (measured: 8 % of a dome's volume).
    // A wall chord is also no longer than the ring's own longest
    // sample step: the wall's span is a zero-length pole run that the
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
        (v_rim, v_pole),
        chord_height,
        longest_step,
    );
    out.uv.push((u_end, v_rim));
    out.pos.push(ring.pos[0]);
    out.sources.push(pole_span);
    for &v in &ladder {
        out.uv.push((u_end, v));
        out.pos.push(surface.point_at(u_end, v));
        out.sources.push(pole_span);
    }
    out.uv.push((u_end, v_pole));
    out.pos.push(pole);
    out.sources.push(pole_span);
    // The return run walks from u_end back to u0 strictly
    // monotonically: a ring that reaches a pole itself (a rim through
    // the other pole, PK-646) holds runs of near-equal u, and repeating
    // them left zero-length pole-row segments the coverage walk could
    // not cover.
    let toward = (ring.uv[0].0 - u_end).signum();
    let step_floor = 1e-6 * period;
    let mut last_u = u_end;
    for k in (1..n).rev() {
        let u = ring.uv[k].0;
        if (u - last_u) * toward < step_floor || (ring.uv[0].0 - u) * toward < step_floor {
            continue;
        }
        out.uv.push((u, v_pole));
        out.pos.push(pole);
        out.sources.push(pole_span);
        last_u = u;
    }
    out.uv.push((ring.uv[0].0, v_pole));
    out.pos.push(pole);
    out.sources.push(pole_span);
    for &v in ladder.iter().rev() {
        out.uv.push((ring.uv[0].0, v));
        out.pos.push(surface.point_at(ring.uv[0].0, v));
        out.sources.push(pole_span);
    }
    Some(out)
}

/// Interior v stations of the meridian `u` from `v_rim` toward
/// `v_pole` (exclusive of both), bisected until every chord's
/// midpoint sits within `chord_height` of the surface and no chord is
/// longer than `longest_step`.
fn wall_stations(
    surface: &SurfaceGeometry,
    u: f64,
    (v_rim, v_pole): (f64, f64),
    chord_height: f64,
    longest_step: f64,
) -> Vec<f64> {
    let mut vs = vec![v_rim, v_pole];
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
