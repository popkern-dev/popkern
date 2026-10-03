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
//! pole the face contains follows from the loop convention: a loop
//! runs counter-clockwise about the face normal, so in the chart the
//! face lies to the LEFT of the walk when the face normal agrees with
//! `∂u × ∂v` and to the right otherwise; a walk toward +u has +v on its
//! left. A side without a pole there (the far side of a cone's apex)
//! keeps the decline.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::BoundaryRing;

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
    // Chart handedness against the face normal at a ring sample.
    let (u, v) = ring.uv[n / 2];
    let chart_normal = surface.partial_u(u, v).cross(&surface.partial_v(u, v));
    let face_normal = surface.normal_at(u, v) * if orientation { 1.0 } else { -1.0 };
    let agrees = chart_normal.dot(&face_normal);
    if agrees.abs() <= 1e-12 * chart_normal.length().max(1e-300) {
        return None;
    }
    let face_toward_larger_v = (agrees > 0.0) != (winding > 0);
    let v_pole = match surface {
        SurfaceGeometry::Sphere { .. } => {
            if face_toward_larger_v {
                std::f64::consts::PI
            } else {
                0.0
            }
        }
        SurfaceGeometry::Cone { .. } if !face_toward_larger_v => 0.0,
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
    let ladder = wall_stations(surface, ring.uv[0].0, v_rim, v_pole, chord_height);
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
    for k in (1..n).rev() {
        out.uv.push((ring.uv[k].0, v_pole));
        out.pos.push(pole);
        out.sources.push(pole_span);
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
/// midpoint sits within `chord_height` of the surface.
fn wall_stations(
    surface: &SurfaceGeometry,
    u: f64,
    v_rim: f64,
    v_pole: f64,
    chord_height: f64,
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
            if m.distance_to(&mid) > chord_height {
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
