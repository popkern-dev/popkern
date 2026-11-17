//! Periodic-branch unwrapping for boundary-first rings (split out
//! of `boundary_first.rs` for the 500-line ratchet).
//!
//! Pcurves on periodic charts can be authored in different period
//! branches; the unwrap walks the ring and shifts each sample to the
//! branch nearest its predecessor (positions are canonical 3-D points
//! and unaffected). Chart-degenerate samples (poles/apexes, arbitrary
//! azimuth) are skipped as anchors and rewritten from their flanking
//! regular samples. Returns `false` on a genuine net full wrap — the
//! caller routes those to the structured-revolution or uniform-grid paths.

use popkern_core::geometry::surface::{SurfaceEval as _, SurfaceGeometry};

use super::BoundaryRing;

/// Which ring samples are chart-degenerate in u (empty when u is not
/// periodic).
///
/// A sample AT a chart pole (sphere pole, cone apex) has
/// undefined azimuth: shifting u by half a period leaves
/// the 3D point fixed. Its pcurve u is arbitrary, and if
/// it participates in the consecutive unwrap it chains
/// the rest of the ring into a spurious ±period branch —
/// a phantom "full wrap" (flagship corner caps hold one
/// vertex exactly at the cap sphere's pole; the phantom
/// winding declined 2 of 4 caps). Degenerate samples are
/// skipped as unwrap anchors and their u is rewritten
/// from the flanking regular samples afterwards.
/// Gate at 1e-4 of the ring scale: a sample whose
/// half-period u-shift moves it by less than a
/// ten-thousandth of a boundary chord has no meaningful
/// azimuth for charting, and rewriting its u moves
/// geometry sub-chord. A pole VERTEX projected at
/// v ≈ 1e-6 shifts by ~2·r·v ≈ 1e-8 — above the old
/// 1e-7·scale gate — and its arbitrary atan2 u chained
/// PK-325's corner-cap rings into a phantom full wrap
/// (all four wavy-plate top caps declined to uniform-grid
/// sampling and their chains tore off the strips).
pub(super) fn chart_degenerate(
    ring: &BoundaryRing,
    surface: &SurfaceGeometry,
    period_u: Option<f64>,
    options: &crate::TessellationOptions,
) -> Vec<bool> {
    let ring_scale = ring
        .pos
        .iter()
        .zip(ring.pos.iter().cycle().skip(1))
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max)
        .max(1e-9);
    let Some(p) = period_u else {
        return Vec::new();
    };
    // Each sample's half-period shift: twice its distance from the pole.
    let shifts: Vec<f64> = ring
        .uv
        .iter()
        .zip(ring.pos.iter())
        .map(|(&(u, v), pos)| surface.point_at(u + 0.5 * p, v).distance_to(pos))
        .collect();
    let gate = 1.0e-4 * ring_scale;
    let gate = gate.max(pole_pass_gate(&shifts, gate, options.chord_height));
    shifts.iter().map(|&s| s < gate).collect()
}

/// A wider gate for a ring passing beside a pole within a twentieth of
/// the chord height: its samples within half the chord height of the pole
/// (a half-period shift under the chord height) are charted at the
/// pole, the chart's error there under the chord height, positions
/// kept. Such a pass turns its azimuth by up to half a period within
/// its distance `d` of the pole, at a chart height `d / R` against a
/// lattice cell; at distance `r` it has all but `asin(d / r)` of that
/// turn behind it, so past ten times `d` (half the chord height here)
/// the run along the pole row
/// meets the curve leaving it within a tenth of a radian. Gated at the
/// chord fraction alone, a ring passing 1e-5 beside a sphere's pole
/// with chords up to 1 charted only the samples within 5e-5 of it, and
/// the curve's turn just past them (0.33 rad at a chart height of
/// 1e-5) folded against the pole row's run: the face fell to the
/// fallback grid, covering part of itself twice (PK-1188). A pass is
/// a dip in the shift along the ring: a sample off the pole (a shift
/// past a millionth of the gate, rounding at the pole itself reading
/// 1e-17) no farther from it than either neighbour. An edge running
/// into the pole falls to its pole sample without a dip, its samples
/// halved toward the pole (a die's pip caps, their seams meeting the
/// pole at 3e-6 then 1e-18): gated at its nearest sample it charted
/// that edge's last stretch at the pole and left the cap open. A pass
/// farther out keeps the chord-fraction gate (its sweep the lattice
/// follows), as does a ring meeting the pole only at a vertex.
fn pole_pass_gate(shifts: &[f64], gate: f64, chord: f64) -> f64 {
    let n = shifts.len();
    let nearest = (0..n)
        .filter(|&i| {
            let s = shifts[i];
            s > 1.0e-6 * gate && s <= shifts[(i + n - 1) % n] && s <= shifts[(i + 1) % n]
        })
        .map(|i| shifts[i])
        .fold(f64::INFINITY, f64::min);
    if 10.0 * nearest <= chord { chord } else { 0.0 }
}

pub(super) fn unwrap_dir(
    ring: &mut BoundaryRing,
    period: f64,
    pick: fn(&mut (f64, f64)) -> &mut f64,
    degenerate: &[bool],
) -> bool {
    let is_degen = |i: usize| degenerate.get(i).copied().unwrap_or(false);
    let mut prev: Option<f64> = None;
    for i in 0..ring.uv.len() {
        if is_degen(i) {
            continue;
        }
        let cur = pick(&mut ring.uv[i]);
        if let Some(p) = prev {
            while *cur - p > 0.5 * period {
                *cur -= period;
            }
            while *cur - p < -0.5 * period {
                *cur += period;
            }
        }
        prev = Some(*cur);
    }
    // Rewrite each degenerate run's coordinate to the
    // mean of its flanking regular samples (circular):
    // any value evaluates to the same 3D point, and the
    // mean keeps both adjacent UV segments non-zero for
    // the Delaunay coverage walk.
    let n = ring.uv.len();
    for i in 0..n {
        if !is_degen(i) {
            continue;
        }
        let before = (1..n)
            .map(|k| (i + n - k) % n)
            .find(|&j| !is_degen(j))
            .map(|j| *pick(&mut ring.uv[j].clone()));
        let after = (1..n)
            .map(|k| (i + k) % n)
            .find(|&j| !is_degen(j))
            .map(|j| *pick(&mut ring.uv[j].clone()));
        if let (Some(b), Some(a)) = (before, after) {
            *pick(&mut ring.uv[i]) = 0.5 * (b + a);
        }
    }
    // Net winding across the regular samples = a genuine
    // full wrap; boundary-first cannot chart it.
    let regular: Vec<usize> = (0..n).filter(|&i| !is_degen(i)).collect();
    match (regular.first(), regular.last()) {
        (Some(&f), Some(&l)) => {
            let first = *pick(&mut ring.uv[f].clone());
            let last = *pick(&mut ring.uv[l].clone());
            (last - first).abs() <= 0.5 * period
        }
        _ => true,
    }
}
/// Put every hole ring on the outer ring's period sheet.
///
/// The unwrap makes each ring continuous on its own, but the sheet a
/// ring lands on is an accident of its first sample: a seam sample
/// whose projection reads u = 2π − ε (its ordinate a rounding below
/// the axis) carries the whole outer ring one period over, while the
/// holes, whose pcurves never touch the seam, stay on the base sheet.
/// The holes then fall OUTSIDE the outer polygon, the even-odd region
/// test admits the grid under them, and the face over-fills (measured
/// on a sphere-cube subtraction at the 1e-3 tolerance once its rims
/// were exact arcs: 546 triangles instead of 384, 29.1 units of UV
/// area inside a 16.1-unit ring, 70% of the body's volume lost).
/// Delaunay in the chart needs all rings on ONE sheet, so each hole is
/// shifted by the whole-period offset (u and v, ±2 periods) under
/// which a majority of its samples fall inside the outer ring; zero
/// shift wins ties. A ring no shift places inside keeps its sheet —
/// it may be a second boundary component rather than a hole.
/// Place a face's rings on one chart sheet: closed-NURBS rings move
/// into the knot domain, a sphere ring through both poles takes the
/// side its flag names, and holes join the outer ring's sheet.
pub(super) fn place_rings_on_chart(
    rings: &mut Vec<BoundaryRing>,
    surface: &SurfaceGeometry,
    orientation: bool,
    chord_height: f64,
) {
    use super::boundary_first_pole_cap::{complement_through_poles, enclose_both_poles};
    shift_into_closed_nurbs_domain(rings, surface);
    if let [ring] = rings.as_mut_slice() {
        complement_through_poles(ring, surface, orientation);
    }
    enclose_both_poles(rings, surface, orientation, chord_height);
    align_hole_sheets(rings, surface);
}

pub(super) fn align_hole_sheets(rings: &mut [BoundaryRing], surface: &SurfaceGeometry) {
    let (period_u, period_v) =
        popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(surface);
    let offsets = |period: Option<f64>| -> Vec<f64> {
        match period {
            Some(p) if p > 0.0 => [0.0, -1.0, 1.0, -2.0, 2.0].iter().map(|k| k * p).collect(),
            _ => vec![0.0],
        }
    };
    let (du_offsets, dv_offsets) = (offsets(period_u), offsets(period_v));
    if du_offsets.len() == 1 && dv_offsets.len() == 1 {
        return;
    }
    let Some((outer, holes)) = rings.split_first_mut() else {
        return;
    };
    let outer = std::slice::from_ref(&*outer);
    for hole in holes {
        // A genuine hole has ALL its samples inside the outer ring; a
        // second boundary component that merely shares part of the
        // outer ring's border (a pinch lobe) never has a majority, and
        // must not be dragged onto the outer's sheet (measured: the
        // torus band's upper lobe was shifted onto the lower lobe's
        // equator and coverage failed on every coincident segment).
        let mut best = (hole.uv.len() / 2, 0.0, 0.0);
        for &du in &du_offsets {
            for &dv in &dv_offsets {
                let inside = hole
                    .uv
                    .iter()
                    .filter(|(u, v)| super::point_in_region((u + du, v + dv), outer))
                    .count();
                if inside > best.0 {
                    best = (inside, du, dv);
                }
            }
        }
        let (_, du, dv) = best;
        if du != 0.0 || dv != 0.0 {
            if popkern_core::trace::trace_enabled("boundary-first") {
                eprintln!(
                    "    [bf] hole ring shifted by ({du:.4},{dv:.4}) onto the outer ring's sheet"
                );
            }
            for uv in &mut hole.uv {
                uv.0 += du;
                uv.1 += dv;
            }
        }
    }
}

// On decline, dump the unwrapped ring samples: the jump
// pattern distinguishes genuine winding (monotone drift
// over the ring) from branch-incoherent pcurves (one
// edge's run offset by a full period) — they need
// different fixes.
pub(super) fn dump_ring(ring: &BoundaryRing, dir: &str) {
    if popkern_core::trace::trace_enabled("boundary-first") {
        eprintln!("    [bf] full {dir}-wrap ring — cannot chart");
        for (i, ((u, v), p)) in ring.uv.iter().zip(ring.pos.iter()).enumerate() {
            eprintln!(
                "      [bf] ring[{i}] uv=({u:.4},{v:.4}) pos=({:.4},{:.4},{:.4})",
                p.x, p.y, p.z
            );
        }
    }
}

/// Rings on a NURBS surface CLOSED in a direction (first and last
/// control rows coincide) may come out of the consecutive pcurve walk
/// one or more whole periods off the knot domain (a seamed band read
/// from its seam at the domain end runs [2π, 4π]). The polygon is
/// right, but a NURBS evaluates outside its domain by clamping, so
/// every interior point lands on the domain-end iso-line (measured:
/// jack_in_the_box body 56's extrusion band meshed 84.2 of OCC's
/// 127.5, watertight). The surface has no analytic period to wrap by,
/// so the rings move, together, by the whole domain lengths that put
/// the outer ring's centre inside the domain.
fn shift_into_closed_nurbs_domain(rings: &mut [BoundaryRing], surface: &SurfaceGeometry) {
    let SurfaceGeometry::Nurbs(n) = surface else {
        return;
    };
    let Some(outer) = rings.first() else {
        return;
    };
    let ((u0, u1), (v0, v1)) = n.parameter_range();
    let scale = n
        .control_points
        .iter()
        .fold(0.0_f64, |m, p| {
            m.max(p.x.abs()).max(p.y.abs()).max(p.z.abs())
        })
        .max(1.0);
    let rows_match = |a: usize, b: usize| {
        n.control_points[a].distance_to(&n.control_points[b]) <= 1e-9 * scale
            && (n.weights[a] - n.weights[b]).abs() <= 1e-12 * n.weights[a].abs().max(1.0)
    };
    let closed_u = n.num_u >= 2 && (0..n.num_v).all(|j| rows_match(j, (n.num_u - 1) * n.num_v + j));
    let closed_v =
        n.num_v >= 2 && (0..n.num_u).all(|i| rows_match(i * n.num_v, i * n.num_v + n.num_v - 1));
    let shift = |closed: bool, lo: f64, hi: f64, values: &mut dyn Iterator<Item = f64>| -> f64 {
        let period = hi - lo;
        if !closed || period <= 0.0 {
            return 0.0;
        }
        let (mut a, mut b) = (f64::INFINITY, f64::NEG_INFINITY);
        for x in values {
            a = a.min(x);
            b = b.max(x);
        }
        let eps = 1e-9 * period;
        if a >= lo - eps && b <= hi + eps {
            return 0.0;
        }
        -(((0.5 * (a + b) - lo) / period).floor() * period)
    };
    let du = shift(closed_u, u0, u1, &mut outer.uv.iter().map(|q| q.0));
    let dv = shift(closed_v, v0, v1, &mut outer.uv.iter().map(|q| q.1));
    if du == 0.0 && dv == 0.0 {
        return;
    }
    for ring in rings.iter_mut() {
        for q in &mut ring.uv {
            q.0 += du;
            q.1 += dv;
        }
    }
}

#[cfg(test)]
#[path = "boundary_first_unwrap_tests.rs"]
mod tests;
