//! Pole runs: a sphere ring that passes through BOTH chart poles
//! bounds two regions, and the face flag picks which (PK-646).
//!
//! A sphere face bounded by two pole-to-pole curves (a lune — the
//! sphere minus a thin wedge, written as one face with no seam) has a
//! zero-winding chart ring: down one meridian, up the other. The ring
//! encloses the strip BETWEEN the meridians; the face may be that
//! strip or its complement, and only the flag says which. The kernel
//! flag says the loop runs counter-clockwise in the chart (the face
//! lies left of the walk), so a ring whose chart winding disagrees
//! bounds the complement (measured: jack_in_the_box body 33, a
//! 2.85° wedge cut from a ball, meshed the 1.82-volume wedge where
//! OpenCASCADE's solid holds 90.07).
//!
//! At a pole every u is the same 3-D point, so the chart may run along
//! the pole row for any whole number of periods at no 3-D cost. The
//! complement is charted by making each pole visit an explicit pole
//! row and moving the section between the two visits one period: the
//! strip `[a, b]` becomes `[b − p, a]`. Pole-row samples carry the
//! pole's 3-D point and a zero-length span of the edge that leaves the
//! pole, exactly as the pole cap's return run does.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::BoundaryRing;

/// Re-chart `ring` as the complement of the strip it encloses when it
/// visits each sphere pole once and runs against the face flag.
/// Leaves the ring unchanged (returns `false`) otherwise.
pub(in crate::boundary_first) fn complement_through_poles(
    ring: &mut BoundaryRing,
    surface: &SurfaceGeometry,
    orientation: bool,
) -> bool {
    let SurfaceGeometry::Sphere { .. } = surface else {
        return false;
    };
    let n = ring.uv.len();
    if n < 4 {
        return false;
    }
    let area = signed_area(&ring.uv);
    if area == 0.0 || (area > 0.0) == orientation {
        return false;
    }
    let period = std::f64::consts::TAU;
    let poles = [
        (surface.point_at(0.0, 0.0), 0.0),
        (
            surface.point_at(0.0, std::f64::consts::PI),
            std::f64::consts::PI,
        ),
    ];
    let scale = ring
        .pos
        .iter()
        .zip(ring.pos.iter().cycle().skip(1))
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max)
        .max(1e-12);
    let pole_of = |p: &Point3| {
        poles
            .iter()
            .position(|(q, _)| q.distance_to(p) <= 1e-6 * scale)
    };
    let tags: Vec<Option<usize>> = ring.pos.iter().map(pole_of).collect();
    // Start just after a pole run so runs never straddle the wrap.
    let Some(start) = (0..n).find(|&i| tags[i].is_none() && tags[(i + n - 1) % n].is_some()) else {
        return false;
    };
    let order: Vec<usize> = (0..n).map(|k| (start + k) % n).collect();
    // Sections between pole runs, and the runs themselves.
    let mut runs: Vec<(usize, Vec<usize>)> = Vec::new(); // (pole, indices)
    let mut sections: Vec<Vec<usize>> = vec![Vec::new()];
    for &i in &order {
        match tags[i] {
            Some(pole) => match runs.last_mut() {
                Some((p, idx)) if *p == pole && sections.last().is_some_and(Vec::is_empty) => {
                    idx.push(i)
                }
                _ => {
                    runs.push((pole, vec![i]));
                    sections.push(Vec::new());
                }
            },
            None => {
                if let Some(last) = sections.last_mut() {
                    last.push(i);
                }
            }
        }
    }
    // The ring ends on a run (it started just after one), so the
    // trailing section is empty: exactly two runs, one per pole.
    if runs.len() != 2 || runs[0].0 == runs[1].0 || sections.len() != 3 {
        return false;
    }
    if !sections[2].is_empty() || sections[0].is_empty() || sections[1].is_empty() {
        return false;
    }
    let mut best: Option<(f64, BoundaryRing)> = None;
    for shift in [-period, period] {
        let candidate = assemble(ring, &sections, &runs, &poles, shift, period);
        let a = signed_area(&candidate.uv);
        // The complement has the flag's winding and lies within one
        // chart period; the other shift overlaps the strip itself.
        let bound = period * std::f64::consts::PI;
        if (a > 0.0) == orientation
            && a.abs() < bound
            && best.as_ref().is_none_or(|b| a.abs() < b.0)
        {
            best = Some((a.abs(), candidate));
        }
    }
    match best {
        Some((_, re_charted)) => {
            *ring = re_charted;
            true
        }
        None => false,
    }
}

/// The ring rebuilt as: section 0 as walked, the first pole run as a
/// pole row, section 1 moved by `shift`, the second pole run as a pole
/// row back to section 0's start.
fn assemble(
    ring: &BoundaryRing,
    sections: &[Vec<usize>],
    runs: &[(usize, Vec<usize>)],
    poles: &[(Point3, f64); 2],
    shift: f64,
    period: f64,
) -> BoundaryRing {
    let mut out = BoundaryRing {
        uv: Vec::new(),
        pos: Vec::new(),
        sources: Vec::new(),
    };
    let push_section = |out: &mut BoundaryRing, section: &[usize], du: f64| {
        for &i in section {
            out.uv.push((ring.uv[i].0 + du, ring.uv[i].1));
            out.pos.push(ring.pos[i]);
            out.sources.push(ring.sources[i]);
        }
    };
    let push_pole_row =
        |out: &mut BoundaryRing, run: &[usize], pole: usize, (from, to): (f64, f64)| {
            let (point, v) = poles[pole];
            // The run's last sample starts the segment that leaves the
            // pole; its span begins at the pole vertex.
            let leaving = ring.sources[run[run.len() - 1]];
            let span = (leaving.0, leaving.1, leaving.1);
            let steps = ((to - from).abs() / (period / 128.0)).ceil().max(1.0) as usize;
            for k in 0..=steps {
                out.uv
                    .push((from + (to - from) * k as f64 / steps as f64, v));
                out.pos.push(point);
                out.sources.push(if k == steps { leaving } else { span });
            }
        };
    let u_of = |i: usize, du: f64| ring.uv[i].0 + du;
    let (s0, s1) = (&sections[0], &sections[1]);
    push_section(&mut out, s0, 0.0);
    push_pole_row(
        &mut out,
        &runs[0].1,
        runs[0].0,
        (u_of(s0[s0.len() - 1], 0.0), u_of(s1[0], shift)),
    );
    push_section(&mut out, s1, shift);
    push_pole_row(
        &mut out,
        &runs[1].1,
        runs[1].0,
        (u_of(s1[s1.len() - 1], shift), u_of(s0[0], 0.0)),
    );
    out
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

/// A sphere face whose outer ring runs against the face flag and
/// touches neither pole is the sphere OUTSIDE that ring — a ball with
/// a stem, a spherical shell with a hole — and contains both poles,
/// which no polygon in the chart encloses (measured: the rest of a
/// radius-10 sphere beyond a cap clear of both poles meshed the cap,
/// 904 of 3284). Its chart is the whole period box: a seam ladder
/// down the meridian `u0` (met from both sides, the same v stations
/// on each wall), the two pole rows, and the rings as holes. `u0`
/// sits in the middle of the azimuth band the rings leave free.
/// Leaves the rings unchanged (returns `false`) otherwise.
pub(in crate::boundary_first) fn enclose_both_poles(
    rings: &mut Vec<BoundaryRing>,
    surface: &SurfaceGeometry,
    orientation: bool,
    chord_height: f64,
) -> bool {
    use std::f64::consts::{PI, TAU};
    let SurfaceGeometry::Sphere { .. } = surface else {
        return false;
    };
    let Some(outer) = rings.first() else {
        return false;
    };
    let area = signed_area(&outer.uv);
    if outer.uv.len() < 3 || area == 0.0 || (area > 0.0) == orientation {
        return false;
    }
    let north = surface.point_at(0.0, 0.0);
    let south = surface.point_at(0.0, PI);
    let longest_step = outer
        .pos
        .iter()
        .zip(outer.pos.iter().cycle().skip(1))
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max);
    let near_pole = |p: &Point3| {
        p.distance_to(&north) <= 1e-6 * longest_step.max(1e-12)
            || p.distance_to(&south) <= 1e-6 * longest_step.max(1e-12)
    };
    if rings.iter().any(|r| r.pos.iter().any(near_pole)) {
        return false; // pole-touching rings are complement_through_poles' class
    }
    let (lo, hi) = rings
        .iter()
        .flat_map(|r| r.uv.iter().map(|q| q.0))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), u| {
            (a.min(u), b.max(u))
        });
    if !hi.is_finite() || !lo.is_finite() || hi - lo >= TAU {
        return false;
    }
    let u0 = hi + 0.5 * (TAU - (hi - lo)) - TAU;
    let u1 = u0 + TAU;
    let source = outer.sources[0];
    let span = (source.0, source.1, source.1);
    let ladder = super::wall_stations(surface, u0, (PI, 0.0), chord_height, longest_step);
    let row = |from: f64, to: f64| -> Vec<f64> {
        let steps = ((to - from).abs() / (TAU / 128.0)).ceil().max(1.0) as usize;
        (0..steps)
            .map(|k| from + (to - from) * k as f64 / steps as f64)
            .collect()
    };
    let mut boxed = BoundaryRing {
        uv: Vec::new(),
        pos: Vec::new(),
        sources: Vec::new(),
    };
    let mut push = |u: f64, v: f64| {
        boxed.uv.push((u, v));
        boxed.pos.push(surface.point_at(u, v));
        boxed.sources.push(span);
    };
    // Counter-clockwise: the south row toward +u, up the u1 wall, the
    // north row back, down the u0 wall.
    for u in row(u0, u1) {
        push(u, PI);
    }
    push(u1, PI);
    for &v in &ladder {
        push(u1, v);
    }
    for u in row(u1, u0) {
        push(u, 0.0);
    }
    push(u0, 0.0);
    for &v in ladder.iter().rev() {
        push(u0, v);
    }
    rings.insert(0, boxed);
    true
}
