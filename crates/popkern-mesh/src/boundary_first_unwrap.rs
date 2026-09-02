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
pub(super) fn align_hole_sheets(rings: &mut [BoundaryRing], surface: &SurfaceGeometry) {
    let (period_u, period_v) = surface.parameter_periodicity();
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
