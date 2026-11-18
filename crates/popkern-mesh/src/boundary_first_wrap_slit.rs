//! Seam slitting for FULL-WRAP boundary rings on singly-periodic
//! charts (PK-468).
//!
//! A face whose boundary rings each wind a full chart period — an
//! annular band on a cone or cylinder whose seam the merge dissolved —
//! cannot be charted as simple polygons by the unwrap, and the
//! boundary-first tessellator used to decline it to the uniform-grid
//! fallback. A fallback face inside a boundary-first body tears along
//! its ENTIRE rim (its boundary sampling matches no neighbour):
//! measured on the m5 countersunk tip, the cone band fell back and the
//! served, geometrically-correct body meshed with 154 open segments.
//!
//! The slit cuts the band at one azimuth: each wrapped ring is rotated
//! and unwrapped into the window [u*, u* + p], exact crossing points
//! are interpolated at both window walls, and the two ring paths are
//! joined through two seam ladders (the same v stations on both walls,
//! whose 3D positions agree to floating-point periodicity and weld).
//! The result is ONE simple polygon the standard Delaunay covers.
//!
//! Scope: exactly two wrapped rings with opposite net winding, each
//! crossing the chosen azimuth exactly once, all samples within the
//! window. Anything else keeps the honest decline.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::BoundaryRing;

#[path = "boundary_first_wrap_slit_walls.rs"]
mod walls;
use walls::wall_crossings;
pub(super) use walls::{collapse_stubs, net_winding};

/// One cut ring path: window uv points, matched 3-D positions, and
/// per-point canonical edge spans.
type CutPath = (
    Vec<(f64, f64)>,
    Vec<Point3>,
    Vec<(popkern_core::topology::EdgeHandle, f64, f64)>,
);

/// Build the slit polygon from wrapped rings: cut every ring at one
/// sample-clear azimuth into wall-to-wall pieces, then walk the piece
/// graph — at each wall exit, continue along the wall (a seam ladder)
/// to the v-nearest unused piece endpoint on the same wall, reversing
/// pieces as needed — until the cycle closes. Handles the annular
/// band (2 pieces), the two-start helical band (one ring, 2 pieces),
/// and multi-turn helical bands (2·turns pieces). Returns `None`
/// (caller keeps its decline) when no clear azimuth exists, a piece
/// leaves the window, or the walk cannot consume every piece in one
/// cycle.
pub(super) fn slit_band(
    (wrapped, others): (&[BoundaryRing], &[BoundaryRing]),
    period: f64,
    surface: &SurfaceGeometry,
    transposed: bool,
    (chord_height, azimuth): (f64, Option<f64>),
) -> Option<Vec<BoundaryRing>> {
    // In a transposed chart (rings wrapping in v) the first coordinate
    // is the surface's v.
    let point_at = |a: f64, b: f64| {
        if transposed {
            surface.point_at(b, a)
        } else {
            surface.point_at(a, b)
        }
    };
    let ring0 = wrapped.first()?;
    // A ring's wall crossings at `u_star` (a ring running along a wall
    // counts as uncuttable: more than any ring could cross).
    let crossings_of = |ring: &BoundaryRing, u_star: f64| -> usize {
        wall_crossings(ring, u_star, period).map_or(usize::MAX / 8, |c| c.len())
    };
    // The face's other rings: the slit runs wall to wall, so it may cross
    // none, nor graze their samples (a dome band's slit fell 7e-4 from a
    // loop). A transposed band's rings are transposed; others are not read.
    let others: &[BoundaryRing] = if transposed { &[] } else { others };
    let clear_of_samples = |u_star: f64| -> bool {
        let clearance = 1.0e-5 * period;
        if others.iter().any(|r| crossings_of(r, u_star) > 0) {
            return false;
        }
        for ring in wrapped.iter().chain(others) {
            for q in &ring.uv {
                let d = popkern_math::periodic::offset_mod(q.0, u_star, period);
                if d < clearance || period - d < clearance {
                    return false;
                }
            }
        }
        true
    };
    // Fewest total crossings wins (>= 2, even), then the azimuth farthest
    // from the other rings' samples.
    use popkern_core::geometry::chart::atlas::along_offset as offset;
    let other_clearance = |u_star: f64| -> f64 {
        let dist = |q: &(f64, f64)| offset(Some(period), u_star, q.0, 1.0);
        let near = |d: f64| d.min(period - d);
        others
            .iter()
            .flat_map(|r| r.uv.iter())
            .map(|q| near(dist(q)))
            .fold(f64::INFINITY, f64::min)
    };
    let mut best: Option<(usize, f64, f64)> = None;
    for k in 0..37 {
        let cand = ring0.uv[0].0 + period * (0.031 + k as f64 / 37.0);
        if !clear_of_samples(cand) {
            continue;
        }
        let total: usize = wrapped.iter().map(|r| crossings_of(r, cand)).sum();
        let clearance = other_clearance(cand);
        let better = |(b, _, bc): (usize, f64, f64)| total < b || (total == b && clearance > bc);
        if total >= 2 && total.is_multiple_of(2) && best.is_none_or(better) {
            best = Some((total, cand, clearance));
        }
    }
    // An atlas band's meridian is chosen by the atlas (its real-edge
    // crossings are samples on it): cut there or not at all.
    if let Some(az) = azimuth {
        let total: usize = wrapped.iter().map(|r| crossings_of(r, az)).sum();
        let clear = others.iter().all(|r| crossings_of(r, az) == 0);
        best = (clear && total >= 2 && total.is_multiple_of(2)).then_some((total, az, 0.0));
    }
    let Some((_, u_star, _)) = best else {
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!("      [bf-slit] no clear even-crossing azimuth found");
        }
        return None;
    };
    let hi_wall = u_star + period;
    // Both walls are one meridian in space: evaluated at u* alike, so the
    // two copies of a wall point agree to the bit (pieces weld by bits).
    let wall_point = |u: f64, v: f64| {
        point_at(
            if (u - hi_wall).abs() < 1e-9 {
                u_star
            } else {
                u
            },
            v,
        )
    };

    // Cut every ring into pieces between successive wall crossings.
    let mut pieces: Vec<CutPath> = Vec::new();
    for ring in wrapped {
        let n = ring.uv.len();
        let mut cross = wall_crossings(ring, u_star, period)?;
        if cross.is_empty() {
            return None;
        }
        cross.sort_by(|a, b| a.segment.cmp(&b.segment).then(a.t.total_cmp(&b.t)));
        let m = cross.len();
        for c in 0..m {
            let (start, end) = (cross[c], cross[(c + 1) % m]);
            let (from, v_start) = (start.segment, start.v);
            let (to, v_end) = (end.segment, end.v);
            // Samples strictly after `from`'s segment up to and
            // including `to`'s segment start.
            let mut uv: Vec<(f64, f64)> = Vec::new();
            let mut pos: Vec<Point3> = Vec::new();
            let mut sources = Vec::new();
            let mut prev_u: Option<f64> = None;
            let mut k = (from + 1) % n;
            loop {
                // A wall sample ending the piece is its wall point.
                if k == to && end.at_sample {
                    break;
                }
                let (mut u, v) = ring.uv[k];
                if let Some(pv) = prev_u {
                    while u - pv > 0.5 * period {
                        u -= period;
                    }
                    while u - pv < -0.5 * period {
                        u += period;
                    }
                }
                prev_u = Some(u);
                uv.push((u, v));
                pos.push(ring.pos[k]);
                sources.push(ring.sources[k]);
                if k == to {
                    break;
                }
                k = (k + 1) % n;
                if uv.len() > n + 2 {
                    return None;
                }
            }
            if uv.is_empty() {
                return None;
            }
            // Shift into the window; a piece between successive
            // crossings spans at most one period.
            let lo = uv.iter().map(|q| q.0).fold(f64::MAX, f64::min);
            let shift = ((lo - u_star) / period).floor() * period;
            for q in &mut uv {
                q.0 -= shift;
            }
            if uv
                .iter()
                .any(|q| q.0 < u_star - 0.02 * period || q.0 > hi_wall + 0.02 * period)
            {
                if popkern_core::trace::trace_enabled("boundary-first") {
                    let (lo2, hi2) = uv
                        .iter()
                        .fold((f64::MAX, f64::MIN), |(a, b), q| (a.min(q.0), b.max(q.0)));
                    eprintln!(
                        "      [bf-slit] piece leaves window: u [{lo2:.3},{hi2:.3}] vs [{u_star:.3},{hi_wall:.3}]"
                    );
                }
                return None;
            }
            let near_wall = |u: f64| {
                if (u - u_star).abs() <= (u - hi_wall).abs() {
                    u_star
                } else {
                    hi_wall
                }
            };
            let sw = near_wall(uv[0].0);
            let ew = near_wall(uv.last()?.0);
            // Same-wall (u-turn) pieces are legal: a boundary that
            // dips past a wall and returns contributes a chunk whose
            // both ends sit on one wall; the graph walk traverses it
            // and continues along the same wall.
            // A wall sample is its own wall point (a canonical sample the
            // neighbour face shares); else the wall point is interpolated.
            let src0 = sources[0];
            let src_end = *sources.last()?;
            uv.insert(0, (sw, v_start));
            let (p0, s0) = match start.at_sample {
                true => (ring.pos[from], ring.sources[from]),
                false => (wall_point(sw, v_start), (src0.0, src0.1, src0.1)),
            };
            pos.insert(0, p0);
            sources.insert(0, s0);
            uv.push((ew, v_end));
            let (p1, s1) = match end.at_sample {
                true => (ring.pos[to], ring.sources[to]),
                false => (wall_point(ew, v_end), (src_end.0, src_end.2, src_end.2)),
            };
            pos.push(p1);
            sources.push(s1);
            pieces.push((uv, pos, sources));
        }
    }
    if pieces.len() < 2 || !pieces.len().is_multiple_of(2) {
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!("      [bf-slit] {} piece(s) (need even >= 2)", pieces.len());
        }
        return None;
    }

    // Walk the piece graph with DETERMINISTIC wall pairing: at each
    // wall, sort every piece endpoint by v; consecutive pairs
    // (lowest two, next two, ...) bound the interior wall intervals —
    // below the lowest crossing lies the exterior. A greedy
    // nearest-endpoint walk paired across intervals and the ladders
    // overlapped on the wall (measured: the two-start shank's slit
    // polygon self-overlapped and coverage failed).
    let n_pieces = pieces.len();
    #[derive(Clone, Copy, PartialEq, Eq, Debug)]
    struct EndRef {
        piece: usize,
        at_start: bool,
    }
    let end_uv = |e: EndRef| -> (f64, f64) {
        let p = &pieces[e.piece];
        if e.at_start {
            p.0[0]
        } else {
            *p.0.last().unwrap_or(&p.0[0])
        }
    };
    let mut lo_wall_ends: Vec<EndRef> = Vec::new();
    let mut hi_wall_ends: Vec<EndRef> = Vec::new();
    for pi in 0..pieces.len() {
        for at_start in [true, false] {
            let e = EndRef {
                piece: pi,
                at_start,
            };
            let (u, _) = end_uv(e);
            if (u - u_star).abs() < 1e-9 {
                lo_wall_ends.push(e);
            } else if (u - hi_wall).abs() < 1e-9 {
                hi_wall_ends.push(e);
            } else {
                return None;
            }
        }
    }
    let mut partner: popkern_math::hash::FastHashMap<(usize, bool), EndRef> =
        popkern_math::hash::FastHashMap::default();
    for ends in [&mut lo_wall_ends, &mut hi_wall_ends] {
        if !ends.len().is_multiple_of(2) {
            return None;
        }
        ends.sort_by(|a, b| end_uv(*a).1.total_cmp(&end_uv(*b).1));
        // Same-crossing entries (a shared crossing between two
        // adjacent pieces of one ring — exactly equal v) MUST pair
        // with each other; an exact tie left to the alternating
        // chunking pairs across cycles and the cycles cross
        // (measured: three coincident wall copies on the two-start
        // shank and the coverage saw unresolvable duplicates).
        let mut taken = vec![false; ends.len()];
        for i in 0..ends.len() {
            if taken[i] {
                continue;
            }
            if let Some(jj) = (i + 1..ends.len())
                .find(|&j| !taken[j] && (end_uv(ends[j]).1 - end_uv(ends[i]).1).abs() < 1e-9)
            {
                taken[i] = true;
                taken[jj] = true;
                partner.insert((ends[i].piece, ends[i].at_start), ends[jj]);
                partner.insert((ends[jj].piece, ends[jj].at_start), ends[i]);
            }
        }
        let rest: Vec<EndRef> = ends
            .iter()
            .zip(taken.iter())
            .filter(|&(_, &t)| !t)
            .map(|(e, _)| *e)
            .collect();
        if !rest.len().is_multiple_of(2) {
            return None;
        }
        for pair in rest.chunks(2) {
            partner.insert((pair[0].piece, pair[0].at_start), pair[1]);
            partner.insert((pair[1].piece, pair[1].at_start), pair[0]);
        }
    }
    let mut used = vec![false; n_pieces];
    let push_piece =
        |target: &mut (Vec<(f64, f64)>, Vec<Point3>, Vec<_>), piece: &CutPath, forward: bool| {
            if forward {
                target.0.extend(piece.0.iter().copied());
                target.1.extend(piece.1.iter().copied());
                target.2.extend(piece.2.iter().copied());
            } else {
                target.0.extend(piece.0.iter().rev().copied());
                target.1.extend(piece.1.iter().rev().copied());
                target.2.extend(piece.2.iter().rev().copied());
            }
        };
    let ladder_into = |bundle: &mut (Vec<(f64, f64)>, Vec<Point3>, Vec<_>),
                       from: (f64, f64),
                       to: (f64, f64),
                       src: (popkern_core::topology::EdgeHandle, f64, f64)| {
        let steps = ((from.1 - to.1).abs() / (0.02 * period).max(1e-9))
            .ceil()
            .clamp(1.0, 64.0) as usize;
        let mut vs: Vec<f64> = (0..=steps)
            .map(|k| from.1 + (to.1 - from.1) * k as f64 / steps as f64)
            .collect();
        // A wall along a CURVED iso-line (a torus band's major circle in
        // a transposed chart) also refines until each chord sits within
        // a quarter of the chord height (PK-641: elbows sagged 0.025 at
        // 0.005): its stations are the tensor lattice's rows across the
        // band, whose cell diagonals sag up to twice a side's and whose
        // over-chord cells can only split the rails (a sphere band on the
        // atlas kept 7 rows at the chord height and sagged 5.7e-3 while
        // its rails split every round).
        let wall_chord = 0.25 * chord_height;
        for _ in 0..6 {
            let mut next = vec![vs[0]];
            for w in vs.windows(2) {
                let (a, b) = (wall_point(from.0, w[0]), wall_point(from.0, w[1]));
                let m = wall_point(from.0, 0.5 * (w[0] + w[1]));
                let mid = Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
                if m.distance_to(&mid) > wall_chord && next.len() < 512 {
                    next.push(0.5 * (w[0] + w[1]));
                }
                next.push(w[1]);
            }
            let split = next.len() > vs.len();
            vs = next;
            if !split {
                break;
            }
        }
        for &v in &vs[1..vs.len() - 1] {
            bundle.0.push((from.0, v));
            bundle.1.push(wall_point(from.0, v));
            bundle.2.push((src.0, src.1, src.1));
        }
    };
    // Every closed cycle is its own window polygon: a complex region
    // (a shank wall minus groove strips) legitimately decomposes into
    // several disjoint pieces under the cut; even-odd containment
    // sorts interior from holes downstream.
    let mut rings_out: Vec<CutPath> = Vec::new();
    while let Some(seed) = used.iter().position(|&u| !u) {
        used[seed] = true;
        let mut cycle: (Vec<(f64, f64)>, Vec<Point3>, Vec<_>) =
            (Vec::new(), Vec::new(), Vec::new());
        push_piece(&mut cycle, &pieces[seed], true);
        let start_end = EndRef {
            piece: seed,
            at_start: true,
        };
        let mut cur_exit = EndRef {
            piece: seed,
            at_start: false,
        };
        let mut closed = false;
        for _ in 0..n_pieces {
            let next_entry = *partner.get(&(cur_exit.piece, cur_exit.at_start))?;
            let src = *cycle.2.last()?;
            ladder_into(&mut cycle, end_uv(cur_exit), end_uv(next_entry), src);
            if next_entry == start_end {
                closed = true;
                break;
            }
            if used[next_entry.piece] {
                return None; // pairing inconsistency
            }
            used[next_entry.piece] = true;
            push_piece(&mut cycle, &pieces[next_entry.piece], next_entry.at_start);
            cur_exit = EndRef {
                piece: next_entry.piece,
                at_start: !next_entry.at_start,
            };
        }
        if !closed {
            return None;
        }
        rings_out.push(cycle);
    }
    let mut out: Vec<BoundaryRing> = Vec::new();
    for (uv, pos, sources) in rings_out {
        let n_all = uv.len();
        let mut keep = vec![true; n_all.max(1)];
        for i in 0..n_all {
            let j = (i + 1) % n_all;
            let d2 = (uv[i].0 - uv[j].0).powi(2) + (uv[i].1 - uv[j].1).powi(2);
            if d2 < 1.0e-18 {
                keep[j] = false;
            }
        }
        let uv: Vec<(f64, f64)> = uv
            .into_iter()
            .zip(keep.iter())
            .filter(|&(_, &k)| k)
            .map(|(x, _)| x)
            .collect();
        let pos: Vec<Point3> = pos
            .into_iter()
            .zip(keep.iter())
            .filter(|&(_, &k)| k)
            .map(|(x, _)| x)
            .collect();
        let sources = sources
            .into_iter()
            .zip(keep.iter())
            .filter(|&(_, &k)| k)
            .map(|(x, _)| x)
            .collect();
        if uv.len() >= 4 {
            out.push(BoundaryRing { uv, pos, sources });
        }
    }
    if out.is_empty() {
        return None;
    }
    Some(out)
}
