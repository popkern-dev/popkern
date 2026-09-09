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

/// One cut ring path: window uv points, matched 3-D positions, and
/// per-point canonical edge spans.
type CutPath = (
    Vec<(f64, f64)>,
    Vec<Point3>,
    Vec<(popkern_core::topology::EdgeHandle, f64, f64)>,
);

/// Net u-winding of a consecutively-unwrapped ring, in periods
/// (rounded): +1/-1 for a full wrap, 0 for a chartable ring.
pub(super) fn net_winding(ring: &BoundaryRing, period: f64) -> i64 {
    match (ring.uv.first(), ring.uv.last()) {
        (Some(f), Some(l)) => ((l.0 - f.0) / period).round() as i64,
        _ => 0,
    }
}

/// Collapse short out-and-back STUBS (a dissolved seam's remnant: the
/// ring walks out to a point and straight back within one or two
/// samples). `dissolve_pinches` cannot cut these — both lobes of a
/// legal cut need three points — and a surviving stub plants exact
/// duplicate points that the Delaunay coverage reads as unresolvable
/// zero-distance T-points.
pub(super) fn collapse_stubs(ring: &mut BoundaryRing) {
    loop {
        let n = ring.uv.len();
        if n < 5 {
            return;
        }
        let mut cut: Option<(usize, usize)> = None;
        'scan: for i in 0..n {
            for dj in 2..=3 {
                let j = (i + dj) % n;
                let (pu, pv) = ring.uv[i];
                let (qu, qv) = ring.uv[j];
                if (pu - qu).hypot(pv - qv) < 1.0e-9 {
                    cut = Some((i, dj));
                    break 'scan;
                }
            }
        }
        let Some((i, dj)) = cut else {
            return;
        };
        // Delete points i+1 ..= i+dj (cyclic), keeping one pinch copy.
        let dels: Vec<usize> = (1..=dj).map(|k| (i + k) % n).collect();
        let keep = |k: &usize| !dels.contains(k);
        let idx: Vec<usize> = (0..n).filter(keep).collect();
        ring.uv = idx.iter().map(|&k| ring.uv[k]).collect();
        ring.pos = idx.iter().map(|&k| ring.pos[k]).collect();
        ring.sources = idx.iter().map(|&k| ring.sources[k]).collect();
    }
}

/// Build the slit polygon from two opposite-winding wrapped rings.
/// Returns `None` (caller keeps its decline) when the configuration
/// is outside the certified scope.
pub(super) fn slit_band(
    ring_a: &BoundaryRing,
    ring_b: &BoundaryRing,
    period: f64,
    surface: &SurfaceGeometry,
) -> Option<BoundaryRing> {
    // Candidate azimuths: pick one that each ring's cyclic path
    // crosses exactly once (a wiggly boundary crosses some azimuths
    // three or more times; the single-crossing cut keeps the polygon
    // simple).
    let crossings = |ring: &BoundaryRing, u_star: f64| -> usize {
        let mut count = 0usize;
        let n = ring.uv.len();
        for i in 0..n {
            let ua = ring.uv[i].0;
            let ub = if i + 1 < n {
                ring.uv[i + 1].0
            } else {
                // Closure jump: the ring is closed on the quotient;
                // its unwrapped closure lands a whole winding away.
                ring.uv[0].0
                    + (ring.uv[n - 1].0 - ring.uv[0].0 + 0.5 * period).div_euclid(period) * period
            };
            let k_lo = ((ua.min(ub) - u_star) / period).ceil() as i64;
            let k_hi = ((ua.max(ub) - u_star) / period).floor() as i64;
            if k_hi >= k_lo {
                count += (k_hi - k_lo + 1) as usize;
            }
        }
        count
    };
    // u* must also stand clear of every boundary sample's azimuth:
    // a wall that grazes a sample duplicates it on insertion, and the
    // Delaunay coverage sees an unresolvable zero-distance T-point
    // (measured: the retry rounds exhausted and the face fell back).
    let clear_of_samples = |u_star: f64| -> bool {
        let clearance = 1.0e-3 * period;
        for ring in [ring_a, ring_b] {
            for q in &ring.uv {
                let d = (q.0 - u_star).rem_euclid(period);
                if d < clearance || period - d < clearance {
                    return false;
                }
            }
        }
        true
    };
    let mut u_star = None;
    for k in 0..37 {
        let cand = ring_a.uv[0].0 + period * (0.031 + k as f64 / 37.0);
        if crossings(ring_a, cand) == 1 && crossings(ring_b, cand) == 1 && clear_of_samples(cand) {
            u_star = Some(cand);
            break;
        }
    }
    let u_star = u_star?;
    if popkern_core::trace::trace_enabled("boundary-first") {
        eprintln!(
            "      [bf-slit] u* {u_star:.6} (ring sizes {} / {})",
            ring_a.uv.len(),
            ring_b.uv.len()
        );
    }

    // Cut one ring open at its single crossing of u* (mod period):
    // rotate so the crossing sits between the last and first sample,
    // shift every sample into [u*, u* + p], and interpolate exact
    // points at both walls. `rising` = the path runs from the u* wall
    // to the u*+p wall.
    let cut = |ring: &BoundaryRing, rising: bool| -> Option<CutPath> {
        let n = ring.uv.len();
        // Find the crossing segment in the cyclic walk.
        let mut cross_at = None;
        for i in 0..n {
            let ua = ring.uv[i].0;
            let ub = ring.uv[(i + 1) % n].0
                + if i + 1 == n {
                    ((ring.uv[n - 1].0 - ring.uv[0].0 + 0.5 * period) / period).floor() * period
                } else {
                    0.0
                };
            let k = ((ua - u_star) / period).floor();
            let wall = u_star + (k + 1.0) * period;
            if (ua < wall && ub >= wall) || (ua >= wall && ub < wall) {
                cross_at = Some(i);
                break;
            }
        }
        let start = (cross_at? + 1) % n;
        // Walk the whole ring from `start`, re-unwrapping
        // consecutively, then shift into the window.
        let mut uv: Vec<(f64, f64)> = Vec::with_capacity(n + 2);
        let mut pos: Vec<Point3> = Vec::with_capacity(n + 2);
        let mut sources = Vec::with_capacity(n + 2);
        let mut prev_u: Option<f64> = None;
        for k in 0..n {
            let i = (start + k) % n;
            let (mut u, v) = ring.uv[i];
            if let Some(p) = prev_u {
                while u - p > 0.5 * period {
                    u -= period;
                }
                while u - p < -0.5 * period {
                    u += period;
                }
            }
            prev_u = Some(u);
            uv.push((u, v));
            pos.push(ring.pos[i]);
            sources.push(ring.sources[i]);
        }
        // Shift so the path occupies [u*, u* + p].
        let lo = uv.iter().map(|q| q.0).fold(f64::MAX, f64::min);
        let shift = ((lo - u_star) / period).floor() * period;
        for q in &mut uv {
            q.0 -= shift;
        }
        // Everything must stay inside the window (small excursions
        // past a wall would cross the seam ladders).
        let hi_wall = u_star + period;
        if uv
            .iter()
            .any(|q| q.0 < u_star - 0.02 * period || q.0 > hi_wall + 0.02 * period)
        {
            return None;
        }
        // Exact wall points. The single cut segment runs last ->
        // first-across-the-closure and crosses exactly one wall; the
        // crossing's v serves BOTH path ends (they are the same 3-D
        // point seen from the two walls).
        let first = uv[0];
        let last = *uv.last()?;
        let (start_wall, end_wall) = if last.0 >= first.0 {
            (u_star, hi_wall)
        } else {
            (hi_wall, u_star)
        };
        let cut_b_u = first.0 + if last.0 >= first.0 { period } else { -period };
        let t = if (cut_b_u - last.0).abs() < 1e-30 {
            0.0
        } else {
            ((end_wall - last.0) / (cut_b_u - last.0)).clamp(0.0, 1.0)
        };
        let v_cross = last.1 + (first.1 - last.1) * t;
        let end_pt = (end_wall, v_cross);
        let start_pt = (start_wall, v_cross);
        let eval = |q: (f64, f64)| surface.point_at(q.0, q.1);
        let src0 = sources[0];
        let src_end = *sources.last()?;
        uv.insert(0, start_pt);
        pos.insert(0, eval(start_pt));
        sources.insert(0, (src0.0, src0.1, src0.1));
        uv.push(end_pt);
        pos.push(eval(end_pt));
        sources.push((src_end.0, src_end.2, src_end.2));
        // Orient the path in the requested direction.
        let runs_rising = uv.last()?.0 >= uv[0].0;
        if runs_rising != rising {
            uv.reverse();
            pos.reverse();
            sources.reverse();
        }
        Some((uv, pos, sources))
    };

    let (a_uv, a_pos, a_src) = cut(ring_a, true)?;
    let (b_uv, b_pos, b_src) = cut(ring_b, false)?;

    // Seam ladders: shared v stations on both walls so the two sides
    // of the cut weld. Station count from the v spans.
    let ladder = |from: (f64, f64),
                  to: (f64, f64),
                  src: (popkern_core::topology::EdgeHandle, f64, f64)|
     -> CutPath {
        let steps = ((from.1 - to.1).abs() / (0.02 * period).max(1e-9))
            .ceil()
            .clamp(1.0, 64.0) as usize;
        let mut uv = Vec::new();
        let mut pos = Vec::new();
        let mut sources = Vec::new();
        for k in 1..steps {
            let v = from.1 + (to.1 - from.1) * k as f64 / steps as f64;
            uv.push((from.0, v));
            pos.push(surface.point_at(from.0, v));
            sources.push((src.0, src.1, src.1));
        }
        (uv, pos, sources)
    };
    let a_end = *a_uv.last()?;
    let b_start = b_uv[0];
    let b_end = *b_uv.last()?;
    let a_start = a_uv[0];
    let (l1_uv, l1_pos, l1_src) = ladder(a_end, b_start, *a_src.last()?);
    let (l2_uv, l2_pos, l2_src) = ladder(b_end, a_start, *b_src.last()?);

    let mut uv = a_uv;
    let mut pos = a_pos;
    let mut sources = a_src;
    uv.extend(l1_uv);
    pos.extend(l1_pos);
    sources.extend(l1_src);
    uv.extend(b_uv);
    pos.extend(b_pos);
    sources.extend(b_src);
    uv.extend(l2_uv);
    pos.extend(l2_pos);
    sources.extend(l2_src);

    // Drop coincident consecutive points (closure pair included) —
    // zero-length segments are unresolvable T-points to the coverage
    // check.
    let mut keep = vec![true; uv.len()];
    for i in 0..uv.len() {
        let j = (i + 1) % uv.len();
        let d2 = (uv[i].0 - uv[j].0).powi(2) + (uv[i].1 - uv[j].1).powi(2);
        if d2 < 1.0e-18 {
            keep[j] = false;
        }
    }
    let filter3 = |v: Vec<(f64, f64)>| -> Vec<(f64, f64)> {
        v.into_iter()
            .zip(keep.iter())
            .filter(|&(_, &k)| k)
            .map(|(x, _)| x)
            .collect()
    };
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
    let uv = filter3(uv);
    if uv.len() < 4 {
        return None;
    }
    Some(BoundaryRing { uv, pos, sources })
}

/// Pair and join a face's wrapped rings (moved out of
/// `boundary_first.rs` under its size ratchet): exactly two
/// |winding| = 1 rings slit into one chartable polygon; anything else
/// keeps the honest decline with the ring dump.
pub(super) fn join_wrapped_rings(
    face_handle: popkern_core::topology::FaceHandle,
    wrapped_rings: &[BoundaryRing],
    p: f64,
    surface: &SurfaceGeometry,
) -> Option<BoundaryRing> {
    // Stored walk direction decides the SIGN; the slit orients both
    // paths explicitly, so only |winding| == 1 matters.
    if wrapped_rings.len() == 2
        && net_winding(&wrapped_rings[0], p).abs() == 1
        && net_winding(&wrapped_rings[1], p).abs() == 1
        && let Some(joined) = slit_band(&wrapped_rings[0], &wrapped_rings[1], p, surface)
    {
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!(
                "    [bf] face {face_handle:?} annular-band slit: joined {} pts",
                joined.uv.len()
            );
            let n = joined.uv.len();
            for i in (0..3).chain(n.saturating_sub(6)..n) {
                eprintln!(
                    "      [bf-slit] pt[{i}] uv=({:.4},{:.4})",
                    joined.uv[i].0, joined.uv[i].1
                );
            }
        }
        return Some(joined);
    }
    if popkern_core::trace::trace_enabled("boundary-first") {
        let winds: Vec<i64> = wrapped_rings.iter().map(|r| net_winding(r, p)).collect();
        eprintln!(
            "    [bf] face {face_handle:?} u-wrap decline ({} wrapped ring(s), windings {winds:?})",
            wrapped_rings.len()
        );
    }
    for r in wrapped_rings {
        super::boundary_first_unwrap::dump_ring(r, "u");
    }
    None
}

/// Dissolve a wrap-failed ring's slits, pinches, and stubs, then
/// stash the surviving lobes for pairing. The ring arrives
/// consecutively unwrapped, so slit dissolution is UNAMBIGUOUS (the
/// periodic scoping in `dissolve_slits` guards raw-branch seam
/// coincidence, which the unwrap already resolved), and pinches must
/// dissolve here because the original decline returned before the
/// pinch stage — a dissolved-seam stub left in place collapses the
/// slit polygon's coverage with zero-width T-points (measured: an
/// out-and-back spike at u = 5π/2 on the m5 tip band's bottom ring).
pub(super) fn stash_wrapped_ring(ring: BoundaryRing, wrapped_rings: &mut Vec<BoundaryRing>) {
    for piece in super::boundary_first_slit::dissolve_slits(ring, (None, None)) {
        for mut lobe in super::boundary_first_slit::dissolve_pinches(piece) {
            collapse_stubs(&mut lobe);
            if lobe.uv.len() >= 3 {
                wrapped_rings.push(lobe);
            }
        }
    }
}
