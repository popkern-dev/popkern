//! The atlas band's cut meridian (chart-service.md §6a, PK-1267).
//!
//! A band piece whose rings wind the along period is charted by the
//! wrap slit, which cuts it at one meridian and interpolates its wall
//! points in the chart. A wall point on a row is on the row (an
//! iso-line), but one on a real edge lies off its curve, and the
//! neighbour face, which shares the edge's canonical samples, never sees
//! it: a T-junction, open in the welded mesh (a dome's band slit crossed
//! its cut circle and left four edges open). So the atlas chooses the
//! meridian itself and makes every real-edge crossing a canonical sample,
//! exactly as it does at its rows: a crossing between two samples is
//! requested at its curve fraction (bisected on the edge in space), and
//! one within the encroachment dedup of a sample snaps that sample's
//! along coordinate onto the meridian (its position is kept). The slit
//! then cuts at that meridian, at those samples.
//!
//! The meridian is chosen from a fixed list of azimuths, never from the
//! samples, so the choice holds as requested samples land on it: the
//! first, by fewest crossings, crossing no ring that does not wind and
//! clear of every ring's turning points in the along coordinate (where a
//! crossing is ill-conditioned and refinement flips it); the clearest
//! valid one when none is clear.

use popkern_core::geometry::chart::atlas::along_offset;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle};

use super::{BoundaryRing, Cut, SplitRequest};

/// Candidate azimuths per period.
const CANDIDATES: usize = 24;

/// The meridian's outcome for a band piece.
pub(super) enum Meridian {
    /// No ring winds the period: the band needs no slit.
    Unneeded,
    /// Every real-edge crossing is a sample on this azimuth.
    At(f64),
    /// Samples were requested at crossings (the face declines this round).
    Requested,
    /// No azimuth crosses the winding rings alone: the slit decides.
    None,
}

/// `b - a` reduced to `(-p/2, p/2]`.
fn reduced(p: f64, a: f64, b: f64) -> f64 {
    let d = along_offset(Some(p), a, b, 1.0);
    if d > 0.5 * p { d - p } else { d }
}

/// A ring's net winding of the along period, from its consecutive steps.
fn winding(cut: &Cut, p: f64, ring: &BoundaryRing) -> i64 {
    let n = ring.uv.len();
    let total: f64 = (0..n)
        .map(|i| reduced(p, cut.along(ring.uv[i]), cut.along(ring.uv[(i + 1) % n])))
        .sum();
    (total / p).round() as i64
}

/// The ring segments crossing the meridian at `az`: `(segment, sample on
/// it)` with the sample index when a sample lies exactly on it.
fn crossings(cut: &Cut, p: f64, ring: &BoundaryRing, az: f64) -> Vec<usize> {
    let n = ring.uv.len();
    let mut out = Vec::new();
    for i in 0..n {
        let a = reduced(p, az, cut.along(ring.uv[i]));
        let b = a + reduced(p, cut.along(ring.uv[i]), cut.along(ring.uv[(i + 1) % n]));
        // Half-open at the segment's end, so a sample on the meridian is
        // counted once, by the segment it starts.
        if (a <= 0.0 && b > 0.0) || (a >= 0.0 && b < 0.0) {
            out.push(i);
        }
    }
    out
}

/// How far the nearest turning point of the along coordinate on any ring
/// lies from `az`, read from positions.
fn turning_clearance(
    cut: &Cut,
    p: f64,
    s: &SurfaceGeometry,
    rings: &[BoundaryRing],
    az: f64,
) -> f64 {
    let mut clearance = f64::INFINITY;
    for ring in rings {
        let along: Vec<f64> = ring
            .uv
            .iter()
            .zip(&ring.pos)
            .map(|(&uv, q)| cut.along(s.closest_point(q).unwrap_or(uv)))
            .collect();
        let n = along.len();
        for i in 0..n {
            let back = reduced(p, along[(i + n - 1) % n], along[i]);
            let ahead = reduced(p, along[i], along[(i + 1) % n]);
            if back * ahead <= 0.0 {
                clearance = clearance.min(reduced(p, az, along[i]).abs());
            }
        }
    }
    clearance
}

/// Choose the band's meridian and make its real-edge crossings samples.
pub(super) fn band_meridian(
    brep: &BRep,
    (surface, cut): (&SurfaceGeometry, &Cut),
    rings: &mut [BoundaryRing],
    split_requests: &mut Vec<SplitRequest>,
) -> Meridian {
    let Some(p) = cut.period else {
        return Meridian::Unneeded;
    };
    let winds: Vec<bool> = rings.iter().map(|r| winding(cut, p, r) != 0).collect();
    if !winds.iter().any(|&w| w) {
        return Meridian::Unneeded;
    }
    // (azimuth, crossings, clearance) of each valid candidate.
    let valid: Vec<(f64, usize, f64)> = (0..CANDIDATES)
        .filter_map(|k| {
            let az = cut.start + p * (k as f64 + 0.5) / CANDIDATES as f64;
            let mut total = 0usize;
            for (ring, &w) in rings.iter().zip(&winds) {
                let c = crossings(cut, p, ring, az).len();
                if !w && c > 0 {
                    return None;
                }
                total += c;
            }
            (total >= 2 && total.is_multiple_of(2))
                .then(|| (az, total, turning_clearance(cut, p, surface, rings, az)))
        })
        .collect();
    let margin = p / (2.0 * CANDIDATES as f64);
    let clear = valid
        .iter()
        .filter(|c| c.2 > margin)
        .min_by_key(|c| c.1)
        .or_else(|| valid.iter().max_by(|a, b| a.2.total_cmp(&b.2)));
    let Some(&(az, _, _)) = clear else {
        return Meridian::None;
    };
    let mut requested = false;
    for ring in rings.iter_mut() {
        let n = ring.uv.len();
        for i in crossings(cut, p, ring, az) {
            let j = (i + 1) % n;
            let (edge, fa, fb) = ring.sources[i];
            // A row's samples lie on an iso-line: its wall point is exact.
            if edge == EdgeHandle::default() {
                continue;
            }
            let on = |uv: (f64, f64)| reduced(p, az, cut.along(uv)) == 0.0;
            if on(ring.uv[i]) {
                continue;
            }
            let Some(f) = crossing_fraction(brep, (surface, cut, p, az), edge, (fa, fb)) else {
                continue;
            };
            let dedup = crate::encroachment::SAMPLE_DEDUP;
            let snap = |uv: &mut (f64, f64)| {
                let along = cut.along(*uv) - reduced(p, az, cut.along(*uv));
                match cut.cap.singular.axis {
                    popkern_core::geometry::chart::Axis::V => uv.0 = along,
                    popkern_core::geometry::chart::Axis::U => uv.1 = along,
                }
            };
            if (f - fa).abs() < dedup {
                snap(&mut ring.uv[i]);
            } else if (f - fb).abs() < dedup {
                snap(&mut ring.uv[j]);
            } else {
                split_requests.push((edge, f));
                requested = true;
            }
        }
    }
    if requested {
        Meridian::Requested
    } else {
        Meridian::At(az)
    }
}

/// The curve fraction in `(fa, fb)` where `edge` crosses the meridian.
fn crossing_fraction(
    brep: &BRep,
    (surface, cut, p, az): (&SurfaceGeometry, &Cut, f64, f64),
    edge: EdgeHandle,
    (fa, fb): (f64, f64),
) -> Option<f64> {
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let gap = |f: f64| {
        let q = curve
            .geometry
            .point_at(curve.t_start + (curve.t_end - curve.t_start) * f);
        surface
            .closest_point(&q)
            .map(|uv| reduced(p, az, cut.along(uv)))
    };
    let (mut lo, mut hi) = (fa, fb);
    let (mut glo, ghi) = (gap(lo)?, gap(hi)?);
    if glo == 0.0 {
        return Some(lo);
    }
    if ghi == 0.0 {
        return Some(hi);
    }
    if glo.signum() == ghi.signum() {
        return None;
    }
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let g = gap(mid)?;
        if g.signum() == glo.signum() {
            lo = mid;
            glo = g;
        } else {
            hi = mid;
        }
        if (hi - lo).abs() <= 1e-15 * (fa - fb).abs().max(1e-300) {
            break;
        }
    }
    Some(0.5 * (lo + hi))
}
