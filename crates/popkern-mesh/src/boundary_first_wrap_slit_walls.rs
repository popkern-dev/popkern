//! Where a consecutively-unwrapped ring crosses the slit's walls
//! (`u* + k·p`): between samples, interpolated, or at a sample lying on a
//! wall, which an atlas band's meridian puts there (its real-edge
//! crossings are canonical samples, `boundary_first_atlas_meridian.rs`).
//! A wall sample crosses when its neighbours lie on opposite sides of the
//! wall, and touches (no crossing) when they lie on one side.

use popkern_math::PortableMath;

use super::BoundaryRing;

/// A wall crossing: the ring segment it lies on, its fraction along the
/// segment, its second coordinate, and whether it is the segment's first
/// sample itself (fraction 0, a wall sample).
#[derive(Debug, Clone, Copy)]
pub(super) struct WallCrossing {
    pub(super) segment: usize,
    pub(super) t: f64,
    pub(super) v: f64,
    pub(super) at_sample: bool,
}

/// Relative distance in periods within which a sample is on a wall.
const ON_WALL: f64 = 1e-12;

/// The wall index a first coordinate lies on, if any.
fn wall_of(u: f64, u_star: f64, period: f64) -> Option<i64> {
    let w = (u - u_star) / period;
    let k = w.round();
    ((w - k).abs() < ON_WALL).then_some(k as i64)
}

/// Every wall crossing of `ring`, in ring order; `None` when a segment
/// runs along a wall (the slit cannot cut there).
pub(super) fn wall_crossings(
    ring: &BoundaryRing,
    u_star: f64,
    period: f64,
) -> Option<Vec<WallCrossing>> {
    let n = ring.uv.len();
    if n == 0 {
        return Some(Vec::new());
    }
    // The ring closes by its winding: the sample after the last is the
    // first shifted by it.
    let closing = ((ring.uv[n - 1].0 - ring.uv[0].0 + 0.5 * period) / period).floor() * period;
    let after = |i: usize| -> f64 {
        if i + 1 == n {
            ring.uv[0].0 + closing
        } else {
            ring.uv[i + 1].0
        }
    };
    let before = |i: usize| -> f64 {
        if i == 0 {
            ring.uv[n - 1].0 - closing
        } else {
            ring.uv[i - 1].0
        }
    };
    let mut out = Vec::new();
    for i in 0..n {
        let (ua, va) = ring.uv[i];
        let ub = after(i);
        let vb = ring.uv[(i + 1) % n].1;
        let (ka, kb) = (wall_of(ua, u_star, period), wall_of(ub, u_star, period));
        if let Some(k) = ka {
            if kb == Some(k) {
                return None;
            }
            let wall = u_star + k as f64 * period;
            let (prev, next) = (before(i) - wall, ub - wall);
            if wall_of(before(i), u_star, period) == Some(k) {
                return None;
            }
            if prev * next < 0.0 {
                out.push(WallCrossing {
                    segment: i,
                    t: 0.0,
                    v: va,
                    at_sample: true,
                });
            }
        }
        let k_lo = ((ua.min(ub) - u_star) / period).ceil() as i64;
        let k_hi = ((ua.max(ub) - u_star) / period).floor() as i64;
        for k in k_lo..=k_hi {
            if Some(k) == ka || Some(k) == kb {
                continue;
            }
            let wall = u_star + k as f64 * period;
            let t = if (ub - ua).abs() < 1e-30 {
                0.0
            } else {
                ((wall - ua) / (ub - ua)).clamp(0.0, 1.0)
            };
            out.push(WallCrossing {
                segment: i,
                t,
                v: va + (vb - va) * t,
                at_sample: false,
            });
        }
    }
    Some(out)
}

/// Net u-winding of a consecutively-unwrapped ring, in periods
/// (rounded): +1/-1 for a full wrap, 0 for a chartable ring.
pub(in super::super) fn net_winding(ring: &BoundaryRing, period: f64) -> i64 {
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
pub(in super::super) fn collapse_stubs(ring: &mut BoundaryRing) {
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
                if (pu - qu).phypot(pv - qv) < 1.0e-9 {
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
