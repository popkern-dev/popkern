//! A cone face holding its apex, cut along one ruling to the apex for
//! its development (split from `boundary_first_cone_development.rs`
//! under the file-size rule).
//!
//! A face whose outer ring turns once round the axis without passing the
//! apex holds the apex inside it (a plane missing the apex by a few
//! tolerances cuts an ellipse whose near end sweeps round it; the face
//! is the cone less the thin cap that ellipse bounds). Unrolled, such a
//! face spans `2π sin α` round the origin and is not one sheet until cut.
//! The cut is the ruling from one ring sample down to the apex, met from
//! both sides: the ring is closed one turn on to that sample's copy, down
//! the ruling to the apex and back up it to the sample, the two walls
//! carrying the same stations, so the face is the ring closing round the
//! axis through the apex, which the development serves (two rays from
//! the origin, as a seam's two uses). The face chart's apex cap does the
//! same along a degenerate apex row, where a sweep a few tolerances from
//! the apex runs nearly along that row and the meshing fails.
//!
//! The cut sample is the one whose ruling crosses no ring between it and
//! the apex, farthest from every other sample nearer the apex than it
//! (the walls are interior to the face; a wall grazing the ring would
//! leave slivers). The wall stations are uniform along the ruling, as
//! many as the ring's mean step needs: a ruling is straight, so its
//! stations only set the lattice density beside it. Their spans are
//! zero-length at the cut sample (the segment refinement cannot split
//! them, as the face chart's apex-cap walls).

use std::f64::consts::TAU;

use popkern_math::PortableMath;

use popkern_core::geometry::Point3;

use super::super::BoundaryRing;

/// The most stations a wall carries.
const MAX_WALL_STATIONS: usize = 512;

/// The most samples tried as the cut.
const MAX_CANDIDATES: usize = 256;

/// The apex, `sin α` and a position's axial height above the apex.
type Axis<'h> = (Point3, f64, &'h dyn Fn(&Point3) -> f64);

/// A ring segment: its outer-ring index (`usize::MAX` on a hole), its
/// end azimuths and slants.
type Segment = (usize, (f64, f64), (f64, f64));

/// The face's outer ring cut along a ruling to the apex and closed back
/// through it, or `None` when every ruling from a sample to the apex
/// crosses a ring. `az` are the outer ring's unwrapped azimuths (turning
/// `turns`, ±1, from its first sample back to it), `holes` the holes'
/// rings and azimuths, `slant(p)` a position's slant distance from the
/// apex and `height(p)` its axial height.
pub(super) fn cut_to_apex(
    ring: &BoundaryRing,
    az: &[f64],
    turns: f64,
    holes: &[(&BoundaryRing, &[f64])],
    (apex, sin_a, height): Axis<'_>,
) -> Option<BoundaryRing> {
    let n = ring.pos.len();
    if n < 3 || az.len() != n {
        return None;
    }
    let slant = |p: &Point3| p.distance_to(&apex);
    // Every ring segment with its end azimuths and slants, the outer
    // ring's closing one a turn on.
    let mut segments: Vec<Segment> = Vec::new();
    for j in 0..n {
        let k = (j + 1) % n;
        let ak = if k == 0 { az[0] + turns * TAU } else { az[k] };
        segments.push((j, (az[j], ak), (slant(&ring.pos[j]), slant(&ring.pos[k]))));
    }
    for (h, haz) in holes {
        let m = h.pos.len();
        for j in 0..m {
            let k = (j + 1) % m;
            segments.push((
                usize::MAX,
                (haz[j], haz[k]),
                (slant(&h.pos[j]), slant(&h.pos[k])),
            ));
        }
    }
    let samples: Vec<(f64, f64)> = (0..n)
        .map(|j| (az[j], slant(&ring.pos[j])))
        .chain(
            holes
                .iter()
                .flat_map(|(h, haz)| haz.iter().zip(&h.pos).map(|(&a, p)| (a, slant(p)))),
        )
        .collect();
    let stride = n.div_ceil(MAX_CANDIDATES).max(1);
    let mut best: Option<(usize, f64)> = None;
    for c in (0..n).step_by(stride) {
        let (ac, sc) = (az[c], slant(&ring.pos[c]));
        if sc.is_nan() || sc <= 0.0 {
            continue;
        }
        let crosses = segments.iter().any(|&(j, (a0, a1), (s0, s1))| {
            if j == c || j == (c + n - 1) % n {
                return false;
            }
            let (lo, hi) = (a0.min(a1), a0.max(a1));
            let k = ((lo - ac) / TAU).ceil();
            let at = ac + k * TAU;
            if at > hi || hi - lo <= 0.0 {
                return false;
            }
            let f = (at - a0) / (a1 - a0);
            s0 + f * (s1 - s0) < sc
        });
        if crosses {
            continue;
        }
        // The nearest other sample nearer the apex, measured in the
        // development from the wall (either copy).
        let clearance = samples
            .iter()
            .enumerate()
            .filter(|&(j, _)| j >= n || (j != c && j != (c + 1) % n && j != (c + n - 1) % n))
            .filter(|&(_, &(_, s))| s < sc)
            .map(|(_, &(a, s))| {
                let d = a - ac;
                let turn = (d - TAU * (d / TAU).round()).abs() * sin_a;
                if turn < 0.5 * std::f64::consts::PI {
                    s * turn.psin()
                } else {
                    s
                }
            })
            .fold(f64::INFINITY, f64::min);
        if best.is_none_or(|(_, b)| clearance > b) {
            best = Some((c, clearance));
        }
    }
    let (c, _) = best?;
    let p = ring.pos[c];
    let mean_step = (0..n)
        .map(|j| ring.pos[j].distance_to(&ring.pos[(j + 1) % n]))
        .sum::<f64>()
        / n as f64;
    let wall_len = slant(&p);
    if mean_step.is_nan() || wall_len.is_nan() || mean_step <= 0.0 || wall_len <= 0.0 {
        return None;
    }
    let m = ((wall_len / mean_step).ceil() as usize).clamp(1, MAX_WALL_STATIONS);
    let mut out = BoundaryRing {
        uv: Vec::with_capacity(n + 2 * m + 1),
        pos: Vec::with_capacity(n + 2 * m + 1),
        sources: Vec::with_capacity(n + 2 * m + 1),
    };
    for j in (c..n).chain(0..c) {
        out.uv.push(ring.uv[j]);
        out.pos.push(ring.pos[j]);
        out.sources.push(ring.sources[j]);
    }
    // The walls' spans: zero-length at the cut sample, on the edge whose
    // segment enters it.
    let entering = ring.sources[(c + n - 1) % n];
    let wall = (entering.0, entering.2, entering.2);
    let (u_cut, u_copy) = (ring.uv[c].0, ring.uv[c].0 + turns * TAU);
    let h_cut = height(&p);
    let station = |k: usize| {
        let f = k as f64 / m as f64;
        (apex + (p - apex) * f, h_cut * f)
    };
    // The cut sample's copy a turn on, down the ruling, the apex, back up.
    out.uv.push((u_copy, ring.uv[c].1));
    out.pos.push(p);
    out.sources.push(wall);
    for k in (1..m).rev() {
        let (q, h) = station(k);
        out.uv.push((u_copy, h));
        out.pos.push(q);
        out.sources.push(wall);
    }
    out.uv.push((u_cut, 0.0));
    out.pos.push(apex);
    out.sources.push(wall);
    for k in 1..m {
        let (q, h) = station(k);
        out.uv.push((u_cut, h));
        out.pos.push(q);
        out.sources.push(wall);
    }
    Some(out)
}
