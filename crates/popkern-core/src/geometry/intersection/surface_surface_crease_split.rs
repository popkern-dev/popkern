//! Splitting marched sample runs at C0 crease crossings (PK-442).
//!
//! A record that runs across a faceted chart's C0 knot line has a
//! genuine tangent corner there; a single smooth NURBS fit rounds it
//! (measured 1.1e-2 on a unit prism crease — four orders of magnitude
//! over tolerance). Cutting the run at each crossing, with the exact
//! crossing refined by the pinned Newton (pin u = the knot value),
//! gives every piece a smooth interior that fits to tolerance, and
//! the pieces chain through the shared refined junction downstream —
//! the same contract as seam splitting.

use crate::geometry::intersection::marching::marching_newton::newton_boundary_crossing;
use crate::geometry::intersection::marching::{BoundaryPin, MarchSample};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// Split each sample run wherever a chart coordinate crosses a C0
/// interior knot of its surface. Runs that never cross a crease come
/// back unchanged.
pub(crate) fn split_samples_at_c0_creases(
    segments: Vec<Vec<MarchSample>>,
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Vec<Vec<MarchSample>> {
    let (c1u, c1v) = crate::geometry::chart_periods::c0_interior_knots(surface1);
    let (c2u, c2v) = crate::geometry::chart_periods::c0_interior_knots(surface2);
    if c1u.is_empty() && c1v.is_empty() && c2u.is_empty() && c2v.is_empty() {
        return segments;
    }
    let knots = [&c1u, &c1v, &c2u, &c2v];
    let coords = |s: &MarchSample| [s.u1, s.v1, s.u2, s.v2];
    // First crease strictly between two consecutive samples, as
    // (fraction, axis, knot value).
    let first_crossing = |a: &MarchSample, b: &MarchSample| -> Option<(f64, usize, f64)> {
        let ca = coords(a);
        let cb = coords(b);
        let mut best: Option<(f64, usize, f64)> = None;
        for axis in 0..4 {
            for &k in knots[axis] {
                let (lo, hi) = (ca[axis].min(cb[axis]), ca[axis].max(cb[axis]));
                if k <= lo + 1.0e-12 || k >= hi - 1.0e-12 {
                    continue;
                }
                let t = (k - ca[axis]) / (cb[axis] - ca[axis]);
                if best.is_none_or(|(bt, _, _)| t < bt) {
                    best = Some((t, axis, k));
                }
            }
        }
        best
    };
    let mut out: Vec<Vec<MarchSample>> = Vec::with_capacity(segments.len());
    for seg in segments {
        if seg.len() < 2 {
            out.push(seg);
            continue;
        }
        let mut current: Vec<MarchSample> = vec![seg[0].clone()];
        for w in seg.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            if let Some((t, axis, k)) = first_crossing(a, b) {
                let lerp = |x: f64, y: f64| x + t * (y - x);
                let g = [
                    lerp(a.u1, b.u1),
                    lerp(a.v1, b.v1),
                    lerp(a.u2, b.u2),
                    lerp(a.v2, b.v2),
                ];
                let pin = BoundaryPin {
                    on_surface_one: axis < 2,
                    axis_u: axis % 2 == 0,
                    value: k,
                };
                // Sanity: the crossing must stay in the straddling
                // pair's neighbourhood; otherwise keep the run uncut
                // here (a fit across one unrefined crease beats a cut
                // at a wrong point).
                let refined = newton_boundary_crossing(
                    surface1, g[0], g[1], surface2, g[2], g[3], pin, tolerance,
                )
                .filter(|s| {
                    let gap = a.point.distance_to(&b.point);
                    s.point.distance_to(&a.point) <= gap * 4.0 + tolerance.linear()
                        && s.point.distance_to(&b.point) <= gap * 4.0 + tolerance.linear()
                });
                if let Some(s) = refined {
                    current.push(s.clone());
                    if current.len() >= 2 {
                        out.push(std::mem::take(&mut current));
                    } else {
                        current.clear();
                    }
                    current.push(s);
                }
            }
            current.push(b.clone());
        }
        if current.len() >= 2 {
            out.push(current);
        }
    }
    out
}
