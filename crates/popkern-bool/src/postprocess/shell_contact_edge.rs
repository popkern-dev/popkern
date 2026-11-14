//! An edge of one shell touching a face of another (R5.8 of the
//! tangency foundations' review: an edge can graze a face).
//!
//! The face-pair search refines seeds by alternating projection
//! between the two surfaces, which converges to their unconstrained
//! closest approach. Where the touch is on a face's boundary the
//! surfaces themselves pass through each other beyond it, and that
//! approach lies off the trim: a ball resting on a box's edge was
//! never found, and the union served box and ball as two shells
//! touching at a point (validate_all failing on them). Here every edge
//! of the face is walked against the other face's surface: the
//! smallest distance from the edge to the surface, found from stations
//! along it and refined by ternary search, is a contact when it is
//! within the tolerance and its foot lies on the other face's trim.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

/// Stations per edge for the coarse search.
const STATIONS: usize = 32;
/// Stations refined per edge (the closest by distance).
const REFINED: usize = 3;
/// Ternary-search rounds per refined station.
const ROUNDS: usize = 60;

/// A point's distance to a surface, its foot there and the foot's uv.
type Gap = Option<(f64, Point3, (f64, f64))>;

/// Where an edge of `f` touches face `g` (surface `sg`), if anywhere;
/// an edge `clear` proves apart from `g` is not walked.
pub(super) fn edge_contact(
    brep: &BRep,
    f: FaceHandle,
    (g, sg): (FaceHandle, &SurfaceGeometry),
    linear: f64,
    clear: Option<&super::shell_contact::EdgeClear>,
) -> Option<Point3> {
    let face = brep.faces.get(f)?;
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());
    let gap = |p: &Point3| -> Gap {
        let (q, uv) = super::shell_contact::foot(sg, p)?;
        Some((p.distance_to(&q), q, uv))
    };
    for lh in loops {
        let start = brep.loops.get(lh)?.half_edge;
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let h = brep.half_edges.get(he)?;
            if let Some(c) = brep
                .edges
                .get(h.edge)
                .and_then(|e| e.curve.as_ref())
                .filter(|c| !clear.is_some_and(|clear| clear(f, c)))
            {
                let at = |s: f64| c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * s);
                let dist = |s: f64| gap(&at(s)).map_or(f64::INFINITY, |(d, ..)| d);
                let mut coarse: Vec<(f64, usize)> = (0..=STATIONS)
                    .map(|i| (dist(i as f64 / STATIONS as f64), i))
                    .collect();
                if let Some(p) = run_middle(&coarse, linear, &at, &gap, &|uv| {
                    super::shell_contact::on_trim(brep, g, sg, uv, linear)
                }) {
                    return Some(p);
                }
                coarse.sort_by(|a, b| a.0.total_cmp(&b.0));
                for &(_, i) in coarse.iter().take(REFINED) {
                    let step = 1.0 / STATIONS as f64;
                    let (mut lo, mut hi) = (
                        (i as f64 * step - step).max(0.0),
                        (i as f64 * step + step).min(1.0),
                    );
                    for _ in 0..ROUNDS {
                        let (m1, m2) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
                        if dist(m1) <= dist(m2) {
                            hi = m2;
                        } else {
                            lo = m1;
                        }
                    }
                    let p = at(0.5 * (lo + hi));
                    if let Some((d, q, uv)) = gap(&p)
                        && d <= linear
                        && super::shell_contact::on_trim(brep, g, sg, uv, linear)
                    {
                        return Some(p + (q - p) * 0.5);
                    }
                }
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    None
}

/// The middle of the first run of consecutive stations within `linear`
/// of the surface whose middle lies on the trim: an edge lying along
/// the other face touches it along a line, and its middle (not an end,
/// which may be a corner of both faces, where the extent probe reads
/// several directions on both and calls a line an area) names it.
fn run_middle(
    coarse: &[(f64, usize)],
    linear: f64,
    at: &dyn Fn(f64) -> Point3,
    gap: &dyn Fn(&Point3) -> Gap,
    on_trim: &dyn Fn((f64, f64)) -> bool,
) -> Option<Point3> {
    let mut i = 0;
    while i < coarse.len() {
        if coarse[i].0 > linear {
            i += 1;
            continue;
        }
        let start = i;
        while i + 1 < coarse.len() && coarse[i + 1].0 <= linear {
            i += 1;
        }
        if i > start {
            let s = 0.5 * (coarse[start].1 + coarse[i].1) as f64 / STATIONS as f64;
            let p = at(s);
            if let Some((d, q, uv)) = gap(&p)
                && d <= linear
                && on_trim(uv)
            {
                return Some(p + (q - p) * 0.5);
            }
        }
        i += 1;
    }
    None
}
