//! The exact gap between two planar faces bounded by segments, for the
//! self-contact search (PK-1155).
//!
//! Two faces of one valid shell do not cross, so their nearest points
//! lie on a boundary of one, or on both interiors with the segment
//! between normal to both: the planes parallel, and then a boundary
//! point attains the same gap (slide the interior point until it or its
//! foot meets a boundary). So the gap between two such faces is the
//! least gap from an edge of either to the other face's region, and a
//! segment's gap to a planar polygon is closed form: nought where it
//! crosses the plane inside, else the least of its ends' heights where
//! their feet lie inside and its distances to the boundary's segments
//! (where its foot leaves the region it is no nearer the region than
//! the boundary). A faceted body's triangles, each pair cleared by the
//! generic edge test's projections and region queries in half a
//! millisecond, are cleared in microseconds.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use super::clear::{plane_of, straight};

/// A planar face whose every boundary edge is a segment: its plane's
/// point, unit normal and two directions in it, and each boundary
/// segment with its ends' coordinates in the plane.
pub(super) struct Polygon {
    origin: Point3,
    normal: Vec3,
    axes: (Vec3, Vec3),
    segments: Vec<(Point3, Point3)>,
    flat: Vec<((f64, f64), (f64, f64))>,
    /// The box of its segments, which holds the region.
    bx: (Point3, Point3),
}

/// The face `face` as a [`Polygon`], when it is planar and every edge of
/// its boundary loops is a segment; its vertex loops lie inside it and
/// add nothing.
pub(super) fn polygon_of(brep: &BRep, face: FaceHandle) -> Option<Polygon> {
    let f = brep.faces.get(face)?;
    let (origin, normal) = plane_of(&f.surface.as_ref()?.geometry)?;
    let e1 = normal.cross(&Vec3::new(1.0, 0.0, 0.0));
    let e1 = if e1.length() < 0.5 {
        normal.cross(&Vec3::new(0.0, 1.0, 0.0))
    } else {
        e1
    }
    .normalized()?;
    let e2 = normal.cross(&e1);
    let mut segments = Vec::new();
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let start = brep.loops.get(l)?.half_edge;
        if brep.half_edges.get(start).is_none() {
            continue;
        }
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let h = brep.half_edges.get(he)?;
            let c = brep.edges.get(h.edge)?.curve.as_ref()?;
            if !straight(c) {
                return None;
            }
            segments.push((c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end)));
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    let at = |p: &Point3| ((*p - origin).dot(&e1), (*p - origin).dot(&e2));
    let flat = segments.iter().map(|(a, b)| (at(a), at(b))).collect();
    let bx = box_of(&segments)?;
    Some(Polygon {
        origin,
        normal,
        axes: (e1, e2),
        segments,
        flat,
        bx,
    })
}

impl Polygon {
    /// Whether the point of the plane at `(x, y)` lies inside the region
    /// (crossings of a ray along the first axis, every loop together).
    fn holds(&self, (x, y): (f64, f64)) -> bool {
        let mut inside = false;
        for &((x0, y0), (x1, y1)) in &self.flat {
            if (y0 > y) != (y1 > y) && x < x0 + (y - y0) / (y1 - y0) * (x1 - x0) {
                inside = !inside;
            }
        }
        inside
    }

    fn flat_of(&self, p: &Point3) -> (f64, f64) {
        let d = *p - self.origin;
        (d.dot(&self.axes.0), d.dot(&self.axes.1))
    }

    /// The least distance from the segment `a`-`b` to the region.
    #[cfg(test)]
    pub(super) fn segment_gap(&self, a: Point3, b: Point3) -> f64 {
        self.segment_gap_within(a, b, f64::INFINITY)
    }

    /// The least distance from the segment `a`-`b` to the region when it
    /// is at most `cutoff`, else a distance over `cutoff`: boundary
    /// segments whose box lies farther than `cutoff` from the segment's
    /// are passed over.
    pub(super) fn segment_gap_within(&self, a: Point3, b: Point3, cutoff: f64) -> f64 {
        let (ha, hb) = (
            (a - self.origin).dot(&self.normal),
            (b - self.origin).dot(&self.normal),
        );
        if (ha <= 0.0) != (hb <= 0.0) || ha == 0.0 || hb == 0.0 {
            let s = if ha == hb { 0.0 } else { ha / (ha - hb) };
            if self.holds(self.flat_of(&(a + (b - a) * s))) {
                return 0.0;
            }
        }
        let mut best = f64::INFINITY;
        for (p, h) in [(a, ha), (b, hb)] {
            if h.abs() < best && self.holds(self.flat_of(&p)) {
                best = h.abs();
            }
        }
        let near = (lower(a, b), upper(a, b));
        for &(c, d) in &self.segments {
            if box_gap(near, (lower(c, d), upper(c, d))) <= cutoff.min(best) {
                best = best.min(segment_distance(a, b, c, d));
            }
        }
        best
    }
}

/// Whether the regions of `p` and `q` lie farther than `tol` apart:
/// every segment of either farther than that from the other's region,
/// at once when its box is farther than that from the region's box.
pub(super) fn polygons_apart(p: &Polygon, q: &Polygon, tol: f64) -> bool {
    let clear = |from: &Polygon, to: &Polygon| {
        from.segments.iter().all(|&(a, b)| {
            box_gap((lower(a, b), upper(a, b)), to.bx) > tol
                || to.segment_gap_within(a, b, tol) > tol
        })
    };
    clear(p, q) && clear(q, p)
}

/// The box of `segments`, when there are any.
pub(super) fn box_of(segments: &[(Point3, Point3)]) -> Option<(Point3, Point3)> {
    let mut it = segments.iter().flat_map(|(a, b)| [*a, *b]);
    let first = it.next()?;
    Some(it.fold((first, first), |(lo, hi), p| (lower(lo, p), upper(hi, p))))
}

fn lower(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z))
}

fn upper(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z))
}

/// The distance between two boxes, nought when they meet.
fn box_gap((alo, ahi): (Point3, Point3), (blo, bhi): (Point3, Point3)) -> f64 {
    let gap = |a0: f64, a1: f64, b0: f64, b1: f64| (b0 - a1).max(a0 - b1).max(0.0);
    let (x, y, z) = (
        gap(alo.x, ahi.x, blo.x, bhi.x),
        gap(alo.y, ahi.y, blo.y, bhi.y),
        gap(alo.z, ahi.z, blo.z, bhi.z),
    );
    (x * x + y * y + z * z).sqrt()
}

/// The least distance between the segments `p0`-`p1` and `q0`-`q1`
/// (each parameter clamped in turn to its segment, the other refitted).
pub(super) fn segment_distance(p0: Point3, p1: Point3, q0: Point3, q1: Point3) -> f64 {
    let (d1, d2, r) = (p1 - p0, q1 - q0, p0 - q0);
    let (a, e, f) = (d1.dot(&d1), d2.dot(&d2), d2.dot(&r));
    let tiny = 1e-300;
    let (s, t) = if a <= tiny && e <= tiny {
        (0.0, 0.0)
    } else if a <= tiny {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = d1.dot(&r);
        if e <= tiny {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = d1.dot(&d2);
            let denom = a * e - b * b;
            let s = if denom > tiny {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t = (b * s + f) / e;
            if t < 0.0 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else if t > 1.0 {
                (((b - c) / a).clamp(0.0, 1.0), 1.0)
            } else {
                (s, t)
            }
        }
    };
    (p0 + d1 * s).distance_to(&(q0 + d2 * t))
}

#[cfg(test)]
#[path = "self_contact_polygon_tests.rs"]
mod tests;
