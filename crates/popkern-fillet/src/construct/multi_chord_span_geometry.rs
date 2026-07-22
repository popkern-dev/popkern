//! Span-level geometry for multi-chord face splits: segment
//! intersection, sub-edge location, split bookkeeping, and the
//! trim-true sub-span extraction (PK-139 slice 6). Pure helpers for
//! `multi_chord`; nothing here mutates the BRep beyond bookkeeping.

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// Compute the intersection of two 3D line segments.
///
/// Returns `(s, t, point)` where `s ∈ (0,1)` is the parameter
/// along segment `a → b` and `t ∈ (0,1)` the parameter along
/// `c → d`. Returns `None` for parallel, skew, or endpoint-only
/// intersections.
pub(super) fn intersect_segments(
    a: Point3,
    b: Point3,
    c: Point3,
    d: Point3,
) -> Option<(f64, f64, Point3)> {
    let d1 = b - a;
    let d2 = d - c;
    let cross = d1.cross(&d2);
    let cross_sq = cross.dot(&cross);
    if cross_sq < 1e-18 {
        return None;
    }
    let diff = c - a;
    let s = diff.cross(&d2).dot(&cross) / cross_sq;
    let t = diff.cross(&d1).dot(&cross) / cross_sq;
    let interior = |x: f64| x > 1e-4 && x < 1.0 - 1e-4;
    if !interior(s) || !interior(t) {
        return None;
    }
    let pa = Point3::new(a.x + s * d1.x, a.y + s * d1.y, a.z + s * d1.z);
    let pb = Point3::new(c.x + t * d2.x, c.y + t * d2.y, c.z + t * d2.z);
    let gap_sq = (pa.x - pb.x).powi(2) + (pa.y - pb.y).powi(2) + (pa.z - pb.z).powi(2);
    if gap_sq > 1e-10 {
        return None;
    }
    Some((s, t, pa))
}

/// Find the sub-edge in `sub_edges` whose line segment contains
/// `point` (within tolerance), returning the edge and its
/// parameter along the edge in `(0.001, 0.999)`.
pub(super) fn find_sub_edge_containing_point(
    brep: &BRep,
    sub_edges: &[EdgeHandle],
    point: Point3,
) -> Option<(EdgeHandle, f64)> {
    for &eh in sub_edges {
        let edata = brep.edges.get(eh)?;
        let he = brep.half_edges.get(edata.half_edge)?;
        let twin = brep.half_edges.get(he.twin)?;
        let v0 = brep.vertices[he.origin].position;
        let v1 = brep.vertices[twin.origin].position;
        let dx = v1.x - v0.x;
        let dy = v1.y - v0.y;
        let dz = v1.z - v0.z;
        let len_sq = dx * dx + dy * dy + dz * dz;
        if len_sq < 1e-18 {
            continue;
        }
        let t = ((point.x - v0.x) * dx + (point.y - v0.y) * dy + (point.z - v0.z) * dz) / len_sq;
        if t <= 1e-4 || t >= 1.0 - 1e-4 {
            continue;
        }
        let px = v0.x + t * dx;
        let py = v0.y + t * dy;
        let pz = v0.z + t * dz;
        let dist_sq = (point.x - px).powi(2) + (point.y - py).powi(2) + (point.z - pz).powi(2);
        if dist_sq < 1e-10 {
            return Some((eh, t));
        }
    }
    None
}

/// Replace `old` with `[new_a, new_b]` in every `Vec<EdgeHandle>`
/// in `per_trim_sub_edges`. Used after `split_edge` to keep the
/// per-trim sub-edge bookkeeping in sync with the BRep.
pub(super) fn remap_split(
    per_trim_sub_edges: &mut [Vec<EdgeHandle>],
    old: EdgeHandle,
    new_a: EdgeHandle,
    new_b: EdgeHandle,
) {
    for list in per_trim_sub_edges.iter_mut() {
        if let Some(pos) = list.iter().position(|&e| e == old) {
            list.splice(pos..=pos, [new_a, new_b]);
        }
    }
}

/// Sub-span of the fillet's actual 3-D trim curve between two anchor
/// positions, as a uniform degree-1 polyline. `None` when either
/// anchor does not lie on the trim (corner sub-chords) or the span is
/// effectively straight (the chord is already exact).
pub(super) fn trim_subcurve_between(
    trim: &NurbsCurve,
    p_a: &Point3,
    p_b: &Point3,
    tolerance: &Tolerance,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::CurveEval;
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
    let (t0, t1) = trim.parameter_range();
    let on_trim = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(trim.clone()),
        t_start: t0,
        t_end: t1,
    };
    let gate = tolerance.linear().max(1.0e-9) * 10.0;
    let ta =
        popkern_core::geometry::curve_inversion::invert_curve_at_point(&on_trim, *p_a, tolerance)
            .ok()?;
    let tb =
        popkern_core::geometry::curve_inversion::invert_curve_at_point(&on_trim, *p_b, tolerance)
            .ok()?;
    if trim.point_at(ta).distance_to(p_a) > gate || trim.point_at(tb).distance_to(p_b) > gate {
        return None;
    }
    let (lo, hi) = if ta <= tb { (ta, tb) } else { (tb, ta) };
    if (hi - lo).abs() < 1.0e-12 {
        return None;
    }
    let n = 17;
    let mut pts = Vec::with_capacity(n);
    for k in 0..n {
        let t = lo + (hi - lo) * k as f64 / (n - 1) as f64;
        pts.push(trim.point_at(t));
    }
    if ta > tb {
        pts.reverse();
    }
    // Endpoints must be EXACTLY the anchors (split_face matches them
    // to vertices).
    pts[0] = *p_a;
    let last = pts.len() - 1;
    pts[last] = *p_b;
    // Straight span: let the chord path (and its analytic
    // canonicalization) serve it.
    let chord = *p_b - *p_a;
    let chord_len = chord.length();
    if chord_len > 1.0e-12 {
        let dir = chord * (1.0 / chord_len);
        let max_dev = pts
            .iter()
            .map(|p| {
                let d = *p - *p_a;
                (d - dir * d.dot(&dir)).length()
            })
            .fold(0.0_f64, f64::max);
        if max_dev <= gate {
            return None;
        }
    }
    // EXACT SLICE of the trim NURBS — the former 17-point degree-1
    // polyline stand-in carries chord sagitta r·(1−cos(span/32)),
    // which is 0.42 on a 5.87-rad interrupted-rim arc (measured: the
    // z=7 contact edge's boundary-first inversion rejected at
    // 3.8e-1 and both adjacent faces fell to the cracking legacy
    // mesher). Short sub-chords keep sub-tolerance error either way.
    let sliced = {
        let eps = (trim.parameter_range().1 - trim.parameter_range().0) * 1e-9;
        let (t_lo_full, _) = trim.parameter_range();
        let right = if lo > t_lo_full + eps {
            trim.split_at(lo).ok()?.1
        } else {
            trim.clone()
        };
        let (_, r_hi) = right.parameter_range();
        if hi < r_hi - eps {
            right.split_at(hi).ok()?.0
        } else {
            right
        }
    };
    let mut out = sliced;
    if ta > tb {
        // Reverse: mirror control points/weights and the knot vector.
        out.control_points.reverse();
        out.weights.reverse();
        let (k0, k1) = (out.knots[0], *out.knots.last()?);
        out.knots = out.knots.iter().rev().map(|&k| k0 + k1 - k).collect();
    }
    // The slice's endpoints are EXACT trim geometry; the caller
    // snaps the anchor VERTICES onto them (pinning control points to
    // walk-noisy vertices bent the end spans off-surface by the
    // anchor error — measured 1.6e-2 on the z=8 contact edge).
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersect_crossing_segments() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 1.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);
        let d = Point3::new(1.0, 0.0, 0.0);
        let result = intersect_segments(a, b, c, d).unwrap();
        assert!((result.0 - 0.5).abs() < 1e-9);
        assert!((result.1 - 0.5).abs() < 1e-9);
        assert!((result.2.x - 0.5).abs() < 1e-9);
        assert!((result.2.y - 0.5).abs() < 1e-9);
    }

    #[test]
    fn parallel_segments_no_intersection() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);
        let d = Point3::new(1.0, 1.0, 0.0);
        assert!(intersect_segments(a, b, c, d).is_none());
    }

    #[test]
    fn endpoint_touching_segments_no_intersection() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(1.0, 0.0, 0.0);
        let d = Point3::new(1.0, 1.0, 0.0);
        assert!(intersect_segments(a, b, c, d).is_none());
    }

    #[test]
    fn skew_segments_no_intersection() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(0.5, -1.0, 1.0);
        let d = Point3::new(0.5, 1.0, 1.0);
        assert!(intersect_segments(a, b, c, d).is_none());
    }
}
