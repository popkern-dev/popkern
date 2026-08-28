//! Span-level geometry for multi-chord face splits: segment
//! intersection, sub-edge location, split bookkeeping, and the
//! trim-true sub-span extraction (PK-139 slice 6). Pure helpers for
//! `multi_chord`; nothing here mutates the BRep beyond bookkeeping.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};
use popkern_core::topology::{FaceHandle, VertexHandle};

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

/// Closest-approach crossing of two SAMPLED trim curves (PK-141:
/// straight anchor chords miss every crossing on a curved support —
/// a boss-wall chain's vertical rail crosses the horizontal rims'
/// trims on the r=1 cylinder with the chords passing 5e-3 inside
/// it). Samples both curves, finds the closest segment pair, and
/// accepts when the approach distance is within `gate`. Returns
/// `(frac_a, frac_b, point)` with fractions in sample space [0, 1].
pub(super) fn intersect_trim_curves(
    trim_a: &NurbsCurve,
    trim_b: &NurbsCurve,
    gate: f64,
) -> Option<(f64, f64, Point3)> {
    use popkern_core::geometry::curve::CurveEval;
    const N: usize = 33;
    let sample = |c: &NurbsCurve| -> Vec<Point3> {
        let (t0, t1) = c.parameter_range();
        (0..N)
            .map(|k| c.point_at(t0 + (t1 - t0) * k as f64 / (N - 1) as f64))
            .collect()
    };
    let pa = sample(trim_a);
    let pb = sample(trim_b);
    let mut best: Option<(f64, f64, Point3, f64)> = None;
    for i in 0..N - 1 {
        for j in 0..N - 1 {
            let (s, t, d, pt) = closest_segment_pair(pa[i], pa[i + 1], pb[j], pb[j + 1]);
            let fa = (i as f64 + s) / (N - 1) as f64;
            let fb = (j as f64 + t) / (N - 1) as f64;
            if best.as_ref().map(|(_, _, _, bd)| d < *bd).unwrap_or(true) {
                best = Some((fa, fb, pt, d));
            }
        }
    }
    let (fa, fb, pt, d) = best?;
    let interior = |x: f64| x > 1e-3 && x < 1.0 - 1e-3;
    if d <= gate && interior(fa) && interior(fb) {
        Some((fa, fb, pt))
    } else {
        None
    }
}

/// Closest points between two segments; returns `(s, t, distance,
/// midpoint of the closest pair)`.
fn closest_segment_pair(a0: Point3, a1: Point3, b0: Point3, b1: Point3) -> (f64, f64, f64, Point3) {
    let u = a1 - a0;
    let v = b1 - b0;
    let w = a0 - b0;
    let (uu, uv, vv, uw, vw) = (u.dot(&u), u.dot(&v), v.dot(&v), u.dot(&w), v.dot(&w));
    let denom = uu * vv - uv * uv;
    let (mut s, mut t) = if denom.abs() > 1e-18 {
        (
            ((uv * vw - vv * uw) / denom).clamp(0.0, 1.0),
            ((uu * vw - uv * uw) / denom).clamp(0.0, 1.0),
        )
    } else {
        (0.0, (vw / vv.max(1e-18)).clamp(0.0, 1.0))
    };
    // One clamp round to keep the pair consistent.
    if uu > 1e-18 {
        s = ((t * uv - uw) / uu).clamp(0.0, 1.0);
    }
    if vv > 1e-18 {
        t = ((s * uv + vw) / vv).clamp(0.0, 1.0);
    }
    let p = Point3::new(a0.x + s * u.x, a0.y + s * u.y, a0.z + s * u.z);
    let q = Point3::new(b0.x + t * v.x, b0.y + t * v.y, b0.z + t * v.z);
    let mid = Point3::new(0.5 * (p.x + q.x), 0.5 * (p.y + q.y), 0.5 * (p.z + q.z));
    (s, t, p.distance_to(&q), mid)
}

/// PK-141: move a pure mev-spur anchor (`tip`) into the sub-face
/// that contains its chord `partner`. The old spur is `kev`-removed
/// and a fresh spur is rooted at that face's nearest boundary
/// vertex; the tip position is preserved. Returns the replacement
/// vertex and the host face, or `None` when `tip` is not a spur or
/// no host contains the partner.
pub(super) fn rehome_spur_tip(
    brep: &mut BRep,
    tip: VertexHandle,
    partner: VertexHandle,
    sub_faces: &[FaceHandle],
) -> Option<(VertexHandle, FaceHandle)> {
    // The spur out-half-edge: origin = root, dest = tip, and the mev
    // shape (out.next == out.twin). Tip must have no other edges.
    let mut spur_out = None;
    let mut incident = 0usize;
    for (heh, he) in brep.half_edges.iter() {
        let Some(tw) = brep.half_edges.get(he.twin) else {
            continue;
        };
        if tw.origin == tip {
            incident += 1;
            if he.next == he.twin {
                spur_out = Some(heh);
            }
        }
        if he.origin == tip && he.twin != heh {
            // counted via the twin side; nothing extra
        }
    }
    if incident != 1 {
        return None;
    }
    let spur_out = spur_out?;
    let pos = brep.vertices.get(tip)?.position;
    let host = super::face_trim_pipeline::find_face_containing_vertices(
        brep, sub_faces, partner, partner,
    )?;
    // Nearest boundary origin on the host face.
    let fd = brep.faces.get(host)?;
    let mut loops = vec![fd.outer_loop];
    loops.extend(fd.inner_loops.iter().copied());
    let mut best: Option<(popkern_core::topology::HalfEdgeHandle, f64)> = None;
    for lh in loops {
        let l = brep.loops.get(lh)?;
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let he = brep.half_edges.get(cur)?;
            if let Some(v) = brep.vertices.get(he.origin) {
                let d = v.position.distance_to(&pos);
                if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some((cur, d));
                }
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }
    let (root_he, _) = best?;
    let root_pos = brep
        .half_edges
        .get(root_he)
        .and_then(|h| brep.vertices.get(h.origin))
        .map(|v| v.position)?;
    // Transactional: kev-then-mev on a candidate clone, committed
    // only when both succeed — a failed mev after an in-place kev
    // destroyed the spur anchor and left the caller with a dangling
    // vertex handle (the same disease as the anchor transplant;
    // measured as a SlotMap panic on the skewed-rib junction trim).
    let mut candidate = brep.clone();
    candidate.kev(spur_out).ok()?;
    let mev = candidate.mev(root_he, pos).ok()?;
    if let Some(e) = candidate.edges.get_mut(mev.edge) {
        e.curve = popkern_core::geometry::CurveOnEdge::line_segment(root_pos, pos, 1e-9);
    }
    *brep = candidate;
    if std::env::var_os("PK_MULTICHORD_TRACE").is_some() {
        eprintln!(
            "[mc] re-homed spur tip {tip:?} -> {:?} into {host:?} at {pos:?}",
            mev.vertex
        );
    }
    Some((mev.vertex, host))
}
/// Trace helper: the faces each descendant edge currently borders.
pub(super) fn descendant_holders(
    brep: &BRep,
    descendants: &[EdgeHandle],
) -> Vec<(EdgeHandle, Option<FaceHandle>, Option<FaceHandle>)> {
    descendants
        .iter()
        .filter_map(|&de| {
            let e = brep.edges.get(de)?;
            let he = brep.half_edges.get(e.half_edge)?;
            let tw = brep.half_edges.get(he.twin)?;
            Some((
                de,
                brep.loops.get(he.loop_).map(|l| l.face),
                brep.loops.get(tw.loop_).map(|l| l.face),
            ))
        })
        .collect()
}

/// Given two candidate descendant faces `f_a`/`f_b`,
/// pick the one that descends from `desired_parent`.
///
/// "Descends" is checked by surface identity: sub-faces created
/// by `split_face` inherit the surface of the original face, so
/// if `desired_parent`'s surface in the original BRep matches
/// `f_a`'s current surface, `f_a` is the descendant we want.
/// Falls back to `f_a` if neither matches.
pub(super) fn choose_side_face(
    result: &BRep,
    original: &BRep,
    desired_parent: FaceHandle,
    f_a: FaceHandle,
    f_b: FaceHandle,
) -> FaceHandle {
    if f_a == desired_parent {
        return f_a;
    }
    if f_b == desired_parent {
        return f_b;
    }
    let Some(parent_face) = original.faces.get(desired_parent) else {
        return f_a;
    };
    let Some(parent_surf) = parent_face.surface.as_ref() else {
        return f_a;
    };
    let surf_matches = |fh: FaceHandle| -> bool {
        let Some(face) = result.faces.get(fh) else {
            return false;
        };
        let Some(sof) = face.surface.as_ref() else {
            return false;
        };
        match (&parent_surf.geometry, &sof.geometry) {
            (
                SurfaceGeometry::Plane {
                    origin: o1,
                    normal: n1,
                    ..
                },
                SurfaceGeometry::Plane {
                    origin: o2,
                    normal: n2,
                    ..
                },
            ) => {
                let dot_n = n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
                if dot_n.abs() < 0.999 {
                    return false;
                }
                let dx = o2.x - o1.x;
                let dy = o2.y - o1.y;
                let dz = o2.z - o1.z;
                let off = dx * n1.x + dy * n1.y + dz * n1.z;
                off.abs() < 1e-6
            }
            _ => false,
        }
    };
    if surf_matches(f_a) {
        f_a
    } else if surf_matches(f_b) {
        f_b
    } else {
        f_a
    }
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
