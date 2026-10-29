//! Faces a request's bevels consume or sweep past (PK-1006).
//!
//! **The class.** The one-pass pipeline trims each face by the rails
//! of its chamfered edges: every rail must keep a side of the face and
//! every other side must survive. When the rails cross past a side (a
//! cone's slant triangle chamfered on all three sides beyond its
//! inradius, a narrow face whose two rails meet below its base), the
//! face has no region left in that shape. The pipeline either declined
//! somewhere downstream by luck or cut the crossing rails pairwise and
//! served a validator-clean body with material left outside a bevel
//! plane (measured: a 5-sided cone, six edges at d = 0.59, served 0.2%
//! over the exact volume with a vertex 0.28 outside a bevel).
//!
//! **The gate.** On a convex planar face with straight sides and no
//! holes, the region the bevels keep is the face clipped by each rail's
//! half-plane (the side away from its edge). Each rail must keep a span
//! longer than `tol` inside the face and the other rails' half-planes,
//! and each unchamfered side a span longer than `tol` inside the rails'
//! half-planes. Otherwise the request declines with
//! [`BlendFailure::InconsistentGeometry`] naming the face, and bevel
//! growth (PK-1001) serves it from a smaller size, consuming the face
//! or flipping the side at its exact instant.
//!
//! **Scope.** Convex planar faces with straight sides and rails; other
//! faces are not gated here (a rail on a non-convex face acts only over
//! its own span, not as a half-plane). **Tolerance:** the caller's
//! linear tolerance × 10.

use crate::construct::EdgeGeometry;
use crate::error::BlendFailure;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};
use std::collections::BTreeMap;

/// A half-plane of the face's plane: kept where `m · (x − r) >= 0`.
type HalfPlane = (Vec3, Point3);

/// A face side: its edge and its ends in loop order.
type Side = (EdgeHandle, Point3, Point3);

/// A straight curve's ends, when every control point lies on the chord.
fn straight(curve: &NurbsCurve, tol: f64) -> Option<(Point3, Point3)> {
    let (a, b) = (
        *curve.control_points.first()?,
        *curve.control_points.last()?,
    );
    let u = (b - a).normalized()?;
    curve
        .control_points
        .iter()
        .all(|p| {
            let w = *p - a;
            (w - u * w.dot(&u)).length() <= tol
        })
        .then_some((a, b))
}

/// The span of the line `o + s·u` (unit `u`) inside every half-plane,
/// as a length (zero when empty).
fn span(o: Point3, u: Vec3, bounds: &[HalfPlane], tol: f64) -> f64 {
    let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
    for (m, r) in bounds {
        let (a, b) = (m.dot(&u), m.dot(&(o - *r)));
        if a.abs() <= 1e-12 {
            if b < -tol {
                return 0.0;
            }
        } else if a > 0.0 {
            lo = lo.max(-b / a);
        } else {
            hi = hi.min(-b / a);
        }
    }
    (hi - lo).max(0.0)
}

/// The half-plane through `p → q` (in the plane of normal `n`) that
/// keeps `inside`.
fn keeping(p: Point3, q: Point3, n: Vec3, inside: Point3) -> Option<HalfPlane> {
    let m = n.cross(&(q - p)).normalized()?;
    Some(if m.dot(&(inside - p)) >= 0.0 {
        (m, p)
    } else {
        (-m, p)
    })
}

/// The face's straight sides as `(edge, start, end)` in loop order,
/// when the face is planar, holeless and convex.
fn convex_sides(brep: &BRep, face: FaceHandle, tol: f64) -> Option<(Vec3, Point3, Vec<Side>)> {
    let f = brep.faces.get(face)?;
    if !f.inner_loops.is_empty()
        || !matches!(
            f.surface.as_ref()?.geometry,
            popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
        )
    {
        return None;
    }
    let start = brep.loops.get(f.outer_loop)?.half_edge;
    let mut sides = Vec::new();
    let mut he = start;
    for _ in 0..=brep.half_edges.len() {
        let h = brep.half_edges.get(he)?;
        let p = brep.vertices.get(h.origin)?.position;
        let q = brep
            .vertices
            .get(brep.half_edges.get(h.twin)?.origin)?
            .position;
        let c = brep.edges.get(h.edge)?.curve.as_ref()?;
        straight(&c.geometry.to_nurbs_in_range((c.t_start, c.t_end)), tol)?;
        sides.push((h.edge, p, q));
        he = h.next;
        if he == start {
            break;
        }
    }
    if sides.len() < 3 {
        return None;
    }
    // Newell normal and centroid; every turn must agree with it.
    let mut n = Vec3::new(0.0, 0.0, 0.0);
    let mut c = Vec3::new(0.0, 0.0, 0.0);
    let origin = Point3::new(0.0, 0.0, 0.0);
    for &(_, p, q) in &sides {
        n = n + (p - origin).cross(&(q - origin));
        c = c + (p - origin);
    }
    let n = n.normalized()?;
    let centre = origin + c * (1.0 / sides.len() as f64);
    for (i, &(_, p, q)) in sides.iter().enumerate() {
        let r = sides[(i + 1) % sides.len()].2;
        if (q - p).cross(&(r - q)).dot(&n) < -tol * (q - p).length().max((r - q).length()) {
            return None;
        }
    }
    Some((n, centre, sides))
}

/// Why the request's rails leave some face without the region the
/// one-pass pipeline builds, if they do.
pub(crate) fn consumed_face(
    brep: &BRep,
    geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Option<BlendFailure> {
    let tol = (tolerance.linear() * 10.0).max(1e-9);
    let mut rails: BTreeMap<FaceHandle, Vec<(EdgeHandle, &NurbsCurve)>> = BTreeMap::new();
    for g in geometries {
        for (face, rail) in [
            (g.analysis.face_1, &g.trim_3d_1),
            (g.analysis.face_2, &g.trim_3d_2),
        ] {
            rails.entry(face).or_default().push((g.analysis.edge, rail));
        }
    }
    for (face, on_face) in rails {
        let Some((n, centre, sides)) = convex_sides(brep, face, tol) else {
            continue;
        };
        let inner: Vec<HalfPlane> = sides
            .iter()
            .filter_map(|&(_, p, q)| keeping(p, q, n, centre))
            .collect();
        // Each rail's half-plane keeps the side away from its edge.
        let line = |&(edge, rail): &(EdgeHandle, &NurbsCurve)| {
            let (a, b) = straight(rail, tol)?;
            let &(_, p, q) = sides.iter().find(|s| s.0 == edge)?;
            let (m, r) = keeping(a, b, n, p + (q - p) * 0.5)?;
            Some((edge, a, (b - a).normalized()?, (-m, r)))
        };
        let Some(lines) = on_face.iter().map(line).collect::<Option<Vec<_>>>() else {
            continue;
        };
        let fail = |what: String, at: Point3| {
            Some(BlendFailure::InconsistentGeometry {
                reason: format!(
                    "chamfer bevels consume {face:?}: {what} keeps no length (grown from a \
                     smaller size instead)"
                ),
                location: Some(at),
            })
        };
        for (k, &(edge, a, u, _)) in lines.iter().enumerate() {
            let bounds: Vec<HalfPlane> = inner
                .iter()
                .copied()
                .chain(
                    lines
                        .iter()
                        .enumerate()
                        .filter(|&(j, _)| j != k)
                        .map(|(_, l)| l.3),
                )
                .collect();
            if span(a, u, &bounds, tol) <= tol {
                return fail(format!("the rail of {edge:?}"), a);
            }
        }
        let cuts: Vec<HalfPlane> = lines.iter().map(|l| l.3).collect();
        for &(edge, p, q) in sides.iter().filter(|s| lines.iter().all(|l| l.0 != s.0)) {
            let Some(u) = (q - p).normalized() else {
                continue;
            };
            let ends = [(u, p), (-u, q)];
            let bounds: Vec<HalfPlane> = cuts.iter().copied().chain(ends).collect();
            if span(p, u, &bounds, tol) <= tol {
                return fail(format!("its side {edge:?}"), p + (q - p) * 0.5);
            }
        }
    }
    None
}
