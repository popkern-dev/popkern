//! Where an edge of a piece touches a face of the same piece at a point
//! (PK-1155): the vertex [`BRep::glue_touching`] gives both.
//!
//! A boolean whose result rests one part on another inside one shell
//! (a ball on a box's edge, both parts of one solid joined elsewhere)
//! touches itself where an edge meets a face it shares no vertex with.
//! Each edge of a face is cut at stations and measured against the
//! other face's surface; at a station nearer than both neighbours,
//! which lie farther than the tolerance (an edge lying along the face
//! is a curve contact, not a run of points), the
//! nearest point is located as the root of `C′ · N`, `N` the surface's
//! normal at the curve point's foot (a touch is a flat minimum: its
//! value pins the place only to the square root of the rounding; the
//! squared distance's derivative `C′ · (C − S)` is a triple root for a
//! line resting on a ball, pinned to the cube root, measured 2e-5 off;
//! `C′ · N` is a simple root), and a point within the tolerance whose
//! foot lies on the face is a contact, moved onto the edge's end
//! vertex when within the tolerance of it.

use super::{TouchPlan, at};
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, FaceHandle};

/// Stations per edge.
const STATIONS: usize = 32;

impl BRep {
    /// Points where an edge of `f` touches face `g`; at an edge's end
    /// only when `ends_too` (between pieces, a vertex resting on the other
    /// piece is a point contact at its exact place already).
    pub(super) fn edge_touch_points(
        &self,
        f: FaceHandle,
        g: FaceHandle,
        ends_too: bool,
        plan: &mut TouchPlan,
    ) {
        let tol = self.tolerance.linear();
        let Some(surface) = self.faces.get(g).and_then(|x| x.surface.as_ref()) else {
            return;
        };
        let foot = |p: &Point3| -> Option<(Point3, Vec3)> {
            let (u, v) = project_point_to_uv(&surface.geometry, p)?;
            let q = surface.geometry.point_at(u, v);
            // The normal there; at a chart pole, the way to the point.
            let n = surface.geometry.normal_at(u, v);
            let n = if n.length() > 0.5 { n } else { *p - q };
            Some((q, n))
        };
        for e in self.face_edges(f) {
            let Some(c) = self.edges.get(e).and_then(|x| x.curve.as_ref()) else {
                continue;
            };
            let ts: Vec<f64> = (0..=STATIONS).map(|k| at(c, k, STATIONS)).collect();
            let gap = |t: f64| {
                let p = c.geometry.point_at(t);
                foot(&p).map_or(f64::INFINITY, |(q, _)| p.distance_to(&q))
            };
            let gaps: Vec<f64> = ts.iter().map(|&t| gap(t)).collect();
            // The curve's rate along the surface normal at its foot: a
            // simple root at a touch, where the squared distance's own
            // derivative is a triple root for a line on a ball.
            let slope = |t: f64| -> Option<f64> {
                let (_, n) = foot(&c.geometry.point_at(t))?;
                Some(c.geometry.tangent_at(t).dot(&n.normalized()?))
            };
            let ends: Vec<Point3> = self
                .edge_vertices(e)
                .map(|(a, b)| vec![self.vertices[a].position, self.vertices[b].position])
                .unwrap_or_default();
            for k in 0..ts.len() {
                if !ends_too && (k == 0 || k + 1 == ts.len()) {
                    continue;
                }
                let left = if k == 0 { f64::INFINITY } else { gaps[k - 1] };
                let right = gaps.get(k + 1).copied().unwrap_or(f64::INFINITY);
                // An isolated minimum: the gap rises past the tolerance
                // at each neighbouring station (an edge lying along the
                // face is a curve contact, `edge_runs_on_face`).
                if !(gaps[k] <= left && gaps[k] <= right)
                    || gaps[k] > 100.0 * tol
                    || left <= tol
                    || right <= tol
                {
                    continue;
                }
                let lo = ts[k.saturating_sub(1)];
                let hi = ts[(k + 1).min(ts.len() - 1)];
                let t = root(&slope, lo, hi).unwrap_or(ts[k]);
                let p = c.geometry.point_at(t);
                if gap(t) <= tol && self.face_holds(g, &p, tol) {
                    // At an end: the edge's vertex resting on the face.
                    let p = ends
                        .iter()
                        .copied()
                        .find(|q| q.distance_to(&p) <= tol)
                        .unwrap_or(p);
                    plan.add_point(p, tol);
                }
            }
        }
    }
}

/// The root of `g` in `[a, b]` by bisection, when its ends differ in
/// sign (either way round).
fn root(g: &impl Fn(f64) -> Option<f64>, mut a: f64, mut b: f64) -> Option<f64> {
    let (mut ga, gb) = (g(a)?, g(b)?);
    if ga == 0.0 {
        return Some(a);
    }
    if gb == 0.0 {
        return Some(b);
    }
    if ga.signum() == gb.signum() {
        return None;
    }
    for _ in 0..100 {
        let m = 0.5 * (a + b);
        let gm = g(m)?;
        if gm == 0.0 || (b - a).abs() <= f64::EPSILON * m.abs().max(1.0) {
            return Some(m);
        }
        if gm.signum() == ga.signum() {
            (a, ga) = (m, gm);
        } else {
            b = m;
        }
    }
    Some(0.5 * (a + b))
}
