//! The runs of a contact curve inside faces' trims, and the surfaces
//! that touch along a line or at a point.

use std::collections::BTreeMap;

use super::{Box3, TouchPlan, at};
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, EdgeHandle, FaceHandle};

/// Samples along a curve when clipping it to trims.
const SAMPLES: usize = 32;

impl BRep {
    /// The runs of `edge` lying on `face` inside its trim: none unless
    /// every sample of the edge lies on the face's surface.
    pub(in crate::euler::nonmanifold) fn edge_runs_on_face(
        &self,
        edge: EdgeHandle,
        face: FaceHandle,
    ) -> Vec<CurveOnEdge> {
        let tol = self.tolerance.linear();
        if self.face_edges(face).contains(&edge) {
            return Vec::new();
        }
        let Some(curve) = self.edges.get(edge).and_then(|e| e.curve.clone()) else {
            return Vec::new();
        };
        let on_surface = (0..=16).all(|k| {
            let p = curve.geometry.point_at(at(&curve, k, 16));
            self.face_point_trim_containment(face, &p)
                .is_some_and(|c| c.distance <= tol)
        });
        if !on_surface {
            return Vec::new();
        }
        self.runs_inside(&curve, &[face])
            .into_iter()
            .filter(|run| self.face_interior_holds(face, &mid(run), tol))
            .collect()
    }

    /// The runs of `curve` inside every one of `faces`, each end where
    /// it leaves one moved onto the boundary it crosses.
    fn runs_inside(&self, curve: &CurveOnEdge, faces: &[FaceHandle]) -> Vec<CurveOnEdge> {
        let tol = self.tolerance.linear();
        let inside = |t: f64| {
            let p = curve.geometry.point_at(t);
            faces.iter().all(|&f| self.face_holds(f, &p, tol))
        };
        let ts: Vec<f64> = (0..=SAMPLES).map(|k| at(curve, k, SAMPLES)).collect();
        let flags: Vec<bool> = ts.iter().map(|&t| inside(t)).collect();
        let step = (curve.t_end - curve.t_start).abs() / SAMPLES as f64;
        let mut runs = Vec::new();
        let mut start = flags[0].then_some(ts[0]);
        for k in 1..ts.len() {
            if flags[k] == flags[k - 1] {
                continue;
            }
            let (mut a, mut b) = (ts[k - 1], ts[k]);
            for _ in 0..50 {
                let m = 0.5 * (a + b);
                if inside(m) == flags[k - 1] {
                    a = m;
                } else {
                    b = m;
                }
            }
            let t = self.snap_to_boundary(curve, 0.5 * (a + b), step, faces);
            if flags[k] {
                start = Some(t);
            } else if let Some(s) = start.take() {
                runs.push(CurveOnEdge {
                    geometry: curve.geometry.clone(),
                    t_start: s,
                    t_end: t,
                });
            }
        }
        if let Some(s) = start {
            runs.push(CurveOnEdge {
                geometry: curve.geometry.clone(),
                t_start: s,
                t_end: curve.t_end,
            });
        }
        runs
    }

    /// The parameter near `t` (a trim crossing found by bisection,
    /// within about the linear tolerance of the boundary) where `curve`
    /// meets the boundary of `faces` exactly: the nearest vertex within
    /// a few tolerances, else the nearest edge the curve crosses there
    /// (one it runs along on both sides is not crossed), else `t`.
    fn snap_to_boundary(
        &self,
        curve: &CurveOnEdge,
        t: f64,
        step: f64,
        faces: &[FaceHandle],
    ) -> f64 {
        let tol = self.tolerance.linear();
        let reach = 4.0 * tol;
        let point = |s: f64| curve.geometry.point_at(s);
        let h = step * 1e-6;
        let speed = (point(t + h).distance_to(&point(t - h)) / (2.0 * h)).max(1e-300);
        // The parameter window a few tolerances either side of `t`.
        let w = (reach / speed).min(step);
        let (a, b) = (t - w, t + w);
        let here = point(t);
        let vertex = faces
            .iter()
            .flat_map(|&f| self.face_vertices(f))
            .map(|v| self.vertices[v].position)
            .filter(|q| q.distance_to(&here) <= reach)
            .min_by(|p, q| p.distance_to(&here).total_cmp(&q.distance_to(&here)));
        if let Some(q) = vertex {
            return golden_min(|s| point(s).distance_to(&q), a, b);
        }
        let along = 25.0 * w;
        let mut best = (f64::INFINITY, t);
        for &f in faces {
            for e in self.face_edges(f) {
                let gap = |s: f64| self.edge_gap(e, &point(s));
                if gap(t) > reach || (gap(t - along) <= tol && gap(t + along) <= tol) {
                    continue;
                }
                let found = golden_min(gap, a, b);
                if gap(found) < best.0 {
                    best = (gap(found), found);
                }
            }
        }
        if best.0 <= tol { best.1 } else { t }
    }

    /// The tangent ruling or point where `f` and `g` touch, its runs
    /// inside both faces added to `plan`.
    pub(in crate::euler::nonmanifold) fn tangent_contacts(
        &self,
        f: FaceHandle,
        g: FaceHandle,
        boxes: &BTreeMap<FaceHandle, Box3>,
        plan: &mut TouchPlan,
    ) {
        let tol = self.tolerance.linear();
        let ang = self.tolerance.angular();
        let surface = |h: FaceHandle| self.faces[h].surface.as_ref().map(|s| &s.geometry);
        let (Some(sf), Some(sg)) = (surface(f), surface(g)) else {
            return;
        };
        if let Some(p) = tangent_point(sf, sg, tol) {
            if self.face_holds(f, &p, tol) && self.face_holds(g, &p, tol) {
                plan.add_point(p, tol);
            }
            return;
        }
        let Some((origin, direction)) = tangent_line(sf, sg, tol, ang) else {
            return;
        };
        let span = |b: &Box3| {
            let [lo, hi] = b.corners();
            let mut range = (f64::INFINITY, f64::NEG_INFINITY);
            for x in [lo.x, hi.x] {
                for y in [lo.y, hi.y] {
                    for z in [lo.z, hi.z] {
                        let s = (Point3::new(x, y, z) - origin).dot(&direction);
                        range = (range.0.min(s), range.1.max(s));
                    }
                }
            }
            range
        };
        let (a, b) = (span(&boxes[&f]), span(&boxes[&g]));
        let (s0, s1) = (a.0.max(b.0), a.1.min(b.1));
        if s1 - s0 <= tol {
            return;
        }
        let line = CurveOnEdge {
            geometry: CurveGeometry::Line { origin, direction },
            t_start: s0,
            t_end: s1,
        };
        for run in self.runs_inside(&line, &[f, g]) {
            let m = mid(&run);
            if self.face_interior_holds(f, &m, tol) || self.face_interior_holds(g, &m, tol) {
                plan.add_curve(run, tol);
            }
        }
    }
}

/// The point halfway along `run`.
fn mid(run: &CurveOnEdge) -> Point3 {
    run.geometry.point_at(0.5 * (run.t_start + run.t_end))
}

/// The minimum of `f` over `[a, b]` by golden section.
pub(super) fn golden_min(f: impl Fn(f64) -> f64, mut a: f64, mut b: f64) -> f64 {
    let r = 0.5 * (5.0_f64.sqrt() - 1.0);
    for _ in 0..80 {
        let (c, d) = (b - r * (b - a), a + r * (b - a));
        if f(c) < f(d) {
            b = d;
        } else {
            a = c;
        }
    }
    0.5 * (a + b)
}

/// The line two surfaces touch along: a plane and a cylinder whose axis
/// lies parallel to the plane at its radius, two cylinders with
/// parallel axes their radii's sum or difference apart.
fn tangent_line(
    f: &SurfaceGeometry,
    g: &SurfaceGeometry,
    tol: f64,
    ang: f64,
) -> Option<(Point3, Vec3)> {
    use SurfaceGeometry::{Cylinder, Plane};
    match (f, g) {
        (
            Plane { origin, normal, .. },
            Cylinder {
                origin: c,
                axis,
                radius,
            },
        )
        | (
            Cylinder {
                origin: c,
                axis,
                radius,
            },
            Plane { origin, normal, .. },
        ) => {
            let d = (*c - *origin).dot(normal);
            (axis.dot(normal).abs() <= ang && (d.abs() - radius).abs() <= tol)
                .then(|| (*c - *normal * d, *axis))
        }
        (
            Cylinder {
                origin: c1,
                axis: a1,
                radius: r1,
            },
            Cylinder {
                origin: c2,
                axis: a2,
                radius: r2,
            },
        ) => {
            if a1.cross(a2).length() > ang {
                return None;
            }
            let w = *c2 - *c1;
            let across = w - *a1 * w.dot(a1);
            let gap = across.length();
            let u = across.normalized()?;
            if (gap - (r1 + r2)).abs() <= tol {
                Some((*c1 + u * *r1, *a1))
            } else if gap > tol && (gap - (r1 - r2).abs()).abs() <= tol {
                // The smaller inside the larger, touching on the side
                // away from the larger's axis.
                let p = if r1 > r2 {
                    *c1 + u * *r1
                } else {
                    *c2 - u * *r2
                };
                Some((p, *a1))
            } else {
                None
            }
        }
        _ => None,
    }
}

/// The point two surfaces touch at: a sphere on a plane, two spheres
/// their radii's sum or difference apart.
fn tangent_point(f: &SurfaceGeometry, g: &SurfaceGeometry, tol: f64) -> Option<Point3> {
    use SurfaceGeometry::{Plane, Sphere};
    match (f, g) {
        (Plane { origin, normal, .. }, Sphere { center, radius, .. })
        | (Sphere { center, radius, .. }, Plane { origin, normal, .. }) => {
            let d = (*center - *origin).dot(normal);
            ((d.abs() - radius).abs() <= tol).then(|| *center - *normal * d)
        }
        (
            Sphere {
                center: c1,
                radius: r1,
                ..
            },
            Sphere {
                center: c2,
                radius: r2,
                ..
            },
        ) => {
            let w = *c2 - *c1;
            let gap = w.length();
            let u = w.normalized()?;
            if (gap - (r1 + r2)).abs() <= tol {
                Some(*c1 + u * *r1)
            } else if gap > tol && (gap - (r1 - r2).abs()).abs() <= tol {
                Some(if r1 > r2 {
                    *c1 + u * *r1
                } else {
                    *c2 - u * *r2
                })
            } else {
                None
            }
        }
        _ => None,
    }
}
