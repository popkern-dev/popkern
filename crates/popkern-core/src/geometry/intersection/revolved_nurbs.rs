//! Exact intersections of a NURBS surface of revolution with a coaxial
//! plane, cylinder, cone or sphere (PK-1126).
//!
//! A revolve builds its faces as NURBS surfaces of revolution: the
//! profile in `u`, rational quadratic arcs in `v`. Against a plane
//! perpendicular to the axis, or a cylinder, cone or sphere on it, the
//! intersection is a set of circles, one wherever the profile meets
//! the other surface's meridian: the roots of one function `g(u)` along
//! the profile. The pairs went to the marcher, which near a chart pole
//! or where both seams coincide returned records with holes (a ball's
//! equator on a revolved cap came as a circle and two arcs, PK-1128;
//! a paraboloid's tip cut by a plane 1e-4 above it disagreed with its
//! trims by 3.9e-4).
//!
//! Each circle is the surface's own iso-curve at the root, exact and on
//! its chart the line `u = u*`; the other surface's pcurve is its exact
//! affine image on a plane or the co-parameterised chart fit. A touch
//! (`|g|` within the tolerance without a sign change) and a root on the
//! axis (a point) emit nothing, as every arm does: contact discovery
//! types them. A stretch of the profile lying in the other surface (a
//! revolved radial line in its plane) is a coincidence, left to the
//! coincidence detector. A surface the recognition cannot certify, or
//! a pair off the axis, keeps the marching path.

use crate::geometry::curve::CurveEval;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::circle_of_revolution::pcurve_on;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, PortableMath, Vec3};
use crate::tolerance::Tolerance;

/// A parameter window `((u0, u1), (v0, v1))`.
type Window = ((f64, f64), (f64, f64));

/// The revolved NURBS × coaxial analytic answer, in argument order, or
/// `None` when the pair is not of that class.
pub(crate) fn try_revolved_nurbs(
    surface1: &SurfaceGeometry,
    w1: Window,
    surface2: &SurfaceGeometry,
    w2: Window,
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (nurbs, wn, other, wo, swapped) = match (surface1, surface2) {
        (SurfaceGeometry::Nurbs(n), o) if !matches!(o, SurfaceGeometry::Nurbs(_)) => {
            (n, w1, o, w2, false)
        }
        (o, SurfaceGeometry::Nurbs(n)) if !matches!(o, SurfaceGeometry::Nurbs(_)) => {
            (n, w2, o, w1, true)
        }
        _ => return None,
    };
    // The angle in v (the profile in u), or in u: a revolve authors
    // the angle in u (PK-455), read here through the transpose.
    let transposed;
    let (view, profile_window, angle_in_u) = match revolution_axis(nurbs, tolerance) {
        Some(_) => (nurbs, wn.0, false),
        None => {
            transposed = transpose(nurbs);
            (&transposed, wn.1, true)
        }
    };
    let axis = revolution_axis(view, tolerance)?;
    let meridian = Meridian::of(other, axis, tolerance)?;
    let roots = profile_roots(view, profile_window, axis, &meridian, tolerance)?;
    let (v0, v1) = domain(&view.knots_v, view.degree_v);
    let mut out = Vec::with_capacity(roots.len());
    for u in roots {
        let curve_3d = view.iso_curve_at_u(u);
        let ends = [(u, v0), (u, v1)].map(|(p, a)| {
            if angle_in_u {
                Point3::new(a, p, 0.0)
            } else {
                Point3::new(p, a, 0.0)
            }
        });
        let on_nurbs = NurbsCurve {
            degree: 1,
            knots: vec![v0, v0, v1, v1],
            control_points: ends.to_vec(),
            weights: vec![1.0, 1.0],
        };
        let on_other = pcurve_on(other, &curve_3d, wo.0, wo.1, tolerance)?;
        let is_closed =
            curve_3d.point_at(v0).distance_to(&curve_3d.point_at(v1)) <= tolerance.linear();
        let (curve_uv1, curve_uv2) = if swapped {
            (on_other, on_nurbs)
        } else {
            (on_nurbs, on_other)
        };
        out.push(SurfaceSurfaceIntersection::Curve(CurveOnSurface {
            curve_3d,
            curve_uv1,
            curve_uv2,
            is_closed,
        }));
    }
    Some(out)
}

/// The valid parameter domain of a knot vector of `degree`: from
/// `knots[degree]` to `knots[len − 1 − degree]` (a periodic-form ring's
/// vector runs past it, PK-890).
fn domain(knots: &[f64], degree: usize) -> (f64, f64) {
    (knots[degree], knots[knots.len() - 1 - degree])
}

/// `s` with its parameter directions exchanged.
fn transpose(s: &NurbsSurface) -> NurbsSurface {
    let mut control_points = Vec::with_capacity(s.control_points.len());
    let mut weights = Vec::with_capacity(s.weights.len());
    for j in 0..s.num_v {
        for i in 0..s.num_u {
            control_points.push(s.control_points[i * s.num_v + j]);
            weights.push(s.weights[i * s.num_v + j]);
        }
    }
    NurbsSurface {
        degree_u: s.degree_v,
        degree_v: s.degree_u,
        knots_u: s.knots_v.clone(),
        knots_v: s.knots_u.clone(),
        control_points,
        num_u: s.num_v,
        num_v: s.num_u,
        weights,
    }
}

/// A line: a point on it and its unit direction.
#[derive(Clone, Copy)]
pub(crate) struct Axis {
    origin: Point3,
    dir: Vec3,
}

impl Axis {
    /// `(z, r)`: the axial station and the distance from the axis.
    fn meridian_of(&self, p: Point3) -> (f64, f64) {
        let d = p - self.origin;
        let z = d.dot(&self.dir);
        (z, (d - self.dir * z).length())
    }
}

/// The axis `s` revolves about, when it is a surface of revolution in
/// `v` (rational quadratic arcs with separable weights) and sampled
/// points across it reproduce their profile's meridian coordinates.
pub(crate) fn revolution_axis(s: &NurbsSurface, tolerance: &Tolerance) -> Option<Axis> {
    let lin = tolerance.linear();
    if s.degree_v != 2 || s.num_v < 5 || s.num_u < 2 {
        return None;
    }
    let at = |i: usize, j: usize| s.control_points[i * s.num_v + j];
    let wt = |i: usize, j: usize| s.weights[i * s.num_v + j];
    for i in 0..s.num_u {
        for j in 0..s.num_v {
            let ratio = wt(i, j) / wt(i, 0);
            let first = wt(0, j) / wt(0, 0);
            if !ratio.is_finite() || (ratio - first).abs() > 1e-9 * first.abs().max(1.0) {
                return None;
            }
        }
    }
    // The on-circle control points carry the arc's full weight (the
    // arc middles carry cos(φ/2)); the row farthest from the axis fixes
    // it from three of them.
    let full = (0..s.num_v)
        .map(|j| wt(0, j) / wt(0, 0))
        .fold(f64::MIN, f64::max);
    let on_circle: Vec<usize> = (0..s.num_v)
        .filter(|&j| (wt(0, j) / wt(0, 0) - full).abs() <= 1e-9 * full)
        .collect();
    let [j0, j1, j2, ..] = on_circle[..] else {
        return None;
    };
    let spread = |i: usize| at(i, j0).distance_to(&at(i, j1));
    let row = (0..s.num_u).max_by(|&a, &b| spread(a).total_cmp(&spread(b)))?;
    if spread(row) <= lin {
        return None;
    }
    let (a, b, c) = (at(row, j0), at(row, j1), at(row, j2));
    let (ab, ac) = (b - a, c - a);
    let n = ab.cross(&ac);
    let n_sq = n.dot(&n);
    if n_sq <= 1e-300 {
        return None;
    }
    // Circumcentre of a, b, c.
    let centre = a + (n.cross(&ab) * ac.dot(&ac) + ac.cross(&n) * ab.dot(&ab)) * (0.5 / n_sq);
    let axis = Axis {
        origin: centre,
        dir: n.normalized()?,
    };
    // Certify by sampling: every point's (z, r) is its column's.
    let (u0, u1) = domain(&s.knots_u, s.degree_u);
    let (v0, v1) = domain(&s.knots_v, s.degree_v);
    for ku in 0..=6 {
        let u = u0 + (u1 - u0) * f64::from(ku) / 6.0;
        let (z0, r0) = axis.meridian_of(s.point_at(u, v0));
        for kv in 1..=5 {
            let v = v0 + (v1 - v0) * f64::from(kv) / 5.0;
            let (z, r) = axis.meridian_of(s.point_at(u, v));
            if (z - z0).abs() > lin || (r - r0).abs() > lin {
                return None;
            }
        }
    }
    Some(axis)
}

/// The other surface's meridian as a signed distance in the `(z, r)`
/// half-plane about the shared axis.
enum Meridian {
    /// `z = h`.
    Plane { h: f64 },
    /// `r = radius`.
    Cylinder { radius: f64 },
    /// `r = |z − apex| tan α`, weighted by `cos α` to a distance.
    Cone { apex: f64, tan: f64, cos: f64 },
    /// `(z − centre)² + r² = radius²`.
    Sphere { centre: f64, radius: f64 },
}

impl Meridian {
    /// The meridian of `surface` about `axis`, when the two are coaxial.
    fn of(surface: &SurfaceGeometry, axis: Axis, tolerance: &Tolerance) -> Option<Self> {
        let (lin, ang) = (tolerance.linear(), tolerance.angular());
        let parallel = |d: &Vec3| {
            d.normalized()
                .is_some_and(|d| d.cross(&axis.dir).length() <= ang)
        };
        let on_axis = |p: &Point3| axis.meridian_of(*p).1 <= lin;
        match surface {
            SurfaceGeometry::Plane { origin, normal, .. } if parallel(normal) => {
                Some(Self::Plane {
                    h: axis.meridian_of(*origin).0,
                })
            }
            SurfaceGeometry::Cylinder {
                origin,
                axis: a,
                radius,
            } if parallel(a) && on_axis(origin) => Some(Self::Cylinder { radius: *radius }),
            SurfaceGeometry::Cone {
                apex,
                axis: a,
                half_angle,
            } if parallel(a) && on_axis(apex) => {
                let (sin, cos) = (half_angle.psin(), half_angle.pcos());
                (cos.abs() > ang && sin.abs() > ang).then(|| Self::Cone {
                    apex: axis.meridian_of(*apex).0,
                    tan: (sin / cos).abs(),
                    cos: cos.abs(),
                })
            }
            SurfaceGeometry::Sphere { center, radius, .. } if on_axis(center) => {
                Some(Self::Sphere {
                    centre: axis.meridian_of(*center).0,
                    radius: *radius,
                })
            }
            _ => None,
        }
    }

    /// The signed distance of `(z, r)` from the meridian.
    fn gap(&self, (z, r): (f64, f64)) -> f64 {
        match *self {
            Self::Plane { h } => z - h,
            Self::Cylinder { radius } => r - radius,
            Self::Cone { apex, tan, cos } => (r - (z - apex).abs() * tan) * cos,
            Self::Sphere { centre, radius } => {
                ((z - centre) * (z - centre) + r * r).sqrt() - radius
            }
        }
    }
}

/// The profile parameters, within `u_window`, where the profile crosses
/// the meridian off the axis; `None` when a stretch of it lies in the
/// meridian (a coincidence).
fn profile_roots(
    s: &NurbsSurface,
    u_window: (f64, f64),
    axis: Axis,
    meridian: &Meridian,
    tolerance: &Tolerance,
) -> Option<Vec<f64>> {
    let lin = tolerance.linear();
    let (k0, k1) = domain(&s.knots_u, s.degree_u);
    let (lo, hi) = (
        u_window.0.min(u_window.1).max(k0),
        u_window.0.max(u_window.1).min(k1),
    );
    if hi <= lo {
        return Some(Vec::new());
    }
    let v0 = domain(&s.knots_v, s.degree_v).0;
    let g = |u: f64| meridian.gap(axis.meridian_of(s.point_at(u, v0)));
    let n = (32 * (s.num_u - s.degree_u).max(1)).clamp(64, 4096);
    let us: Vec<f64> = (0..=n)
        .map(|k| lo + (hi - lo) * k as f64 / n as f64)
        .collect();
    let gs: Vec<f64> = us.iter().map(|&u| g(u)).collect();
    // Three samples in a row on the meridian: the profile lies in it.
    if gs.windows(3).any(|w| w.iter().all(|x| x.abs() <= lin)) {
        return None;
    }
    let mut roots = Vec::new();
    for k in 0..n {
        let (ga, gb) = (gs[k], gs[k + 1]);
        if ga == 0.0 {
            roots.push(us[k]);
            continue;
        }
        if ga * gb > 0.0 {
            continue;
        }
        let (mut a, mut b, mut g_a) = (us[k], us[k + 1], ga);
        for _ in 0..100 {
            let m = 0.5 * (a + b);
            let gm = g(m);
            if gm * g_a <= 0.0 {
                b = m;
            } else {
                a = m;
                g_a = gm;
            }
        }
        roots.push(0.5 * (a + b));
    }
    if gs[n] == 0.0 {
        roots.push(us[n]);
    }
    // Off the axis only: a root there is a point, discovery's to type.
    roots.retain(|&u| axis.meridian_of(s.point_at(u, v0)).1 > lin);
    roots.dedup_by(|a, b| (*a - *b).abs() <= 1e-12 * (hi - lo).max(1.0));
    Some(roots)
}

#[cfg(test)]
#[path = "revolved_nurbs_tests.rs"]
mod tests;
