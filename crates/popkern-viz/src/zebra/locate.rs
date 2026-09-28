//! Locating the surface point a zebra pixel shows (PK-503).
//!
//! Interpolated `(u, v)` alone is not the pixel's surface point: on a
//! coarse facet the chart's parametric distortion displaces it by a
//! fraction of the facet, which reads as a stripe shift that belongs
//! to the tessellation. The ray hit is the point the eye sees, so the
//! surface view is tessellation-independent up to which facet claimed
//! the pixel.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_mesh::MeshVertex;

use super::{SurfaceHit, dot, unit};

/// Newton iterations for the ray and closest-point refinements.
const REFINE_ITERATIONS: usize = 24;

/// Multiple of the facet reach allowed for a ray root in front of the
/// facet (see [`PixelQuery::ray_hit`]).
const SILHOUETTE_REACH: f64 = 4.0;

/// Least squares `a·x + b·y ≈ r` (3 equations, 2 unknowns) by the
/// normal equations; `None` when `a` and `b` are (relatively) parallel.
fn lsq2(a: [f64; 3], b: [f64; 3], r: [f64; 3]) -> Option<(f64, f64)> {
    let (aa, ab, bb) = (dot(a, a), dot(a, b), dot(b, b));
    let det = aa * bb - ab * ab;
    // NaN-safe: a NaN determinant fails the comparison and declines.
    if det.is_nan() || det <= 1e-24 * aa * bb {
        return None;
    }
    let (ar, br) = (dot(a, r), dot(b, r));
    Some(((bb * ar - ab * br) / det, (aa * br - ab * ar) / det))
}

/// Solve the 3×3 system `[a b c]·x = r` by Cramer's rule; `None` when
/// the columns are (relatively) dependent.
fn solve3(a: [f64; 3], b: [f64; 3], c: [f64; 3], r: [f64; 3]) -> Option<[f64; 3]> {
    let cross = |x: [f64; 3], y: [f64; 3]| {
        [
            x[1] * y[2] - x[2] * y[1],
            x[2] * y[0] - x[0] * y[2],
            x[0] * y[1] - x[1] * y[0],
        ]
    };
    let bc = cross(b, c);
    let det = dot(a, bc);
    let scale = dot(a, a).sqrt() * dot(b, b).sqrt() * dot(c, c).sqrt();
    if det.is_nan() || det.abs() <= 1e-12 * scale {
        return None;
    }
    Some([
        dot(r, bc) / det,
        dot(a, cross(r, c)) / det,
        dot(a, cross(b, r)) / det,
    ])
}

fn to3(v: popkern_core::geometry::Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn point3(surface: &SurfaceGeometry, uv: (f64, f64)) -> [f64; 3] {
    let q = surface.point_at(uv.0, uv.1);
    [q.x, q.y, q.z]
}

/// Clamp non-periodic parameters into the surface's natural range.
fn clamp_uv(surface: &SurfaceGeometry, uv: (f64, f64)) -> (f64, f64) {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let (pu, pv) = surface.parameter_periodicity();
    let c = |x: f64, lo: f64, hi: f64, periodic: bool| {
        if periodic || !(lo.is_finite() && hi.is_finite()) {
            x
        } else {
            x.clamp(lo, hi)
        }
    };
    (c(uv.0, u0, u1, pu.is_some()), c(uv.1, v0, v1, pv.is_some()))
}

/// Scale a parameter step so it stops short of the natural range of
/// any non-periodic coordinate (fraction-to-boundary rule): clamping
/// instead would stall the iteration against a chart bound, which is
/// exactly where sphere/cone poles and collapsed NURBS edges sit.
fn damped_step(surface: &SurfaceGeometry, uv: (f64, f64), du: f64, dv: f64) -> (f64, f64) {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let (pu, pv) = surface.parameter_periodicity();
    let limit = |x: f64, dx: f64, lo: f64, hi: f64, periodic: bool| -> f64 {
        if periodic || dx == 0.0 {
            return 1.0;
        }
        let room = if dx > 0.0 { hi - x } else { x - lo };
        if !room.is_finite() {
            return 1.0;
        }
        (0.9 * room.max(0.0) / dx.abs()).min(1.0)
    };
    let alpha = limit(uv.0, du, u0, u1, pu.is_some()).min(limit(uv.1, dv, v0, v1, pv.is_some()));
    (uv.0 + alpha * du, uv.1 + alpha * dv)
}

/// The query one pixel poses: its facet point, the eye, and how far
/// from the facet an accepted surface point may lie.
pub(super) struct PixelQuery {
    /// Eye position.
    pub eye: [f64; 3],
    /// Perspective-correct facet point on the pixel's ray.
    pub facet: [f64; 3],
    /// Proximity bound (the claiming triangle's longest edge).
    pub reach: f64,
    /// Facing sign of the face's surface chart (`±1`): the visible
    /// point's facing normal is `sense · (S_u × S_v)`.
    pub sense: f64,
}

impl PixelQuery {
    fn near(&self, surface: &SurfaceGeometry, uv: (f64, f64), reach: f64) -> bool {
        let q = point3(surface, uv);
        let r = [
            q[0] - self.facet[0],
            q[1] - self.facet[1],
            q[2] - self.facet[2],
        ];
        dot(r, r).sqrt() <= reach
    }

    /// Whether the surface at `uv` faces the eye along `d`: the eye
    /// sees a surface point only from its outward side, so a ray
    /// root that faces away is the ray's EXIT through the body.
    fn faces_eye(&self, surface: &SurfaceGeometry, uv: (f64, f64), d: [f64; 3]) -> bool {
        popkern_mesh::surface_limit_normal(surface, uv.0, uv.1)
            .is_some_and(|n| self.sense * (n.x * d[0] + n.y * d[1] + n.z * d[2]) < 0.0)
    }

    /// Newton on `F(u, v, t) = S(u, v) − (eye + t·d) = 0` from `seed`;
    /// returns the root and its ray parameter.
    fn ray_hit(&self, surface: &SurfaceGeometry, seed: (f64, f64)) -> Option<((f64, f64), f64)> {
        let e = self.eye;
        let to_facet = [
            self.facet[0] - e[0],
            self.facet[1] - e[1],
            self.facet[2] - e[2],
        ];
        let dist = dot(to_facet, to_facet).sqrt();
        let d = unit(to_facet)?;
        let eps = 1e-12 * (1.0 + dist);
        let (mut uv, mut t) = (seed, dist);
        for _ in 0..REFINE_ITERATIONS {
            let p = point3(surface, uv);
            let f = [
                p[0] - e[0] - t * d[0],
                p[1] - e[1] - t * d[1],
                p[2] - e[2] - t * d[2],
            ];
            if dot(f, f).sqrt() < eps {
                let uv = clamp_uv(surface, uv);
                // A root in FRONT of the facet may sit farther from it
                // than the facet's size: near a silhouette an inscribed
                // facet lies deep inside a convex surface and its ray
                // crosses the surface at a grazing angle. Behind the
                // facet, only the facet's own neighbourhood counts.
                let reach = if t <= dist {
                    SILHOUETTE_REACH * self.reach
                } else {
                    self.reach
                };
                let visible =
                    t > 0.0 && self.near(surface, uv, reach) && self.faces_eye(surface, uv, d);
                return visible.then_some((uv, t));
            }
            let su = to3(surface.partial_u(uv.0, uv.1));
            let sv = to3(surface.partial_v(uv.0, uv.1));
            let step = solve3(su, sv, [-d[0], -d[1], -d[2]], [-f[0], -f[1], -f[2]])?;
            uv = damped_step(surface, clamp_uv(surface, uv), step[0], step[1]);
            t += step[2];
            if !(uv.0.is_finite() && uv.1.is_finite() && t.is_finite()) {
                return None;
            }
        }
        None
    }

    /// Gauss–Newton on `(S − facet)·S_u = (S − facet)·S_v = 0`.
    pub(super) fn closest_point(
        &self,
        surface: &SurfaceGeometry,
        seed: (f64, f64),
    ) -> Option<(f64, f64)> {
        let f = self.facet;
        let eps = 1e-12 * (1.0 + dot(f, f).sqrt());
        let mut uv = seed;
        for _ in 0..REFINE_ITERATIONS {
            let p = point3(surface, uv);
            let r = [p[0] - f[0], p[1] - f[1], p[2] - f[2]];
            let su = to3(surface.partial_u(uv.0, uv.1));
            let sv = to3(surface.partial_v(uv.0, uv.1));
            let (a, b, c) = (dot(su, su), dot(su, sv), dot(sv, sv));
            let (g0, g1) = (dot(r, su), dot(r, sv));
            let det = a * c - b * b;
            if det.is_nan() || det.abs() <= 1e-24 * (a * c).max(1e-300) {
                return None;
            }
            let du = -(c * g0 - b * g1) / det;
            let dv = -(a * g1 - b * g0) / det;
            uv = damped_step(surface, clamp_uv(surface, uv), du, dv);
            if !(uv.0.is_finite() && uv.1.is_finite()) {
                return None;
            }
            if (du * du * a + dv * dv * c).sqrt() < eps {
                return self.near(surface, uv, self.reach).then_some(uv);
            }
        }
        None
    }
}

/// Which non-periodic coordinates of `uv` sit on (or past) a bound of
/// the natural range, and the bound value.
pub(super) fn active_bounds(surface: &SurfaceGeometry, uv: (f64, f64)) -> [Option<f64>; 2] {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let (pu, pv) = surface.parameter_periodicity();
    let bound = |x: f64, lo: f64, hi: f64, periodic: bool| {
        if periodic || !(lo.is_finite() && hi.is_finite()) {
            None
        } else if x <= lo {
            Some(lo)
        } else if x >= hi {
            Some(hi)
        } else {
            None
        }
    };
    [
        bound(uv.0, u0, u1, pu.is_some()),
        bound(uv.1, v0, v1, pv.is_some()),
    ]
}

impl PixelQuery {
    /// Box-constrained Gauss–Newton on `|S(u, v) − (eye + t·d)|²`
    /// (active set: a coordinate a step would push past its bound is
    /// pinned there and the remaining unknowns re-solved). Beside a
    /// seam the pixel's ray meets the face surface's continuation just
    /// outside the domain; the constrained minimizer is the domain
    /// boundary point nearest the ray, whose normal is the face's own
    /// boundary normal. Accepted only when the ray passes within the
    /// facet reach of it and it faces the eye.
    fn boundary_hit(&self, surface: &SurfaceGeometry, seed: (f64, f64)) -> Option<(f64, f64)> {
        let e = self.eye;
        let to_facet = [
            self.facet[0] - e[0],
            self.facet[1] - e[1],
            self.facet[2] - e[2],
        ];
        let dist = dot(to_facet, to_facet).sqrt();
        let d = unit(to_facet)?;
        let eps = 1e-12 * (1.0 + dist);
        let (mut uv, mut t) = (clamp_uv(surface, seed), dist);
        for _ in 0..REFINE_ITERATIONS {
            let p = point3(surface, uv);
            let f = [
                p[0] - e[0] - t * d[0],
                p[1] - e[1] - t * d[1],
                p[2] - e[2] - t * d[2],
            ];
            let su = to3(surface.partial_u(uv.0, uv.1));
            let sv = to3(surface.partial_v(uv.0, uv.1));
            let nd = [-d[0], -d[1], -d[2]];
            let full = solve3(su, sv, nd, [-f[0], -f[1], -f[2]])?;
            let trial = (uv.0 + full[0], uv.1 + full[1]);
            let pinned = active_bounds(surface, trial);
            let (du, dv, dt) = match pinned {
                [None, None] => (full[0], full[1], full[2]),
                [Some(ub), None] => {
                    // Pin u; least squares in (v, t).
                    let du = ub - uv.0;
                    let r = [-f[0] - su[0] * du, -f[1] - su[1] * du, -f[2] - su[2] * du];
                    let (dv, dt) = lsq2(sv, nd, r)?;
                    (du, dv, dt)
                }
                [None, Some(vb)] => {
                    let dv = vb - uv.1;
                    let r = [-f[0] - sv[0] * dv, -f[1] - sv[1] * dv, -f[2] - sv[2] * dv];
                    let (du, dt) = lsq2(su, nd, r)?;
                    (du, dv, dt)
                }
                [Some(ub), Some(vb)] => {
                    let (du, dv) = (ub - uv.0, vb - uv.1);
                    let r = [
                        -f[0] - su[0] * du - sv[0] * dv,
                        -f[1] - su[1] * du - sv[1] * dv,
                        -f[2] - su[2] * du - sv[2] * dv,
                    ];
                    (du, dv, dot(r, nd))
                }
            };
            uv = clamp_uv(surface, (uv.0 + du, uv.1 + dv));
            t += dt;
            if !(uv.0.is_finite() && uv.1.is_finite() && t.is_finite()) {
                return None;
            }
            let moved = (dot(su, su) * du * du + dot(sv, sv) * dv * dv + dt * dt).sqrt();
            if moved < eps {
                let accepted = t > 0.0
                    && pinned.iter().any(Option::is_some)
                    && self.near(surface, uv, self.reach)
                    && self.faces_eye(surface, uv, d)
                    && {
                        let q = point3(surface, uv);
                        let off = [
                            q[0] - e[0] - t * d[0],
                            q[1] - e[1] - t * d[1],
                            q[2] - e[2] - t * d[2],
                        ];
                        dot(off, off).sqrt() <= self.reach
                    };
                return accepted.then_some(uv);
            }
        }
        None
    }

    /// Bracketed ray march for grazing rays, where the ray Jacobian
    /// `[S_u S_v −d]` is near-singular and Newton wanders: the signed
    /// distance `f(t) = (P(t) − S(cp(P(t))))·n` of the ray point
    /// `P(t)` from the surface (closest-point projection, warm-started,
    /// `n` the facing normal) changes sign at the crossing. The facet's
    /// own sign says which way the crossing lies (an inscribed facet
    /// sits below the surface: toward the eye); steps of one facet
    /// reach bracket it, bisection pins it.
    fn marched_hit(&self, surface: &SurfaceGeometry, seed: (f64, f64)) -> Option<(f64, f64)> {
        let e = self.eye;
        let to_facet = [
            self.facet[0] - e[0],
            self.facet[1] - e[1],
            self.facet[2] - e[2],
        ];
        let dist = dot(to_facet, to_facet).sqrt();
        let d = unit(to_facet)?;
        let at = |t: f64| [e[0] + t * d[0], e[1] + t * d[1], e[2] + t * d[2]];
        let signed = |t: f64, warm: (f64, f64)| -> Option<(f64, (f64, f64))> {
            let p = at(t);
            let probe = PixelQuery {
                eye: e,
                facet: p,
                reach: f64::INFINITY,
                sense: self.sense,
            };
            let uv = probe.closest_point(surface, warm)?;
            let n = popkern_mesh::surface_limit_normal(surface, uv.0, uv.1)?;
            let q = point3(surface, uv);
            let r = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
            Some((self.sense * (r[0] * n.x + r[1] * n.y + r[2] * n.z), uv))
        };
        let (f0, uv0) = signed(dist, seed)?;
        // Below the surface (f < 0): the crossing is toward the eye.
        let dir = if f0 < 0.0 { -1.0 } else { 1.0 };
        let (mut a, mut fa, mut uva) = (dist, f0, uv0);
        let mut bracket = None;
        for k in 1..=(2.0 * SILHOUETTE_REACH) as usize {
            let t = dist + dir * self.reach * 0.5 * k as f64;
            if t <= 0.0 {
                break;
            }
            let (ft, uvt) = signed(t, uva)?;
            if (ft < 0.0) != (fa < 0.0) {
                bracket = Some((a, fa, t, uvt));
                break;
            }
            (a, fa, uva) = (t, ft, uvt);
        }
        let (mut lo, flo, mut hi, _) = bracket?;
        let neg_lo = flo < 0.0;
        let mut uv = uva;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            let (fm, uvm) = signed(mid, uv)?;
            uv = uvm;
            if (fm < 0.0) == neg_lo {
                lo = mid;
            } else {
                hi = mid;
            }
            if (hi - lo).abs() < 1e-13 * (1.0 + dist) {
                break;
            }
        }
        let t = 0.5 * (lo + hi);
        let p = at(t);
        let q = point3(surface, uv);
        let gap = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
        (gap <= 1e-9 * (1.0 + dist) && self.faces_eye(surface, uv, d)).then_some(uv)
    }
}

/// Locate the surface point a pixel shows: the nearest front-facing
/// intersection of its view ray with the surface near the facet,
/// over Newton solves from every seed, then a bracketed march along
/// the ray for grazing rays, then the domain-boundary point nearest
/// the ray (pixels beside a seam); failing all three, the closest surface point to the facet
/// point; failing that, the first seed.
///
/// Seeds are the facet's interpolated parameter followed by the
/// claiming triangle's vertex parameters: at a chart pole the
/// interpolated parameter inherits the pole vertex's arbitrary `u`,
/// and the ray hit lies in the direction of a ring vertex instead.
/// Convergence is relative to the facet's distance from the eye (a
/// scale-free criterion, not a model tolerance). A hit farther than
/// `reach` from the facet point belongs to another sheet of the
/// surface — a torus ray crosses the tube twice — and is rejected.
pub(super) fn locate_on_surface(
    surface: &SurfaceGeometry,
    seeds: &[(f64, f64)],
    query: &PixelQuery,
) -> ((f64, f64), SurfaceHit) {
    // The nearest front-facing root is the point the eye sees.
    if let Some((uv, _)) = seeds
        .iter()
        .filter_map(|&s| query.ray_hit(surface, s))
        .min_by(|a, b| a.1.total_cmp(&b.1))
    {
        return (uv, SurfaceHit::Ray);
    }
    if let Some(uv) = seeds.first().and_then(|&s| query.marched_hit(surface, s)) {
        return (uv, SurfaceHit::Ray);
    }
    if let Some(uv) = seeds
        .iter()
        .find_map(|&s| query.boundary_hit(surface, s))
        .filter(|&uv| !super::footprint::on_closed_seam(surface, uv, query.reach))
    {
        return (uv, SurfaceHit::DomainBoundary);
    }
    if let Some(uv) = seeds.iter().find_map(|&s| query.closest_point(surface, s)) {
        return (uv, SurfaceHit::ClosestPoint);
    }
    (
        seeds.first().copied().unwrap_or((0.0, 0.0)),
        SurfaceHit::Interpolated,
    )
}

/// Longest edge of a mesh triangle: the proximity bound for its
/// pixels' surface points.
pub(super) fn triangle_reach(vs: &[&MeshVertex; 3]) -> f64 {
    let d = |a: &MeshVertex, b: &MeshVertex| a.position.distance_to(&b.position);
    d(vs[0], vs[1]).max(d(vs[1], vs[2])).max(d(vs[2], vs[0]))
}

/// Shift `x` by whole periods to lie within half a period of `anchor`.
pub(super) fn unwrap_near(x: f64, anchor: f64, period: Option<f64>) -> f64 {
    match period {
        Some(p) if p > 0.0 => x - p * ((x - anchor) / p).round(),
        _ => x,
    }
}

/// The three vertex parameters of a triangle unwrapped onto one sheet
/// of any periodic direction, so a triangle whose vertices straddle a
/// seam in stored coordinates does not interpolate across the whole
/// period.
pub(super) fn unwrapped_vertex_uvs(
    surface: Option<&SurfaceGeometry>,
    vs: &[&MeshVertex; 3],
) -> [(f64, f64); 3] {
    let (pu, pv) = surface.map_or((None, None), |s| s.parameter_periodicity());
    let (u0, v0) = vs[0].uv;
    let w = |k: usize| {
        (
            unwrap_near(vs[k].uv.0, u0, pu),
            unwrap_near(vs[k].uv.1, v0, pv),
        )
    };
    [(u0, v0), w(1), w(2)]
}

/// Perspective-correct native `(u, v)` of a pixel from unwrapped
/// vertex parameters and perspective-correct barycentrics `c`.
pub(super) fn interpolate_uv(uvs: &[(f64, f64); 3], c: [f64; 3]) -> (f64, f64) {
    (
        c[0] * uvs[0].0 + c[1] * uvs[1].0 + c[2] * uvs[2].0,
        c[0] * uvs[0].1 + c[1] * uvs[1].1 + c[2] * uvs[2].1,
    )
}
