//! Finding where two surfaces are tangent (PK-1126).
//!
//! The intersector traces where surfaces cross; where they only touch
//! it emits nothing, by contract (the tangency-silent view callers rely
//! on). This module finds the touches, beside the curves.
//!
//! A contact is a critical point of the signed gap between the
//! surfaces, read over the first surface's parameters:
//!
//! `g(u, v) = (S₁(u, v) − S₂(c)) · n₂(c)`, `c` the closest point of `S₂`,
//!
//! where `g` is within the tolerance of zero. At an isolated touch `g`
//! has an extremum there, at a crossing of branches a saddle, along a
//! tangent curve a degenerate line of critical points; its Hessian at
//! the contact is the relative curvature the classifier reads
//! (`contact_classify`). Away from contacts, along the curves where the
//! surfaces cross, `g` passes through zero with a nonzero gradient and
//! is no critical point, so traced curves are not rediscovered.
//!
//! The search is a grid over the first surface's window, seeds where
//! the grid's gap gradient is locally least and the gap small, and a
//! damped Newton on `∇g = 0` from each, its derivatives by central
//! differences of `g`. A converged point is a contact when `|g| ≤ ε`
//! and the classifier accepts it (normals parallel); the classifier
//! then types it. A minimum of the gap above `ε` is a near miss, not a
//! contact, and is reported as one for observability.

use crate::geometry::Point3;
use crate::geometry::intersection::contact_classify::{
    ContactClassification, ContactScale, classify_contact,
};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// A window of a surface's parameters: `((u0, u1), (v0, v1))`.
pub type ParamWindow = ((f64, f64), (f64, f64));

/// Where two surfaces are tangent, and what they do there.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceContact {
    /// The contact point (on the first surface).
    pub point: Point3,
    /// Its parameters on the first surface.
    pub uv1: (f64, f64),
    /// Its parameters on the second surface.
    pub uv2: (f64, f64),
    /// The signed gap there (within the tolerance).
    pub gap: f64,
    /// What the surfaces do there.
    pub classification: ContactClassification,
}

/// A local minimum of the gap above the tolerance: the surfaces come
/// close without touching (reported, never acted on).
#[derive(Debug, Clone, PartialEq)]
pub struct NearMiss {
    /// The point on the first surface.
    pub point: Point3,
    /// The gap there.
    pub gap: f64,
}

/// Everything the search found.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ContactSearch {
    /// The contacts, in the first surface's parameter order.
    pub contacts: Vec<SurfaceContact>,
    /// Near misses within `NEAR_MISS_FACTOR` tolerances.
    pub near_misses: Vec<NearMiss>,
}

/// Grid stations per direction over the first window.
const GRID: usize = 24;
/// Newton iterations per seed.
const NEWTON_STEPS: usize = 40;
/// Gaps up to this many tolerances are reported as near misses.
const NEAR_MISS_FACTOR: f64 = 100.0;

/// Find where `s1` over `w1` touches `s2` over `w2`.
pub fn find_surface_contacts(
    s1: &SurfaceGeometry,
    w1: ParamWindow,
    s2: &SurfaceGeometry,
    w2: ParamWindow,
    scale: ContactScale,
) -> ContactSearch {
    let gap = Gap { s2, w2 };
    let ((u0, u1), (v0, v1)) = w1;
    let (du, dv) = ((u1 - u0) / GRID as f64, (v1 - v0) / GRID as f64);
    let station = |i: usize, j: usize| (u0 + (i as f64 + 0.5) * du, v0 + (j as f64 + 0.5) * dv);
    // The gap on the grid.
    let mut g = vec![None; GRID * GRID];
    for i in 0..GRID {
        for j in 0..GRID {
            g[i * GRID + j] = gap.at(s1, station(i, j)).map(|(x, _)| x);
        }
    }
    // Its discrete gradient's size: central differences, wrapping
    // across a direction that spans its whole period (a cylinder's
    // seam), one-sided at a window's edge (a sphere's pole).
    let (wrap_u, wrap_v) = full_periods(s1, w1);
    let index = |i: isize, j: isize| -> Option<usize> {
        let n = GRID as isize;
        let i = if wrap_u { i.rem_euclid(n) } else { i };
        let j = if wrap_v { j.rem_euclid(n) } else { j };
        (i >= 0 && j >= 0 && i < n && j < n).then(|| i as usize * GRID + j as usize)
    };
    let point = |k: usize| {
        let (a, b) = station(k / GRID, k % GRID);
        s1.point_at(a, b)
    };
    let mut grad = vec![None; GRID * GRID];
    for i in 0..GRID as isize {
        for j in 0..GRID as isize {
            let k = i as usize * GRID + j as usize;
            let Some(c) = g[k] else {
                continue;
            };
            let p = point(k);
            let diff = |fwd: Option<usize>, back: Option<usize>| -> Option<(f64, f64)> {
                let f = fwd.and_then(|f| g[f].map(|x| (x, point(f))));
                let b = back.and_then(|b| g[b].map(|x| (x, point(b))));
                match (f, b) {
                    (Some((xf, pf)), Some((xb, pb))) => {
                        let h = 0.5 * (pf.distance_to(&pb)).max(1e-300);
                        Some(((xf - xb) / (2.0 * h), h))
                    }
                    (Some((xf, pf)), None) => {
                        let h = pf.distance_to(&p).max(1e-300);
                        Some(((xf - c) / h, h))
                    }
                    (None, Some((xb, pb))) => {
                        let h = pb.distance_to(&p).max(1e-300);
                        Some(((c - xb) / h, h))
                    }
                    (None, None) => None,
                }
            };
            let (Some((gu, hu)), Some((gv, hv))) = (
                diff(index(i + 1, j), index(i - 1, j)),
                diff(index(i, j + 1), index(i, j - 1)),
            ) else {
                continue;
            };
            grad[k] = Some(((gu * gu + gv * gv).sqrt(), c, hu.max(hv)));
        }
    }
    // Seeds: the gradient locally least among the eight neighbours, and
    // the gap small against the grid's spacing.
    let mut seeds = Vec::new();
    for i in 0..GRID {
        for j in 0..GRID {
            let Some((m, c, h)) = grad[i * GRID + j] else {
                continue;
            };
            if c.abs() > 0.5 * h + 10.0 * scale.linear {
                continue;
            }
            let least = (-1isize..=1).all(|a| {
                (-1isize..=1).all(|b| {
                    if (a, b) == (0, 0) {
                        return true;
                    }
                    index(i as isize + a, j as isize + b)
                        .and_then(|k| grad[k])
                        .is_none_or(|(o, _, _)| m <= o)
                })
            });
            if least {
                seeds.push(station(i, j));
            }
        }
    }
    let mut out = ContactSearch::default();
    let window_step = (du.abs().max(dv.abs())) * 1e-3;
    for seed in seeds {
        let Some((uv1, g0)) = newton_critical(s1, &gap, seed, window_step, w1) else {
            continue;
        };
        let point = s1.point_at(uv1.0, uv1.1);
        if g0.abs() > scale.linear {
            if g0.abs() <= NEAR_MISS_FACTOR * scale.linear
                && !out
                    .near_misses
                    .iter()
                    .any(|m| m.point.distance_to(&point) <= 10.0 * scale.linear)
            {
                out.near_misses.push(NearMiss { point, gap: g0 });
            }
            continue;
        }
        let Some((_, uv2)) = gap.at(s1, uv1) else {
            continue;
        };
        if out
            .contacts
            .iter()
            .any(|c| c.point.distance_to(&point) <= 10.0 * scale.linear)
        {
            continue;
        }
        let Ok(classification) = classify_contact(s1, uv1, s2, uv2, scale) else {
            continue;
        };
        out.contacts.push(SurfaceContact {
            point,
            uv1,
            uv2,
            gap: g0,
            classification,
        });
    }
    out.contacts.sort_by(|a, b| {
        a.uv1
            .0
            .total_cmp(&b.uv1.0)
            .then(a.uv1.1.total_cmp(&b.uv1.1))
    });
    out
}

/// Whether the window spans the whole period of `s`'s u and v.
fn full_periods(s: &SurfaceGeometry, w: ParamWindow) -> (bool, bool) {
    let (pu, pv) = s.parameter_periodicity();
    let spans = |(lo, hi): (f64, f64), p: Option<f64>| {
        p.is_some_and(|p| (hi - lo).abs() >= p * (1.0 - 1e-9))
    };
    (spans(w.0, pu), spans(w.1, pv))
}

/// The signed gap of a point of `s1` from `s2` within its window.
struct Gap<'a> {
    s2: &'a SurfaceGeometry,
    w2: ParamWindow,
}

impl Gap<'_> {
    /// The gap at `uv` of `s1` and the closest parameters on `s2`;
    /// `None` when the closest point lies outside `s2`'s window.
    fn at(&self, s1: &SurfaceGeometry, uv: (f64, f64)) -> Option<(f64, (f64, f64))> {
        let p = s1.point_at(uv.0, uv.1);
        let c = self.s2.closest_point(&p)?;
        let c = into_window(self.s2, c, self.w2)?;
        let q = self.s2.point_at(c.0, c.1);
        let n = self.s2.normal_at(c.0, c.1);
        let n = n.normalized().or_else(|| (p - q).normalized())?;
        Some(((p - q).dot(&n), c))
    }
}

/// `uv` moved by whole periods into `w`, padded by a hundredth of its
/// extent; `None` if it lies outside.
fn into_window(s: &SurfaceGeometry, uv: (f64, f64), w: ParamWindow) -> Option<(f64, f64)> {
    let (pu, pv) = s.parameter_periodicity();
    let fit = |x: f64, (lo, hi): (f64, f64), period: Option<f64>| -> Option<f64> {
        let pad = 0.01 * (hi - lo).abs();
        let (lo, hi) = (lo.min(hi) - pad, lo.max(hi) + pad);
        let mut x = x;
        if let Some(p) = period.filter(|p| *p > 0.0) {
            while x < lo {
                x += p;
            }
            while x > hi {
                x -= p;
            }
        }
        (x >= lo && x <= hi).then_some(x)
    };
    Some((fit(uv.0, w.0, pu)?, fit(uv.1, w.1, pv)?))
}

/// Damped Newton on `∇g = 0` from `seed`, derivatives by central
/// differences of step `h` (parameter units), kept inside `w`.
fn newton_critical(
    s1: &SurfaceGeometry,
    gap: &Gap<'_>,
    seed: (f64, f64),
    h: f64,
    w: ParamWindow,
) -> Option<((f64, f64), f64)> {
    let g = |u: f64, v: f64| gap.at(s1, (u, v)).map(|(x, _)| x);
    // A step leaving the window wraps across a whole period, else
    // stops at the edge (a contact at a pole lies on it).
    let (wrap_u, wrap_v) = full_periods(s1, w);
    let (pu, pv) = s1.parameter_periodicity();
    let settle = |x: f64, (lo, hi): (f64, f64), wrap: bool, period: Option<f64>| {
        let (lo, hi) = (lo.min(hi), lo.max(hi));
        match period.filter(|_| wrap) {
            Some(p) => lo + (x - lo).rem_euclid(p),
            None => x.clamp(lo, hi),
        }
    };
    let (mut u, mut v) = seed;
    let mut mu = 1e-6;
    for _ in 0..NEWTON_STEPS {
        let c = g(u, v)?;
        let (e, wv, n, s) = (g(u + h, v)?, g(u - h, v)?, g(u, v + h)?, g(u, v - h)?);
        let (ne, nw, se, sw) = (
            g(u + h, v + h)?,
            g(u - h, v + h)?,
            g(u + h, v - h)?,
            g(u - h, v - h)?,
        );
        let gu = (e - wv) / (2.0 * h);
        let gv = (n - s) / (2.0 * h);
        let huu = (e - 2.0 * c + wv) / (h * h);
        let hvv = (n - 2.0 * c + s) / (h * h);
        let huv = (ne - nw - se + sw) / (4.0 * h * h);
        // Levenberg–Marquardt on the gradient: (H + μ|H|) Δ = −∇g.
        let scale = huu.abs().max(hvv.abs()).max(1e-300);
        let (a, b, d) = (huu + mu * scale, huv, hvv + mu * scale);
        let det = a * d - b * b;
        if det.abs() <= 1e-300 {
            mu *= 10.0;
            continue;
        }
        let du = -(d * gu - b * gv) / det;
        let dv = -(-b * gu + a * gv) / det;
        let (nu, nv) = (
            settle(u + du, w.0, wrap_u, pu),
            settle(v + dv, w.1, wrap_v, pv),
        );
        // Accept a step that does not grow the gradient.
        let grad_now = (gu * gu + gv * gv).sqrt();
        let next = (
            (g(nu + h, nv)? - g(nu - h, nv)?) / (2.0 * h),
            (g(nu, nv + h)? - g(nu, nv - h)?) / (2.0 * h),
        );
        let grad_next = (next.0 * next.0 + next.1 * next.1).sqrt();
        if grad_next <= grad_now {
            (u, v) = (nu, nv);
            mu = (mu * 0.3).max(1e-12);
            if du.abs().max(dv.abs()) <= 1e-9 * h.max(1e-12) * 1e3 {
                break;
            }
        } else {
            mu *= 10.0;
            if mu > 1e8 {
                break;
            }
        }
    }
    Some(((u, v), g(u, v)?))
}

#[cfg(test)]
#[path = "contact_discovery_tests.rs"]
mod tests;
