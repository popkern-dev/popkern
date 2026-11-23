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

/// Find where `s1` over `w1` touches `s2` over `w2`: the same contacts
/// in either order. The search is seeded from both surfaces (a grid
/// over one surface's window can step over a contact the other's finds:
/// a plane's grid missed a paraboloid's tip resting on it, PK-1179),
/// each contact expressed on `s1` and classified in this order, one per
/// place.
pub fn find_surface_contacts(
    s1: &SurfaceGeometry,
    w1: ParamWindow,
    s2: &SurfaceGeometry,
    w2: ParamWindow,
    scale: ContactScale,
) -> ContactSearch {
    let mut out = one_sided(s1, w1, s2, w2, scale);
    let back = one_sided(s2, w2, s1, w1, scale);
    let near = 10.0 * scale.linear;
    for c in back.contacts {
        // Its foot on s1 is the contact there, read off the search's own
        // parameters (a closest point onto s2 fails at its pole).
        let (uv1, uv2) = (c.uv2, c.uv1);
        let point = s1.point_at(uv1.0, uv1.1);
        // Along s2's normal, or s1's where s2's chart has none (a
        // revolved cup's tip): parallel at a contact.
        let Some(n2) = s2
            .normal_at(uv2.0, uv2.1)
            .normalized()
            .or_else(|| s1.normal_at(uv1.0, uv1.1).normalized())
        else {
            continue;
        };
        let g0 = (point - s2.point_at(uv2.0, uv2.1)).dot(&n2);
        if g0.abs() > scale.linear {
            continue;
        }
        let gap = Gap { s2, w2 };
        let same = |o: &SurfaceContact| {
            o.point.distance_to(&point) <= near || one_zone(s1, &gap, o.uv1, uv1, scale.linear)
        };
        let Some(slot) = slot_for(&out.contacts, g0, same) else {
            continue;
        };
        let Ok(classification) = classify_contact(s1, uv1, s2, uv2, scale) else {
            continue;
        };
        admit(
            &mut out.contacts,
            slot,
            SurfaceContact {
                point,
                uv1,
                uv2,
                gap: g0,
                classification,
            },
        );
    }
    for m in back.near_misses {
        // On s1, at the foot of the point found on s2.
        let Some(point) = s1
            .closest_point(&m.point)
            .and_then(|c| into_window(s1, c, w1))
            .map(|c| s1.point_at(c.0, c.1))
        else {
            continue;
        };
        if !out
            .near_misses
            .iter()
            .any(|o| o.point.distance_to(&point) <= near)
        {
            out.near_misses.push(NearMiss { point, gap: m.gap });
        }
    }
    out.contacts.sort_by(|a, b| {
        a.uv1
            .0
            .total_cmp(&b.uv1.0)
            .then(a.uv1.1.total_cmp(&b.uv1.1))
    });
    out
}

/// [`find_surface_contacts`] seeded from `s1` alone.
fn one_sided(
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
        let same = |c: &SurfaceContact| {
            c.point.distance_to(&point) <= 10.0 * scale.linear
                || one_zone(s1, &gap, c.uv1, uv1, scale.linear)
        };
        let Some(slot) = slot_for(&out.contacts, g0, same) else {
            continue;
        };
        let Ok(classification) = classify_contact(s1, uv1, s2, uv2, scale) else {
            continue;
        };
        admit(
            &mut out.contacts,
            slot,
            SurfaceContact {
                point,
                uv1,
                uv2,
                gap: g0,
                classification,
            },
        );
    }
    out.contacts.sort_by(|a, b| {
        a.uv1
            .0
            .total_cmp(&b.uv1.0)
            .then(a.uv1.1.total_cmp(&b.uv1.1))
    });
    out
}

/// The contact Newton reaches on the gap's gradient from `seed` on `s1`
/// (derivatives by central differences of `step`, parameter units, kept
/// inside `w1`), classified, when its gap is within `scale.linear`: a
/// known meeting place refined to the tangency it stands for (two
/// records of one face pair crossing meet only at a tangent crossing,
/// tangency design §4.8). `None` when Newton leaves the windows or the
/// gap there is past the tolerance.
pub fn contact_from(
    (s1, w1): (&SurfaceGeometry, ParamWindow),
    (s2, w2): (&SurfaceGeometry, ParamWindow),
    seed: (f64, f64),
    step: f64,
    scale: ContactScale,
) -> Option<SurfaceContact> {
    let gap = Gap { s2, w2 };
    let (uv1, g0) = newton_critical(s1, &gap, seed, step, w1)?;
    if g0.abs() > scale.linear {
        return None;
    }
    let (_, uv2) = gap.at(s1, uv1)?;
    let classification = classify_contact(s1, uv1, s2, uv2, scale).ok()?;
    Some(SurfaceContact {
        point: s1.point_at(uv1.0, uv1.1),
        uv1,
        uv2,
        gap: g0,
        classification,
    })
}

/// Whether two contacts found at `a` and `b` on `s1` are one: the gap
/// stays within the tolerance all the way between them (the linear
/// tolerance is the resolution, PK-1131). A flat touch (a paraboloid's
/// tip in a ball of its curvature, the gap quartic) holds its gap
/// within the tolerance far past ten tolerances, and Newton from
/// several seeds stopped at points 2e-5 apart there, each kept as its
/// own touch. Between them along the chart, the shorter way round a
/// period.
fn one_zone(
    s1: &SurfaceGeometry,
    gap: &Gap<'_>,
    a: (f64, f64),
    b: (f64, f64),
    linear: f64,
) -> bool {
    let (pu, pv) = crate::geometry::chart::ChartDescriptor::evaluation_periods_of(s1);
    let toward = |x: f64, y: f64, period: Option<f64>| {
        let d = y - x;
        match period.filter(|p| *p > 0.0) {
            Some(p) => d - p * (d / p).round(),
            None => d,
        }
    };
    let (du, dv) = (toward(a.0, b.0, pu), toward(a.1, b.1, pv));
    (1..ZONE_SAMPLES).all(|k| {
        let t = k as f64 / ZONE_SAMPLES as f64;
        gap.at(s1, (a.0 + t * du, a.1 + t * dv))
            .is_some_and(|(x, _)| x.abs() <= linear)
    })
}

/// Where a contact of gap `g0` goes among `found`: `None` when one it is
/// the same as (`same`) already stands at a gap as small, else
/// `Some(Some(i))` to replace the `i`th, `Some(None)` to add it. Within
/// one zone the touch is where the gap is extreme (a flat touch's tip,
/// not a point beside it where Newton stopped).
fn slot_for(
    found: &[SurfaceContact],
    g0: f64,
    same: impl Fn(&SurfaceContact) -> bool,
) -> Option<Option<usize>> {
    match found.iter().position(same) {
        Some(i) if found[i].gap.abs() <= g0.abs() => None,
        Some(i) => Some(Some(i)),
        None => Some(None),
    }
}

/// Puts `c` in the slot [`slot_for`] chose.
fn admit(found: &mut Vec<SurfaceContact>, slot: Option<usize>, c: SurfaceContact) {
    match slot {
        Some(i) => found[i] = c,
        None => found.push(c),
    }
}

/// Intervals of the segment [`one_zone`] samples.
const ZONE_SAMPLES: usize = 8;

/// Whether the window spans the whole period of `s`'s u and v.
fn full_periods(s: &SurfaceGeometry, w: ParamWindow) -> (bool, bool) {
    let (pu, pv) = crate::geometry::chart::ChartDescriptor::evaluation_periods_of(s);
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
    let (pu, pv) = crate::geometry::chart::ChartDescriptor::evaluation_periods_of(s);
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

#[path = "contact_discovery_newton.rs"]
mod newton;
use newton::newton_critical;

#[cfg(test)]
#[path = "contact_discovery_twin_tests.rs"]
mod twin_tests;

#[cfg(test)]
#[path = "contact_discovery_tests.rs"]
mod tests;
