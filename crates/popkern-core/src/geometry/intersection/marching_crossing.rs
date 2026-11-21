//! Marching through a tangent crossing (tangency design §4.8, W3c): two
//! surfaces tangent at a point with relative curvature of both signs meet
//! in two branches crossing there, and the intersection's Jacobian is
//! singular at the point. Stepping by the normals' cross product, the
//! march slowed to the tolerance and traced residual noise round the
//! point: records of fragments a fraction of a millimetre long, and the
//! booleans' cellular results invalid (the W1.3 corpus's V cells, its
//! largest failure class). The crossing is a node of the intersection:
//! a march approaching one ends its trace with a sample at the exact
//! crossing (`contact_discovery::contact_from`, classified
//! `CrossingPoint`), and the seed loop starts each branch half from it
//! (`crossing_departures`), a sample past the zone where the surfaces lie
//! within the tolerance of each other (`√(2·tol / min|λ|)`, the relative
//! curvature's eigenvalues) solved from parameters projected there.
//! Marching on through the crossing instead let a second branch's march
//! run back over the first's approach: duplicate spans, records parallel
//! within 1e-9 (a cone and a frustum crossing at their tangency).

use super::MarchSample;
use super::marching_newton::march_project_with_iters;
use crate::geometry::intersection::contact_classify::{ContactScale, SurfaceContactKind};
use crate::geometry::intersection::contact_discovery::{SurfaceContact, contact_from};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::project_point_to_uv_near;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// One march's surfaces, windows and tolerance.
pub(super) struct Crossing<'a> {
    s1: (&'a SurfaceGeometry, Window),
    s2: (&'a SurfaceGeometry, Window),
    tolerance: &'a Tolerance,
}

/// A parameter window.
type Window = ((f64, f64), (f64, f64));

impl<'a> Crossing<'a> {
    pub(super) fn new(
        s1: (&'a SurfaceGeometry, Window),
        s2: (&'a SurfaceGeometry, Window),
        tolerance: &'a Tolerance,
    ) -> Self {
        Self { s1, s2, tolerance }
    }

    /// Push `sample` (reached from `prev`, its step's direction
    /// `fallback` when the two coincide) and, when a crossing lies just
    /// ahead, the sample at it: `true` then, and the trace ends there.
    pub(super) fn push(
        &self,
        samples: &mut Vec<MarchSample>,
        sample: MarchSample,
        prev: Point3,
        fallback: Vec3,
        step: f64,
    ) -> bool {
        let travel = (sample.point - prev).normalized().unwrap_or(fallback);
        let at = crossing_ahead(self.s1, self.s2, &sample, travel, step, self.tolerance);
        samples.push(sample);
        match at {
            Some(at) => {
                samples.push(at);
                true
            }
            None => false,
        }
    }

    /// When Newton gives out (`samples`' last sample, travelling along
    /// `fallback` with steps of `step`) at a crossing: its sample pushed,
    /// `true`, and the trace ends there. Marches of a cone and a frustum
    /// crossing at their tangency all ended `newton-exhausted` at it.
    pub(super) fn rescue(&self, samples: &mut Vec<MarchSample>, fallback: Vec3, step: f64) -> bool {
        let Some(last) = samples.last().cloned() else {
            return false;
        };
        let travel = match samples.len() {
            n if n >= 2 => (last.point - samples[n - 2].point).normalized(),
            _ => None,
        }
        .unwrap_or(fallback);
        let Some(at) = crossing_ahead(self.s1, self.s2, &last, travel, step, self.tolerance) else {
            return false;
        };
        if at.point.distance_to(&last.point) > self.tolerance.linear() {
            samples.push(at);
        } else if let Some(end) = samples.last_mut() {
            *end = at;
        }
        true
    }
}

/// The normals' sine below which the march looks for a crossing ahead.
const NEAR_TANGENT: f64 = 0.05;

/// The crossing at or near `sample` (`h` the difference step): the
/// contact, its zone (where the surfaces lie within the tolerance) and
/// the branch directions. `None` where the surfaces are not near tangent
/// or meet in no crossing there.
fn crossing_near(
    (s1, w1): (&SurfaceGeometry, Window),
    (s2, w2): (&SurfaceGeometry, Window),
    sample: &MarchSample,
    extent: f64,
    tolerance: &Tolerance,
) -> Option<(SurfaceContact, f64, Vec<Vec3>)> {
    let n1 = s1.normal_at(sample.u1, sample.v1).normalized()?;
    let n2 = s2.normal_at(sample.u2, sample.v2).normalized()?;
    if n1.cross(&n2).length() >= NEAR_TANGENT {
        return None;
    }
    let width = |w: Window| (w.0.1 - w.0.0).abs().max((w.1.1 - w.1.0).abs());
    let scale = ContactScale {
        linear: tolerance.linear(),
        extent: extent.max(tolerance.linear()),
    };
    let h = 1e-6 * width(w1).max(1e-9);
    let contact = contact_from((s1, w1), (s2, w2), (sample.u1, sample.v1), h, scale)?;
    let SurfaceContactKind::CrossingPoint { branches } = &contact.classification.kind else {
        return None;
    };
    let branches: Vec<Vec3> = branches.iter().filter_map(|b| b.normalized()).collect();
    let lambda = contact
        .classification
        .eigenvalues
        .map(|[l0, l1]| l0.abs().min(l1.abs()))
        .filter(|l| *l > 0.0)?;
    let zone = (2.0 * tolerance.linear() / lambda).sqrt();
    Some((contact, zone, branches))
}

/// The crossing a march at `sample`, travelling along `travel` (unit)
/// with steps of `step`, is about to reach (within four steps or two
/// zones ahead): the sample at it. `None` when the surfaces are not near
/// tangent there or no crossing lies ahead.
fn crossing_ahead(
    s1: (&SurfaceGeometry, Window),
    s2: (&SurfaceGeometry, Window),
    sample: &MarchSample,
    travel: Vec3,
    step: f64,
    tolerance: &Tolerance,
) -> Option<MarchSample> {
    let (contact, zone, _) = crossing_near(s1, s2, sample, 64.0 * step, tolerance)?;
    let ahead = contact.point - sample.point;
    if ahead.length() > (4.0 * step).max(2.0 * zone)
        || ahead.dot(&travel) < -10.0 * tolerance.linear()
    {
        return None;
    }
    Some(sample_at(&contact))
}

/// The march sample at a contact.
fn sample_at(contact: &SurfaceContact) -> MarchSample {
    MarchSample {
        u1: contact.uv1.0,
        v1: contact.uv1.1,
        u2: contact.uv2.0,
        v2: contact.uv2.1,
        point: contact.point,
    }
}

/// Seeds on the branch halves leaving the crossings `samples` (a trace)
/// ends at (`crossing_departures`), its mean step scaling the contact
/// search.
pub(crate) fn departures_at_ends(
    samples: &[MarchSample],
    s1: (&SurfaceGeometry, Window),
    s2: (&SurfaceGeometry, Window),
    tolerance: &Tolerance,
) -> Vec<MarchSample> {
    if samples.len() < 2 {
        return Vec::new();
    }
    let arc: f64 = samples
        .windows(2)
        .map(|w| w[0].point.distance_to(&w[1].point))
        .sum();
    let extent = 64.0 * arc / (samples.len() - 1) as f64;
    [&samples[0], &samples[samples.len() - 1]]
        .into_iter()
        .flat_map(|end| crossing_departures(s1, s2, end, extent, tolerance))
        .collect()
}

/// Where a trace ending at `end` stands at a tangent crossing (within
/// the tolerance): a seed on each branch half leaving it, twice the zone
/// out, solved from parameters projected there (Newton from the
/// crossing's own parameters starts where the Jacobian is singular).
/// Empty when `end` is no crossing. `extent` scales the contact search.
fn crossing_departures(
    s1: (&SurfaceGeometry, Window),
    s2: (&SurfaceGeometry, Window),
    end: &MarchSample,
    extent: f64,
    tolerance: &Tolerance,
) -> Vec<MarchSample> {
    let Some((contact, zone, branches)) = crossing_near(s1, s2, end, extent, tolerance) else {
        return Vec::new();
    };
    if contact.point.distance_to(&end.point) > 10.0 * tolerance.linear() {
        return Vec::new();
    }
    let at = sample_at(&contact);
    let out = 2.0 * zone;
    let mut seeds = Vec::new();
    for b in branches.iter().flat_map(|b| [*b, *b * -1.0]) {
        let target = at.point + b * out;
        let (Some(uv1), Some(uv2)) = (
            project_point_to_uv_near(s1.0, &target, at.u1, at.v1),
            project_point_to_uv_near(s2.0, &target, at.u2, at.v2),
        ) else {
            continue;
        };
        let Some((seed, _)) = march_project_with_iters(
            s1.0, uv1.0, uv1.1, s2.0, uv2.0, uv2.1, &target, &b, tolerance,
        ) else {
            continue;
        };
        // On the branch half, past the zone.
        if (seed.point - at.point).dot(&b) >= zone {
            seeds.push(seed);
        }
    }
    seeds
}
