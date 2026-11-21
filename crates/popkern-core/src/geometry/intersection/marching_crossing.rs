//! Marching through a tangent crossing (tangency design §4.8, W3c): two
//! surfaces tangent at a point with relative curvature of both signs meet
//! in two branches crossing there, and the intersection's Jacobian is
//! singular at the point. Stepping by the normals' cross product, the
//! march slowed to the tolerance and traced residual noise round the
//! point: records of fragments a fraction of a millimetre long, and the
//! booleans' cellular results invalid (the W1.3 corpus's V cells, its
//! largest failure class). Approaching one, the march places a sample at
//! the exact crossing (`contact_discovery::contact_from`, classified
//! `CrossingPoint`) and leaves it along the branch running on from the
//! way it came, past the zone where the surfaces lie within the tolerance
//! of each other (`√(2·tol / min|λ|)`, the relative curvature's
//! eigenvalues), where the Jacobian is regular again.

use super::MarchSample;
use crate::geometry::intersection::contact_classify::{ContactScale, SurfaceContactKind};
use crate::geometry::intersection::contact_discovery::contact_from;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// One march's surfaces, windows and tolerance, and the crossings it has
/// passed.
pub(super) struct Crossing<'a> {
    s1: (&'a SurfaceGeometry, Window),
    s2: (&'a SurfaceGeometry, Window),
    tolerance: &'a Tolerance,
    crossed: Vec<Point3>,
}

/// A parameter window.
type Window = ((f64, f64), (f64, f64));

impl<'a> Crossing<'a> {
    pub(super) fn new(
        s1: (&'a SurfaceGeometry, Window),
        s2: (&'a SurfaceGeometry, Window),
        tolerance: &'a Tolerance,
    ) -> Self {
        Self {
            s1,
            s2,
            tolerance,
            crossed: Vec::new(),
        }
    }

    /// Push `sample` (reached from `prev`, its step's direction
    /// `fallback` when the two coincide) and, when a crossing lies just
    /// ahead, the sample at it: then the direction and step to leave it
    /// by.
    pub(super) fn push(
        &mut self,
        samples: &mut Vec<MarchSample>,
        sample: MarchSample,
        prev: Point3,
        fallback: Vec3,
        step: f64,
    ) -> Option<(Vec3, f64)> {
        let travel = (sample.point - prev).normalized().unwrap_or(fallback);
        let next = crossing_ahead(
            self.s1,
            self.s2,
            &sample,
            travel,
            step,
            self.tolerance,
            &self.crossed,
        );
        samples.push(sample);
        let (at, leave, out) = next?;
        self.crossed.push(at.point);
        samples.push(at);
        Some((leave, out))
    }
}

/// The normals' sine below which the march looks for a crossing ahead.
const NEAR_TANGENT: f64 = 0.05;

/// The crossing a march at `sample`, travelling along `travel` (unit)
/// with steps of `step`, is about to pass through: the sample at it, the
/// direction to leave it in (the branch running on from `travel`) and the
/// step to leave by. `None` when the surfaces are not near tangent there,
/// no crossing lies within a few steps ahead, or it was crossed already
/// (`crossed`).
#[allow(clippy::too_many_arguments)]
fn crossing_ahead(
    (s1, w1): (&SurfaceGeometry, Window),
    (s2, w2): (&SurfaceGeometry, Window),
    sample: &MarchSample,
    travel: Vec3,
    step: f64,
    tolerance: &Tolerance,
    crossed: &[Point3],
) -> Option<(MarchSample, Vec3, f64)> {
    let n1 = s1.normal_at(sample.u1, sample.v1).normalized()?;
    let n2 = s2.normal_at(sample.u2, sample.v2).normalized()?;
    if n1.cross(&n2).length() >= NEAR_TANGENT {
        return None;
    }
    let extent = |w: Window| (w.0.1 - w.0.0).abs().max((w.1.1 - w.1.0).abs());
    let scale = ContactScale {
        linear: tolerance.linear(),
        extent: (64.0 * step).max(tolerance.linear()),
    };
    let h = 1e-6 * extent(w1).max(1e-9);
    let contact = contact_from((s1, w1), (s2, w2), (sample.u1, sample.v1), h, scale)?;
    let SurfaceContactKind::CrossingPoint { branches } = &contact.classification.kind else {
        return None;
    };
    let p = contact.point;
    let ahead = p - sample.point;
    if ahead.length() > 4.0 * step
        || ahead.dot(&travel) < 0.0
        || crossed
            .iter()
            .any(|q| q.distance_to(&p) <= 10.0 * tolerance.linear())
    {
        return None;
    }
    // The branch running on from the way the march came.
    let leave = branches
        .iter()
        .filter_map(|b| b.normalized())
        .map(|b| if b.dot(&travel) < 0.0 { b * -1.0 } else { b })
        .max_by(|a, b| a.dot(&travel).total_cmp(&b.dot(&travel)))?;
    let lambda = contact
        .classification
        .eigenvalues
        .map(|[l0, l1]| l0.abs().min(l1.abs()))
        .filter(|l| *l > 0.0)?;
    let zone = (2.0 * tolerance.linear() / lambda).sqrt();
    let at = MarchSample {
        u1: contact.uv1.0,
        v1: contact.uv1.1,
        u2: contact.uv2.0,
        v2: contact.uv2.1,
        point: p,
    };
    Some((at, leave, step.max(2.0 * zone)))
}
