//! A cone ring's sample azimuths about the axis, unwrapped along the
//! ring for its development (split from
//! `boundary_first_cone_development.rs` under the file-size rule).

use popkern_math::PortableMath;

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::EdgeHandle;

use super::super::BoundaryRing;

/// An edge's point at a fraction of its curve window (the ring sources'
/// convention), or `None` when the edge has no curve.
pub(in super::super) type CurveAt<'c> = &'c dyn Fn(EdgeHandle, f64) -> Option<Point3>;

/// The most curve points one ambiguous step reads.
const STEP_BUDGET: usize = 4096;

/// The axis frame: the apex, the azimuth frame (`u = 0` along `e1`) and
/// the axis.
type Frame = (Point3, Vec3, Vec3, Vec3);

/// A ring's sample azimuths about the axis from their positions,
/// unwrapped along the ring (a sample on the apex takes its
/// predecessor's), or `None` when every sample lies on the apex or the
/// ring winds round the axis other than once through the apex.
///
/// A step between two samples turning less than a quarter period, its
/// chord no longer than either end's distance from the axis, is read
/// from the positions alone (where a ruling's pcurve carries an arbitrary
/// `u` beside the apex, its positions keep one azimuth; one position
/// spread along the apex row has no chord). Any other step is ambiguous
/// from positions (a conic turning round the apex within a few samples
/// steps past half a period, or passes round it between a sample beside
/// it and one far out, its vertex a few tolerances from the apex). Along
/// one edge it is read from the edge's own curve, bisected until each
/// piece is unambiguous; else (no curve, or the curve passes the axis)
/// it takes the whole-period representative nearest the pcurve's own
/// step (along one edge; across edges, its nearest representative). The
/// ring's chart `u` is the face chart's own unwrap, which near the apex
/// may itself have read the positions' short way.
///
/// The walk starts at the first sample off the apex after one on it, so
/// a ring closing round the axis through the apex (a face spanning a
/// whole period between its seam's two uses, which meet there) takes its
/// period where it passes the apex: there the development is continuous,
/// the seam's two uses two rays from the origin.
pub(super) fn ring_azimuths(
    ring: &BoundaryRing,
    frame: Frame,
    on_apex: f64,
    curve_at: CurveAt<'_>,
) -> Option<Vec<f64>> {
    let (az, turns, through_apex) = unwrap_azimuths(ring, frame, on_apex, curve_at)?;
    (turns == 0.0 || (turns.abs() == 1.0 && through_apex)).then_some(az)
}

/// A position's azimuth about the axis, `None` on the apex (within
/// `on_apex`) or the axis.
fn azimuth((apex, e1, e2, axis): Frame, on_apex: f64, p: &Point3) -> Option<f64> {
    let d = *p - apex;
    let radial = d - axis * d.dot(&axis);
    (d.length() > on_apex && radial.length() > 0.0).then(|| radial.dot(&e2).patan2(radial.dot(&e1)))
}

/// A position's distance from the axis.
fn off_axis((apex, _, _, axis): Frame, p: &Point3) -> f64 {
    let d = *p - apex;
    (d - axis * d.dot(&axis)).length()
}

/// The azimuth step from `a` to `b` read from their positions, when it
/// is unambiguous: under a quarter period, the chord no longer than
/// either end's distance from the axis.
fn position_step(frame: Frame, on_apex: f64, (a, b): (&Point3, &Point3)) -> Option<f64> {
    let tau = std::f64::consts::TAU;
    let d = azimuth(frame, on_apex, b)? - azimuth(frame, on_apex, a)?;
    let d = d - tau * (d / tau).round();
    let beside = a.distance_to(b) <= off_axis(frame, a).min(off_axis(frame, b));
    (d.abs() < 0.25 * tau && beside).then_some(d)
}

/// The azimuth an edge's curve turns from fraction `f0` (at `p0`) to `f1`
/// (at `p1`), bisected until each piece's step is unambiguous, or `None`
/// when the curve passes the axis or the budget runs out.
fn turn_along(
    (frame, on_apex, curve_at): (Frame, f64, CurveAt<'_>),
    edge: EdgeHandle,
    (f0, f1): (f64, f64),
    (p0, p1): (Point3, Point3),
    budget: &mut usize,
) -> Option<f64> {
    if let Some(d) = position_step(frame, on_apex, (&p0, &p1)) {
        return Some(d);
    }
    let fm = 0.5 * (f0 + f1);
    if *budget == 0 || !(fm > f0.min(f1) && fm < f0.max(f1)) {
        return None;
    }
    *budget -= 1;
    let pm = curve_at(edge, fm)?;
    azimuth(frame, on_apex, &pm)?;
    let ctx = (frame, on_apex, curve_at);
    Some(
        turn_along(ctx, edge, (f0, fm), (p0, pm), budget)?
            + turn_along(ctx, edge, (fm, f1), (pm, p1), budget)?,
    )
}

/// A ring's unwrapped sample azimuths (as [`ring_azimuths`]), the whole
/// turns it makes round the axis and whether it passes the apex.
pub(super) fn unwrap_azimuths(
    ring: &BoundaryRing,
    frame: Frame,
    on_apex: f64,
    curve_at: CurveAt<'_>,
) -> Option<(Vec<f64>, f64, bool)> {
    let tau = std::f64::consts::TAU;
    let n = ring.pos.len();
    let raw = |p: &Point3| azimuth(frame, on_apex, p);
    let nearest = |x: f64, to: f64| x + tau * ((to - x) / tau).round();
    let regular: Vec<bool> = ring.pos.iter().map(|p| raw(p).is_some()).collect();
    let after_apex = (0..n).find(|&i| regular[i] && !regular[(i + n - 1) % n]);
    let first = after_apex.or_else(|| (0..n).find(|&i| regular[i]))?;
    let mut az = vec![0.0; n];
    let start = raw(&ring.pos[first])?;
    let (mut prev, mut prev_i) = (start, first);
    for step in 1..=n {
        let i = (first + step) % n;
        if let Some(a) = raw(&ring.pos[i]) {
            let (p0, p1) = (ring.pos[prev_i], ring.pos[i]);
            let along = || {
                // One edge's segment: its own curve.
                let (edge, f0, f1) = ring.sources[prev_i];
                (prev_i == (i + n - 1) % n).then_some(())?;
                let mut budget = STEP_BUDGET;
                let ctx = (frame, on_apex, curve_at);
                turn_along(ctx, edge, (f0, f1), (p0, p1), &mut budget)
            };
            let d = position_step(frame, on_apex, (&p0, &p1))
                .or_else(along)
                .unwrap_or_else(|| {
                    // Along one edge the pcurve's u is continuous; across
                    // edges their branches may differ by whole periods.
                    let d = nearest(a - prev, 0.0);
                    let du = ring.uv[i].0 - ring.uv[prev_i].0;
                    let one_edge = ring.sources[prev_i].0 == ring.sources[(i + n - 1) % n].0;
                    nearest(d, if one_edge { du } else { nearest(du, 0.0) })
                });
            prev += d;
            prev_i = i;
        }
        if step < n {
            az[i] = prev;
        }
    }
    az[first] = start;
    // Back round to the first sample: the whole periods the ring turned.
    let turns = ((prev - start) / tau).round();
    Some((az, turns, after_apex.is_some()))
}
