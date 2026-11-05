//! A rounded corner of a swept path (`SweepCorner::Round`): the
//! section turned about the corner from the incoming run to the
//! outgoing one. The axis (through the corner, along `a × b`) crosses
//! the section through its centroid, so only the part of the section
//! on the corner's outer side turns: the inner part sweeps inside the
//! two runs already, and turning the whole section about a line
//! through it would make a body crossing itself. That outer part is
//! the section's boundary between its two crossings of the axis,
//! closed along the axis; it turns as a sector whose cut faces meet
//! on the axis.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;
use popkern_sweep::SweepProfile;

use crate::error::ReplayError;

fn fail(detail: impl Into<String>) -> ReplayError {
    ReplayError::FeatureFailed {
        label: "sweep".into(),
        error: crate::failure::FeatureFailure::sweep("sweep", &detail.into()).into(),
    }
}

/// The straight segment from `a` to `b` as a clamped curve of `degree`.
fn chord(a: Point3, b: Point3, degree: usize) -> NurbsCurve {
    let mut knots = vec![0.0; degree + 1];
    knots.extend(vec![1.0; degree + 1]);
    NurbsCurve {
        degree,
        knots,
        control_points: (0..=degree)
            .map(|i| a + (b - a) * (i as f64 / degree as f64))
            .collect(),
        weights: vec![1.0; degree + 1],
    }
}

/// The piece of `wire` from `t0` to `t1`.
fn piece(wire: &NurbsCurve, t0: f64, t1: f64) -> Result<NurbsCurve, ReplayError> {
    let (lo, hi) = wire.parameter_range();
    let tail = if t0 > lo {
        wire.split_at(t0).map_err(|e| fail(e.to_string()))?.1
    } else {
        wire.clone()
    };
    if t1 < hi {
        Ok(tail.split_at(t1).map_err(|e| fail(e.to_string()))?.0)
    } else {
        Ok(tail)
    }
}

/// The body the section sweeps turning from tangent `a` to `b` about
/// the corner `corner`, the section's local x along `normal` (its
/// local y along `a × normal`), its centroid on the corner.
pub(super) fn knuckle(
    profile: &SweepProfile,
    corner: Point3,
    (a, b): (Vec3, Vec3),
    normal: Vec3,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let SweepProfile::Wire(wire) = profile else {
        return Err(fail("a round corner needs an exact (analytic) profile"));
    };
    let analyzed = popkern_sweep::profile::analyze_profile(profile, &tolerance)
        .map_err(|e| fail(format!("the profile: {e}")))?;
    let centre = analyzed.centroid;
    let binormal = a.cross(&normal);
    // The inner side's direction in the section's local coordinates:
    // the outgoing tangent's part in the section plane.
    let inner = Vec3::new(b.dot(&normal), b.dot(&binormal), 0.0);
    let inner = inner
        .normalized()
        .ok_or_else(|| fail("the corner does not turn"))?;
    let side = |t: f64| (wire.point_at(t) - centre).dot(&inner);
    let (lo, hi) = wire.parameter_range();
    const SAMPLES: usize = 512;
    let mut roots = Vec::new();
    let mut prev = (lo, side(lo));
    for k in 1..=SAMPLES {
        let t = lo + (hi - lo) * k as f64 / SAMPLES as f64;
        let f = side(t);
        if (prev.1 < 0.0) != (f < 0.0) {
            let (mut x0, mut x1) = (prev.0, t);
            for _ in 0..80 {
                let m = 0.5 * (x0 + x1);
                if (side(m) < 0.0) == (prev.1 < 0.0) {
                    x0 = m;
                } else {
                    x1 = m;
                }
            }
            roots.push(0.5 * (x0 + x1));
        }
        prev = (t, f);
    }
    let [t1, t2] = roots.as_slice() else {
        return Err(fail(format!(
            "a round corner needs a section the corner's axis crosses twice; it crosses {} times",
            roots.len()
        )));
    };
    // The outer part: where `side` is negative.
    let outer = if side(0.5 * (t1 + t2)) < 0.0 {
        piece(wire, *t1, *t2)?
    } else {
        piece(wire, *t2, hi)?
            .concatenate_clamped(&piece(wire, lo, *t1)?, tolerance.linear() * 10.0)
            .map_err(|e| fail(e.to_string()))?
    };
    let (s0, s1) = outer.parameter_range();
    let (start, end) = (outer.point_at(s0), outer.point_at(s1));
    let closed = outer
        .concatenate_clamped(&chord(end, start, outer.degree), tolerance.linear() * 10.0)
        .map_err(|e| fail(e.to_string()))?;
    // Built in a canonical frame (the axis +z, the section starting in
    // the +x half of the xz plane) and placed: the same sector revolved
    // in place from another azimuth failed its union with the run it
    // meets flush (PK-1089), while the primitives' frame-placed sectors
    // serve.
    let axis = a
        .cross(&b)
        .normalized()
        .ok_or_else(|| fail("the corner does not turn"))?;
    let outward = (normal * (-inner.x) + binormal * (-inner.y))
        .normalized()
        .ok_or_else(|| fail("the corner does not turn"))?;
    let across = axis.cross(&outward);
    let local = |p: &Point3| {
        let w = normal * (p.x - centre.x) + binormal * (p.y - centre.y);
        Point3::new(w.dot(&outward), w.dot(&across), w.dot(&axis))
    };
    let section = NurbsCurve {
        control_points: closed.control_points.iter().map(local).collect(),
        ..closed
    };
    let turn = a.cross(&b).length().patan2(a.dot(&b));
    let origin = Point3::new(0.0, 0.0, 0.0);
    let built = popkern_sweep::revolve(
        &SweepProfile::Wire(section),
        &origin,
        &Vec3::new(0.0, 0.0, 1.0),
        turn,
        tolerance,
        None,
    )
    .map_err(|e| fail(format!("turning the section round a corner: {e}")))?
    .brep;
    let place =
        popkern_core::geometry::transform::Transform3::from_axis_placement(corner, axis, outward)
            .map_err(|e| fail(format!("placing the corner: {e}")))?;
    built
        .transform(&place)
        .map_err(|e| fail(format!("placing the corner: {e}")))
}
