//! Whether a boundary pcurve can meet a record in the face's chart,
//! decided before the pcurve is converted for the curve intersector
//! (PK-1168).
//!
//! The crossing searches convert every rim they test to NURBS within a
//! thousandth of the tolerance, then intersect; the intersector rejects
//! a pair whose control boxes lie farther apart than 100 tolerances
//! (PK-513), but only after the conversion. On the 64-slot disc's top,
//! 63 holes of four edges each and a circular rim, a record clipped
//! there spent 0.4 ms converting the holes' 252 edges to find none
//! crossed, and 0.2 ms converting the rim it lay well inside.
//!
//! The same rejection made first, from the pcurve itself: each form's
//! points lie in a box known in closed form (a segment's ends, an
//! ellipse's centre and larger radius, a sinusoid's sweep and
//! amplitude, a NURBS of positive weights its control points), and the
//! record, of positive weights, in its control hull. A rim whose box
//! lies farther than the pad from that hull cannot meet the record,
//! nor can a circle whose disc, shrunk by the pad, holds the hull (the
//! disc is convex). Anything else is converted and intersected as
//! before.

use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};

/// The pad, in linear tolerances: the intersector's own hull band.
pub(super) const RIM_REACH_TOLERANCES: f64 = 100.0;

/// A box in the chart: `(u_lo, v_lo, u_hi, v_hi)`.
type UvBox = (f64, f64, f64, f64);

/// The box holding every point of `pc` over its parameter range, or
/// `None` for a form with no closed-form box here.
fn pcurve_box(pc: &Pcurve) -> Option<UvBox> {
    let (t0, t1) = pc.parameter_range();
    if !(t0.is_finite() && t1.is_finite()) {
        return None;
    }
    let b = match &pc.geometry {
        PcurveGeometry::Linear { start, end } => (
            start.0.min(end.0),
            start.1.min(end.1),
            start.0.max(end.0),
            start.1.max(end.1),
        ),
        PcurveGeometry::Elliptical {
            center,
            u_radius,
            v_radius,
            ..
        } => {
            let r = u_radius.abs().max(v_radius.abs());
            (center.0 - r, center.1 - r, center.0 + r, center.1 + r)
        }
        PcurveGeometry::Sinusoidal {
            u_offset,
            v_const,
            v_cos,
            v_sin,
        } => {
            use popkern_core::geometry::PortableMath;
            let a = v_cos.phypot(*v_sin);
            (
                u_offset + t0.min(t1),
                v_const - a,
                u_offset + t0.max(t1),
                v_const + a,
            )
        }
        PcurveGeometry::Nurbs(n) => hull_box(n)?,
        PcurveGeometry::CylinderArcViaRationalQuad { .. } => return None,
    };
    [b.0, b.1, b.2, b.3]
        .iter()
        .all(|x| x.is_finite())
        .then_some(b)
}

/// The control box of a NURBS of positive weights (it lies in it).
fn hull_box(n: &NurbsCurve) -> Option<UvBox> {
    if n.control_points.is_empty() || n.weights.iter().any(|&w| w <= 0.0) {
        return None;
    }
    Some(n.control_points.iter().fold(
        (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ),
        |(a, b, c, d), p| (a.min(p.x), b.min(p.y), c.max(p.x), d.max(p.y)),
    ))
}

/// Whether the boundary pcurve `rim` provably lies farther than `pad`
/// from `query` in the chart: false when it is not proved.
pub(super) fn rim_clear_of_query(rim: &Pcurve, query: &NurbsCurve, pad: f64) -> bool {
    let Some(q) = hull_box(query) else {
        return false;
    };
    let Some(r) = pcurve_box(rim) else {
        return false;
    };
    let du = (r.0 - q.2).max(q.0 - r.2).max(0.0);
    let dv = (r.1 - q.3).max(q.1 - r.3).max(0.0);
    if du * du + dv * dv > pad * pad {
        return true;
    }
    // A circle round the record: every control point inside its disc
    // shrunk by the pad.
    if let PcurveGeometry::Elliptical {
        center,
        u_radius,
        v_radius,
        ..
    } = &rim.geometry
    {
        let radius = u_radius.abs();
        if (u_radius.abs() - v_radius.abs()).abs() <= 1e-12 * radius && radius > pad {
            let inner = radius - pad;
            return query.control_points.iter().all(|p| {
                let (du, dv) = (p.x - center.0, p.y - center.1);
                du * du + dv * dv < inner * inner
            });
        }
    }
    false
}

#[cfg(test)]
#[path = "clip_rim_reach_tests.rs"]
mod tests;
