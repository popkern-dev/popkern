//! Draft edges of oblique sections: where a drafted plane cuts a
//! cylinder (an exact ellipse) or touches a cone (its ruling). The
//! ellipse: A hole's rim in a face drafted about the floor (a box drilled
//! through its ends, the ends drafted 8 degrees) had no recompute path,
//! kept its circle and left the face's loop open by the lean.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_math::PortableMath;

use super::edge::EdgeCurveOutcome;

/// The section of the cylinder `(origin, axis, radius)` by the plane
/// `(point, normal)`, run as the old curve ran from `start` (to `end`;
/// round the whole ellipse when `closed`). Its semi-axes: across the
/// plane's tilt (radius, square to the axis and the normal) and along
/// it (radius over the cosine of the tilt).
pub(super) fn plane_cylinder_ellipse(
    (origin, axis, radius): (Point3, Vec3, f64),
    (point, normal): (Point3, Vec3),
    old: Option<&CurveOnEdge>,
    start: Point3,
    end: Point3,
    closed: bool,
    tolerance: &Tolerance,
) -> EdgeCurveOutcome {
    let (Some(a), Some(n)) = (axis.normalized(), normal.normalized()) else {
        return EdgeCurveOutcome::Degenerate;
    };
    let an = a.dot(&n);
    if an.abs() < 1.0e-9 {
        // A plane along the axis cuts the cylinder in lines.
        return EdgeCurveOutcome::UnsupportedPair;
    }
    let centre = origin + a * ((point - origin).dot(&n) / an);
    let tilt = n - a * an;
    let e1 = tilt.normalized().unwrap_or_else(|| {
        let seed = if a.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        (seed - a * seed.dot(&a)).normalized().unwrap_or(seed)
    });
    let e2 = a.cross(&e1);
    let u = (e1 - a * (e1.dot(&n) / an)) * radius;
    let mut v = e2 * radius;
    // Keep the old curve's way round.
    if let Some(old) = old {
        let tm = 0.5 * (old.t_start + old.t_end);
        let p = old.geometry.point_at(tm);
        let tangent = old.geometry.derivative_at(tm, 1);
        let d = p - centre;
        let theta = (d.dot(&v) / v.dot(&v)).patan2(d.dot(&u) / u.dot(&u));
        let along = v * theta.pcos() - u * theta.psin();
        if along.dot(&tangent) < 0.0 {
            v = v * -1.0;
        }
    }
    let angle_of = |p: Point3| {
        let d = p - centre;
        (d.dot(&v) / v.dot(&v)).patan2(d.dot(&u) / u.dot(&u))
    };
    let on = |p: Point3, t: f64| {
        (centre + u * t.pcos() + v * t.psin()).distance_to(&p) <= 10.0 * tolerance.linear()
    };
    let a0 = angle_of(start);
    let a1 = if closed {
        a0 + std::f64::consts::TAU
    } else {
        let mut t = angle_of(end);
        while t <= a0 {
            t += std::f64::consts::TAU;
        }
        t
    };
    if !on(start, a0) || !on(end, a1) {
        return EdgeCurveOutcome::UnsupportedPair;
    }
    EdgeCurveOutcome::Recomputed(CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center: centre,
            u_basis: u,
            v_basis: v,
            start_angle: a0,
            end_angle: a1,
        },
        t_start: a0,
        t_end: a1,
    })
}

/// A plane through a cone's apex and touching it (a drafted side beside
/// the quarter cone a drafted round becomes) meets it in one ruling: the
/// line from the apex along the axis's projection onto the plane. The
/// edge between them was the line where the plane met the cylinder; left
/// stale, it kept its top end h·tan(angle) off both drafted faces
/// (a rounded box drafted 5 degrees opened by 0.8749 at every round).
/// Any other oblique plane meets the cone in a conic and keeps the
/// unsupported-pair fallback.
#[allow(clippy::too_many_arguments)]
pub(super) fn tangent_ruling(
    apex: popkern_core::geometry::Point3,
    axis: Vec3,
    half_angle: f64,
    plane_origin: popkern_core::geometry::Point3,
    n: Vec3,
    start_pos: popkern_core::geometry::Point3,
    end_pos: popkern_core::geometry::Point3,
    tolerance: &Tolerance,
) -> EdgeCurveOutcome {
    let Some(a) = axis.normalized() else {
        return EdgeCurveOutcome::Degenerate;
    };
    let scale = 1.0 + (start_pos - apex).length();
    if (apex - plane_origin).dot(&n).abs() > tolerance.linear() * scale {
        return EdgeCurveOutcome::UnsupportedPair;
    }
    // Touching, not cutting: the plane leans from the axis by the
    // cone's half angle (its normal by the complement, |a·n| = sin α).
    if (a.dot(&n).abs() - half_angle.psin()).abs() > 1.0e-7 {
        return EdgeCurveOutcome::UnsupportedPair;
    }
    let Some(direction) = (a - n * a.dot(&n)).normalized() else {
        return EdgeCurveOutcome::Degenerate;
    };
    EdgeCurveOutcome::Recomputed(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: apex,
            direction,
        },
        t_start: (start_pos - apex).dot(&direction),
        t_end: (end_pos - apex).dot(&direction),
    })
}
