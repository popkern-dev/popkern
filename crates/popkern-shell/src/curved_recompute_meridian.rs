//! A sphere face's seam recomputed on its offset sphere (PK-1205).
//!
//! The seam of a revolved sphere cap runs along a meridian from the
//! opening rim to the far pole. Offset, the rim vertex stays on the
//! opening's plane while the pole moves along the axis, so the new seam
//! is not the old one scaled about the centre (the parallel map then
//! declined: a cap opened at z = 3 put its inner rim vertex at
//! (2.646, 0, 3), where the scaled seam starts at (3.2, 0, 2.4)), and
//! the stale seam left the inner skin's loop open by the thickness. It
//! is the arc of the offset sphere's great circle through both new ends
//! in the old seam's plane: built exact, as a rational quadratic in
//! spans of at most a quarter turn.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_math::PortableMath;

/// The meridian arc of the offset sphere `surf` from `start_pos` to
/// `end_pos`, in the plane of the old seam `old`; `None` when the surface
/// is no sphere, the ends are off it, or the old seam is not a meridian
/// arc of the same centre.
pub(crate) fn sphere_meridian_seam_curve(
    old: &CurveOnEdge,
    surf: &SurfaceGeometry,
    start_pos: Point3,
    end_pos: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    let SurfaceGeometry::Sphere { center, radius, .. } = surf else {
        return None;
    };
    let lin = tolerance.linear();
    // A seam runs between two places; a closed edge (a free rim round
    // the sphere) is no seam, and the parallel map serves it.
    if start_pos.distance_to(&end_pos) <= lin {
        return None;
    }
    let on = |p: &Point3| ((*p - *center).length() - radius).abs() <= lin;
    if !on(&start_pos) || !on(&end_pos) {
        return None;
    }
    // The old seam's plane through the centre: its normal from two
    // spread samples of the old curve.
    let (t0, t1) = (old.t_start, old.t_end);
    let a = old.geometry.point_at(t0) - *center;
    let b = old.geometry.point_at(t0 + 0.61 * (t1 - t0)) - *center;
    let n = a.cross(&b).normalized()?;
    let (u, v) = (start_pos - *center, end_pos - *center);
    if u.dot(&n).abs() > lin || v.dot(&n).abs() > lin {
        return None;
    }
    // Sweep from u to v about n, the way round the old seam ran: of the
    // two arcs, the one whose middle lies toward the old seam's middle.
    let e1 = u.normalized()?;
    let e2 = n.cross(&e1);
    let ccw = {
        let s = v.dot(&e2).patan2(v.dot(&e1));
        if s < 0.0 {
            s + std::f64::consts::TAU
        } else {
            s
        }
    };
    let old_mid = old.geometry.point_at(0.5 * (t0 + t1)) - *center;
    let toward = |sweep: f64| {
        let h = 0.5 * sweep;
        (e1 * h.pcos() + e2 * h.psin()).dot(&old_mid)
    };
    let cw = ccw - std::f64::consts::TAU;
    let sweep = if toward(ccw) >= toward(cw) { ccw } else { cw };
    if sweep.abs() <= lin / radius {
        return None;
    }
    let spans = (sweep.abs() / std::f64::consts::FRAC_PI_2).ceil().max(1.0) as usize;
    let step = sweep / spans as f64;
    let w = (0.5 * step).pcos();
    let at = |angle: f64| *center + (e1 * angle.pcos() + e2 * angle.psin()) * *radius;
    let mut control_points = vec![start_pos];
    let mut weights = vec![1.0];
    let mut knots = vec![0.0; 3];
    for k in 0..spans {
        let mid = (k as f64 + 0.5) * step;
        let shoulder = *center + (e1 * mid.pcos() + e2 * mid.psin()) * (*radius / w);
        control_points.push(shoulder);
        weights.push(w);
        control_points.push(if k + 1 == spans {
            end_pos
        } else {
            at((k + 1) as f64 * step)
        });
        weights.push(1.0);
        let knot = (k + 1) as f64 / spans as f64;
        knots.extend([knot, knot]);
    }
    knots.pop();
    knots.push(1.0);
    knots.push(1.0);
    let curve = NurbsCurve {
        degree: 2,
        knots,
        control_points,
        weights,
    };
    Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve),
        t_start: 0.0,
        t_end: 1.0,
    })
}

/// The circle a plane through a torus's axis cuts from the offset torus
/// (PK-1218): such a plane meets a torus in its meridian circles, each
/// about a point of the core circle with the minor radius. The coaxial
/// arm takes planes square to the axis only, so a half torus shelled
/// open at its cut ends kept its inner rims at the outer minor radius
/// and the wall opened by the thickness. The old rim (a meridian circle
/// of the original torus) keeps its centre, which the offset leaves on
/// the core circle, and its angles in its own plane; its radius becomes
/// the offset's minor radius.
pub(crate) fn torus_meridian_edge_curve(
    old: &CurveOnEdge,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    start_pos: Point3,
    end_pos: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    let lin = tolerance.linear();
    let (torus, plane) = match (surf1, surf2) {
        (t @ SurfaceGeometry::Torus { .. }, p @ SurfaceGeometry::Plane { .. })
        | (p @ SurfaceGeometry::Plane { .. }, t @ SurfaceGeometry::Torus { .. }) => (t, p),
        _ => return None,
    };
    let (
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        },
        SurfaceGeometry::Plane { origin, normal, .. },
    ) = (torus, plane)
    else {
        return None;
    };
    let (a, n) = (axis.normalized()?, normal.normalized()?);
    // The plane holds the axis: square to it, through the centre.
    if a.dot(&n).abs() > 1e-9 || (*center - *origin).dot(&n).abs() > lin {
        return None;
    }
    let CurveGeometry::Arc {
        center: c0,
        normal: arc_normal,
        start_angle,
        end_angle,
        ..
    } = old.geometry
    else {
        return None;
    };
    // The old centre on the core circle, in the plane.
    let rel = c0 - *center;
    let radial = rel - a * rel.dot(&a);
    if (radial.length() - major_radius).abs() > lin
        || rel.dot(&a).abs() > lin
        || arc_normal.cross(&n).length() > 1e-9
    {
        return None;
    }
    let candidate = CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center: c0,
            normal: arc_normal,
            radius: *minor_radius,
            start_angle,
            end_angle,
        },
        t_start: old.t_start,
        t_end: old.t_end,
    };
    let ends_meet = candidate
        .geometry
        .point_at(candidate.t_start)
        .distance_to(&start_pos)
        <= lin
        && candidate
            .geometry
            .point_at(candidate.t_end)
            .distance_to(&end_pos)
            <= lin;
    ends_meet.then_some(candidate)
}
