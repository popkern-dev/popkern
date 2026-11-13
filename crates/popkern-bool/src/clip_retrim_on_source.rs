//! A clipped end refined onto a junction lying on its source curve is a
//! re-trim of the source, keeping the piece's parameter domain
//! (PK-1182).
//!
//! The endpoint refinement moves a clipped record piece's end onto its
//! exact junction. At a shallow crossing that move runs along the
//! record: dragging the end control point there bent two balls' overlap
//! circle 5.5e-6 off both spheres (a 4.75e-5 move), and re-authoring
//! around a pinned end rings the pin into the next segments by design.
//! A junction on the source is reached by taking the source itself over
//! the corrected parameter span, re-authored onto the piece's own
//! domain (downstream window ranges and node parameters read it), the
//! three curves alike (they share the source's parameter).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval as _, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;

use crate::clip_endpoint_far::project_onto_curve_dense;
use crate::clip_reauthor::reauthor_pinned_run;

/// Whether moving `curve`'s end control point (`end`: its last) onto `r`
/// takes the curve farther than `limit` from itself anywhere.
pub(crate) fn drag_bends(curve: &NurbsCurve, end: bool, r: Point3, limit: f64) -> bool {
    let mut dragged = curve.clone();
    let slot = if end {
        dragged.control_points.last_mut()
    } else {
        dragged.control_points.first_mut()
    };
    if let Some(cp) = slot {
        *cp = r;
    }
    let (t0, t1) = dragged.parameter_range();
    (1..32).any(|k| {
        let q = dragged.point_at(t0 + (t1 - t0) * f64::from(k) / 32.0);
        project_onto_curve_dense(curve, &q)
            .is_none_or(|t| curve.point_at(t).distance_to(&q) > limit)
    })
}

/// `piece` re-made as `source` from the other end's parameter to
/// `t_star` (the junction `r`'s), on the piece's domain; `false` (piece
/// untouched) when the other end or the piece's middle does not lie on
/// that span of the source, or a fit fails.
pub(crate) fn retrim_on_source(
    piece: &mut CurveOnSurface,
    source: &CurveOnSurface,
    (end, t_star, r): (bool, f64, Point3),
    (surf_a, surf_b): (&SurfaceGeometry, &SurfaceGeometry),
    tolerance: &Tolerance,
) -> bool {
    let linear = tolerance.linear();
    let (p0, p1) = piece.curve_3d.parameter_range();
    let on_source = |q: Point3| {
        project_onto_curve_dense(&source.curve_3d, &q)
            .filter(|&t| source.curve_3d.point_at(t).distance_to(&q) <= linear)
    };
    let other = piece.curve_3d.point_at(if end { p0 } else { p1 });
    let (Some(t_other), Some(t_mid)) = (
        on_source(other),
        on_source(piece.curve_3d.point_at(0.5 * (p0 + p1))),
    ) else {
        return false;
    };
    let (mut s0, mut s1) = if end {
        (t_other, t_star)
    } else {
        (t_star, t_other)
    };
    // The span must hold the piece; a closed source's span may run round
    // its closure (a piece from 0.752 through 0 to 0.748), taken a
    // period on.
    let (a, b) = source.curve_3d.parameter_range();
    let period = b - a;
    let holds = |s0: f64, s1: f64, t: f64| (s0.min(s1)..=s0.max(s1)).contains(&t);
    if (p1 - p0).abs() < 1e-300 {
        return false;
    }
    if !holds(s0, s1, t_mid) {
        let closed = source.is_closed
            || source
                .curve_3d
                .point_at(a)
                .distance_to(&source.curve_3d.point_at(b))
                <= linear;
        if !closed {
            return false;
        }
        if s0 < s1 {
            s0 += period;
        } else {
            s1 += period;
        }
        if !holds(s0, s1, t_mid) && !holds(s0, s1, t_mid + period) {
            return false;
        }
    }
    let map = move |u: f64| s0 + (s1 - s0) * (u - p0) / (p1 - p0);
    // A source curve at `t`, continued past its closure by its own
    // end-to-start offset (a whole period of a chart curve, none in
    // space), so a run across the closure stays continuous.
    let at = move |c: &NurbsCurve, t: f64| {
        let jump = c.point_at(b) - c.point_at(a);
        if t > b {
            c.point_at(t - period) + jump
        } else if t < a {
            c.point_at(t + period) - jump
        } else {
            c.point_at(t)
        }
    };
    let budget = 0.1 * linear;
    let fit = |src: &NurbsCurve, pin_to: Point3, to_space: &dyn Fn(&Point3) -> Point3| {
        reauthor_pinned_run(
            (p0, p1),
            &|ts: &[f64]| ts.iter().map(|&u| at(src, map(u))).collect(),
            &|pts: &mut [Point3]| {
                let slot = if end { pts.last_mut() } else { pts.first_mut() };
                if let Some(q) = slot {
                    *q = pin_to;
                }
            },
            to_space,
            budget,
        )
    };
    let t_end = if end { s1 } else { s0 };
    let uv_at = |c: &NurbsCurve| at(c, t_end);
    let (Some(c3), Some(u1), Some(u2)) = (
        fit(&source.curve_3d, r, &|p: &Point3| *p),
        fit(
            &source.curve_uv1,
            uv_at(&source.curve_uv1),
            &|q: &Point3| surf_a.point_at(q.x, q.y),
        ),
        fit(
            &source.curve_uv2,
            uv_at(&source.curve_uv2),
            &|q: &Point3| surf_b.point_at(q.x, q.y),
        ),
    ) else {
        return false;
    };
    (piece.curve_3d, piece.curve_uv1, piece.curve_uv2) = (c3, u1, u2);
    true
}
