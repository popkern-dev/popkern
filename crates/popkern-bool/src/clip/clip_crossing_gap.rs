//! Whether a reported uv crossing is one in space (PK-1180), split
//! from the crossing finder under the 500-line file cap.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::HalfEdgeHandle;

use super::clip_uv::UvCrossing;

/// Whether a crossing reported at `(t_boundary, t_curve)` lies past the
/// linear reach in space by [`crossing_gap`] (`false` with no surface
/// to read it on).
pub(super) fn afar(
    surface: Option<&SurfaceGeometry>,
    curves: (&NurbsCurve, &NurbsCurve),
    params: (f64, f64),
    shift: f64,
    tolerance: &Tolerance,
) -> bool {
    surface.is_some_and(|s| {
        let gap = crossing_gap(s, curves, params, shift);
        gap.is_finite() && !tolerance.within_linear(gap)
    })
}

/// Drop the crossings from index `from` on that lie past the linear
/// reach in space ([`afar`]): the lifted pass's, whose parameters are
/// the stored boundary's and the record's, a lift by a period being the
/// same place.
pub(super) fn keep_one_place(
    surface: Option<&SurfaceGeometry>,
    boundary: &[(HalfEdgeHandle, NurbsCurve)],
    curve_uv: &NurbsCurve,
    crossings: &mut Vec<UvCrossing>,
    from: usize,
    tolerance: &Tolerance,
) {
    let mut k = 0;
    crossings.retain(|c| {
        k += 1;
        k <= from
            || boundary
                .iter()
                .find(|(h, _)| *h == c.boundary_he)
                .is_none_or(|(_, pc)| {
                    !afar(
                        surface,
                        (pc, curve_uv),
                        (c.t_boundary, c.t_curve),
                        0.0,
                        tolerance,
                    )
                })
    });
}

/// Gauss-Newton steps refining a reported crossing.
const GAP_STEPS: usize = 24;

/// The least distance in space between a boundary pcurve and a record's
/// pcurve near a reported crossing `(t_boundary, t_curve)`, both read
/// on `surface` with `shift` taken off their u. The intersector decides
/// in the chart, to its tolerance there: where the chart stretches, a
/// near miss under it lies past the tolerance in space (a circle
/// passing 1.5e-6 inside a sphere of radius 2's rim latitude, 7.7e-7
/// rad in uv, came back as two crossings, and the clip kept the 3e-4
/// between them and dropped the circle), and a real crossing comes
/// back up to that far from its true place. Refined in space, where a
/// lift by a period is the same place, a crossing converges to no gap
/// however shallow, and a near miss stalls at its own.
pub(super) fn crossing_gap(
    surface: &SurfaceGeometry,
    (boundary, curve): (&NurbsCurve, &NurbsCurve),
    (t_boundary, t_curve): (f64, f64),
    shift: f64,
) -> f64 {
    use popkern_core::geometry::surface::SurfaceEval as _;
    let at = |c: &NurbsCurve, t: f64| {
        let p = c.point_at(t);
        surface.point_at(p.x - shift, p.y)
    };
    let (rb, rq) = (boundary.parameter_range(), curve.parameter_range());
    let clamp = |t: f64, (lo, hi): (f64, f64)| t.clamp(lo.min(hi), lo.max(hi));
    let (hb, hq) = (
        1e-7 * (rb.1 - rb.0).abs().max(1e-12),
        1e-7 * (rq.1 - rq.0).abs().max(1e-12),
    );
    let (mut tb, mut tq) = (t_boundary, t_curve);
    let mut best = f64::INFINITY;
    for _ in 0..GAP_STEPS {
        let g = at(boundary, tb) - at(curve, tq);
        best = best.min(g.length());
        if best == 0.0 || best.is_nan() {
            break;
        }
        let db = (at(boundary, tb + hb) - at(boundary, tb - hb)) * (0.5 / hb);
        let dq = (at(curve, tq + hq) - at(curve, tq - hq)) * (-0.5 / hq);
        let (a, b, c) = (db.dot(&db), db.dot(&dq), dq.dot(&dq));
        let det = a * c - b * b;
        if det.is_nan() || det.abs() <= 1e-30 * (a * c).max(1e-300) {
            break;
        }
        let (rb_, rq_) = (-db.dot(&g), -dq.dot(&g));
        let (sb, sq) = ((c * rb_ - b * rq_) / det, (a * rq_ - b * rb_) / det);
        let (nb, nq) = (
            clamp(
                tb + sb.clamp(-0.05 * (rb.1 - rb.0).abs(), 0.05 * (rb.1 - rb.0).abs()),
                rb,
            ),
            clamp(
                tq + sq.clamp(-0.05 * (rq.1 - rq.0).abs(), 0.05 * (rq.1 - rq.0).abs()),
                rq,
            ),
        );
        if nb == tb && nq == tq {
            break;
        }
        (tb, tq) = (nb, nq);
    }
    best
}

#[cfg(test)]
#[path = "clip_crossing_gap_tests.rs"]
mod tests;
