//! Survival-rule measurement for pcurve recomputation (see the
//! caller in `pcurve_recompute`): a refit may only replace an
//! existing pcurve when it satisfies the trim contract at least as
//! well.

/// Worst fraction-matched trim-contract deviation of a pcurve against
/// the half-edge's traversal of the 3-D curve: `S(pc(f))` vs
/// `C(t_lo + f·(t_hi − t_lo))` over interior fractions. Used by the
/// recompute survival rule to compare an existing pcurve with a
/// freshly fitted candidate.
pub(in crate::topology) fn pcurve_fraction_contract_error(
    surface: &crate::geometry::surface::SurfaceGeometry,
    pcurve: &crate::geometry::pcurve::Pcurve,
    curve: &crate::geometry::curve::CurveOnEdge,
    (t_lo, t_hi): (f64, f64),
) -> f64 {
    use crate::geometry::curve::CurveEval;
    use crate::geometry::surface::SurfaceEval;
    let eval: &dyn CurveEval = &curve.geometry;
    let mut worst = 0.0_f64;
    for k in 0..9 {
        let f = 0.05 + 0.90 * k as f64 / 8.0;
        let uv = pcurve.point_at(pcurve.t_start + f * (pcurve.t_end - pcurve.t_start));
        let s = surface.point_at(uv.x, uv.y);
        let c = eval.point_at(t_lo + f * (t_hi - t_lo));
        worst = worst.max(s.distance_to(&c));
    }
    worst
}
