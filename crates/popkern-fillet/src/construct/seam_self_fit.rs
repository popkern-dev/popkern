//! Same-face (seam/slit) edge curve refit (PK-177).
//!
//! A self-loop edge has no pair intersection; its curve's honest
//! home is ON the one surface. Direct-projection refit, witness-
//! based sampling, transactional acceptance (strict improvement +
//! endpoint preservation).

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::topology::{BRep, EdgeHandle};

/// Refit a SAME-FACE (seam/slit) edge's curve by direct projection
/// onto its one surface: sample the current curve, project each
/// sample (witness-based), fit degree 3 through the projections,
/// install when the worst residual strictly improves and the
/// endpoints stay put. No pair ambiguity exists for this class.
pub(in crate::construct) fn try_same_face_seam_fit(brep: &mut BRep, edge: EdgeHandle) -> bool {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::fitting::fit_curve;
    use popkern_core::geometry::surface::SurfaceEval;
    let Some(edge_data) = brep.edges.get(edge) else {
        return false;
    };
    let Some(curve) = edge_data.curve.clone() else {
        return false;
    };
    let Some(he) = brep.half_edges.get(edge_data.half_edge) else {
        return false;
    };
    let Some(face) = brep.loops.get(he.loop_).map(|l| l.face) else {
        return false;
    };
    let Some(sof) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .cloned()
    else {
        return false;
    };
    if !matches!(
        sof.geometry,
        popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
    ) {
        return false;
    }
    let band = brep
        .edges
        .get(edge)
        .and_then(|e| e.tolerance)
        .unwrap_or_else(|| brep.tolerance.linear())
        .max(1e-9)
        * 10.0;
    let measure = |c: &popkern_core::geometry::curve::CurveOnEdge| -> f64 {
        let mut worst = 0.0_f64;
        let mut prev: Option<(f64, f64)> = None;
        for i in 0..=8 {
            let t = c.t_start + (c.t_end - c.t_start) * (i as f64) / 8.0;
            let pt = c.geometry.point_at(t);
            let cold =
                popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &pt);
            let seeded = prev.and_then(|(pu, pv)| {
                popkern_core::geometry::uv_projection::project_point_to_uv_near(
                    &sof.geometry,
                    &pt,
                    pu,
                    pv,
                )
            });
            if let Some((d, uv)) = [cold, seeded]
                .into_iter()
                .flatten()
                .map(|(u, v)| (pt.distance_to(&sof.geometry.point_at(u, v)), (u, v)))
                .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal))
            {
                prev = Some(uv);
                worst = worst.max(d);
            }
        }
        worst
    };
    let before = measure(&curve);
    if before <= band.max(3.3e-5) {
        return false;
    }
    // Project 9 samples onto the surface and fit through them.
    let mut points = Vec::with_capacity(9);
    let mut params = Vec::with_capacity(9);
    let mut prev: Option<(f64, f64)> = None;
    for i in 0..=8 {
        let frac = (i as f64) / 8.0;
        let t = curve.t_start + (curve.t_end - curve.t_start) * frac;
        let pt = curve.geometry.point_at(t);
        let cold = popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &pt);
        let seeded = prev.and_then(|(pu, pv)| {
            popkern_core::geometry::uv_projection::project_point_to_uv_near(
                &sof.geometry,
                &pt,
                pu,
                pv,
            )
        });
        let Some((_, uv)) = [cold, seeded]
            .into_iter()
            .flatten()
            .map(|(u, v)| (pt.distance_to(&sof.geometry.point_at(u, v)), (u, v)))
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal))
        else {
            return false;
        };
        prev = Some(uv);
        points.push(sof.geometry.point_at(uv.0, uv.1));
        params.push(frac);
    }
    let Ok(fitted) = fit_curve(&points, &params, 3, points.len()) else {
        return false;
    };
    let (t0, t1) = fitted.parameter_range();
    let candidate = popkern_core::geometry::curve::CurveOnEdge {
        geometry: CurveGeometry::Nurbs(fitted),
        t_start: t0,
        t_end: t1,
    };
    // Endpoints must stay at the edge's vertices.
    let p_start = curve.geometry.point_at(curve.t_start);
    let p_end = curve.geometry.point_at(curve.t_end);
    let e0 = candidate.geometry.point_at(t0).distance_to(&p_start);
    let e1 = candidate.geometry.point_at(t1).distance_to(&p_end);
    if e0.max(e1) > (p_start.distance_to(&p_end) * 0.05).max(1e-6) {
        return false;
    }
    let after = measure(&candidate);
    if after >= before {
        return false;
    }
    if let Some(edge_mut) = brep.edges.get_mut(edge) {
        edge_mut.curve = Some(candidate);
        if std::env::var_os("PK_RECONCILE_TRACE").is_some() {
            eprintln!("[fit-tol] {edge:?} same-face seam refit {before:.3e} -> {after:.3e}");
        }
        return true;
    }
    false
}
