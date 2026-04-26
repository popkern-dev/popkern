//! Plane × cone oblique-slice analytic dispatch (parabola sub-case).
//!
//! When the cutting plane is parallel to exactly one cone generator
//! line, the discriminant `B² − 4AC` of the conic-equation in
//! plane-uv coordinates vanishes, and the intersection is a parabola.
//!
//! Output:
//! - `curve_3d`: exact `CurveGeometry::Parabola { vertex, linear,
//!   quadratic }` derived in closed form from the conic algebra.
//! - pcurve on cone, pcurve on plane: degree-5 NURBS interpolation
//!   over τ-samples, with knot domain `[τ_start, τ_end]`. Same
//!   pattern as the ellipse sub-case.
//! - `ParameterMap::Identity` over τ.

use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, SurfaceSide};
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, ParameterMap,
};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Oblique plane-cone slice — parabola sub-case (`B² − 4AC ≈ 0`).
/// Returns `None` for ellipse / hyperbola / degenerate cases so the
/// caller can fall through to the appropriate sibling dispatcher.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_cone_oblique_parabola(
    plane_origin: Point3,
    plane_normal: Vec3,
    plane_u: Vec3,
    plane_v: Vec3,
    plane_is_b: bool,
    cone_apex: Point3,
    cone_axis: Vec3,
    half_angle: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Reject perpendicular case (handled elsewhere).
    let cn = plane_normal.dot(&cone_axis);
    if cn.abs() >= 1.0 - tol.angular() {
        return None;
    }

    // Conic coefficients in plane-uv (s, r) — same as the ellipse
    // dispatcher. See analytic_trig_plane_cone_oblique.rs for the
    // derivation.
    let delta = plane_origin - cone_apex;
    let delta_a = delta.dot(&cone_axis);
    let mu = plane_u.dot(&cone_axis);
    let nu = plane_v.dot(&cone_axis);
    let ch = half_angle.cos();
    let ch2 = ch * ch;
    let pu2 = plane_u.dot(&plane_u);
    let pv2 = plane_v.dot(&plane_v);
    let pupv = plane_u.dot(&plane_v);
    let dpu = delta.dot(&plane_u);
    let dpv = delta.dot(&plane_v);
    let d2 = delta.dot(&delta);

    let coef_a = mu * mu - ch2 * pu2;
    let coef_b = 2.0 * (mu * nu - ch2 * pupv);
    let coef_c = nu * nu - ch2 * pv2;
    let coef_d = 2.0 * (delta_a * mu - ch2 * dpu);
    let coef_e = 2.0 * (delta_a * nu - ch2 * dpv);
    let coef_f = delta_a * delta_a - ch2 * d2;

    // Discriminant gate. The parabola sub-case is at exactly
    // `B² − 4AC = 0`; admit a small numerical band of width
    // `tol.angular()` since the conic-coefficient magnitude can be
    // O(1) for typical inputs.
    let disc = coef_b * coef_b - 4.0 * coef_a * coef_c;
    if disc.abs() > tol.angular() {
        return None;
    }

    // Eigenvector (alpha, beta) of the 2×2 quadratic-form matrix
    //   M = [[A, B/2], [B/2, C]]
    // for the non-zero eigenvalue. Since one eigenvalue is zero
    // (parabola condition), the other is `trace = A + C`. Use the
    // larger-norm row to dodge FP cancellation, mirroring the
    // perpendicular and ellipse dispatchers.
    let trace_m = coef_a + coef_c;
    if trace_m.abs() < tol.solver() {
        return None; // doubly-degenerate: line pair / single line
    }
    let m_b = coef_b * 0.5;
    let v_a = (m_b, coef_c); // proportional to eigenvector for trace
    let v_b = (coef_a, m_b);
    let n_a = (v_a.0 * v_a.0 + v_a.1 * v_a.1).sqrt();
    let n_b = (v_b.0 * v_b.0 + v_b.1 * v_b.1).sqrt();
    let (alpha, beta) = if n_a >= n_b && n_a > 1.0e-30 {
        (v_a.0 / n_a, v_a.1 / n_a)
    } else if n_b > 1.0e-30 {
        (v_b.0 / n_b, v_b.1 / n_b)
    } else {
        return None;
    };

    // Coefficient of the linear `v` term in the rotated equation
    //   trace·u² + (D·alpha + E·beta)·u + (E·alpha − D·beta)·v + F = 0
    let v_coef = coef_e * alpha - coef_d * beta;
    if v_coef.abs() < tol.solver() {
        // No `v` linear term ⇒ the equation collapses to a pair of
        // parallel lines along v, not a true parabola. Out of scope.
        return None;
    }
    let u_lin = coef_d * alpha + coef_e * beta;

    // Vertex in (u, v) rotated coords: derivative wrt u is zero at
    // the vertex.
    let u_vertex = -u_lin / (2.0 * trace_m);
    let v_vertex = -(trace_m * u_vertex * u_vertex + u_lin * u_vertex + coef_f) / v_coef;

    // Inverse rotation: (s, r) = (alpha·u − beta·v, beta·u + alpha·v).
    let s_vertex = alpha * u_vertex - beta * v_vertex;
    let r_vertex = beta * u_vertex + alpha * v_vertex;

    // 3-D parabola: P(τ) = vertex + τ · linear + τ² · quadratic.
    // From the algebra (substitute s = α·u − β·v, r = β·u + α·v
    // back into P = po + s·pu + r·pv with v(τ) = v_vertex −
    // trace·τ²/v_coef):
    //   linear    = α·pu + β·pv
    //   quadratic = (β·pu − α·pv) · trace / v_coef
    let vertex_3d = plane_origin + plane_u * s_vertex + plane_v * r_vertex;
    let linear_3d = plane_u * alpha + plane_v * beta;
    let quadratic_3d = (plane_u * beta - plane_v * alpha) * (trace_m / v_coef);

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    // τ from anchor: project (anchor − vertex) onto the
    // (linear, quadratic) basis. linear and quadratic are orthogonal
    // when (pu, pv) are orthogonal — typical case in popkern. We
    // compute τ from the linear component and validate against the
    // quadratic component.
    let lin_len_sq = linear_3d.dot(&linear_3d);
    let quad_len_sq = quadratic_3d.dot(&quadratic_3d);
    if lin_len_sq < tol.solver() || quad_len_sq < tol.solver() {
        return None;
    }
    let tau_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let dp = p - vertex_3d;
        let tau = dp.dot(&linear_3d) / lin_len_sq;
        let predicted = vertex_3d + linear_3d * tau + quadratic_3d * (tau * tau);
        let resid = predicted.distance_to(&p);
        if resid > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: resid,
            });
        }
        Ok(tau)
    };
    let tau_start = match tau_of(start, SurfaceSide::A) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let tau_end = match tau_of(end, SurfaceSide::B) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    // The eigenvector sign is ambiguous, so the natural τ may run
    // backwards relative to the (start → end) anchor direction. Flip
    // by the substitution τ → −τ, which negates `linear` and keeps
    // `quadratic` (since τ² is even). Now τ_start < τ_end.
    let (linear_3d, tau_start, tau_end) = if tau_end < tau_start {
        (linear_3d * -1.0, -tau_start, -tau_end)
    } else {
        (linear_3d, tau_start, tau_end)
    };
    if (tau_end - tau_start).abs() < tol.solver() {
        return None; // anchors coincide
    }

    let curve_3d = CurveGeometry::Parabola {
        vertex: vertex_3d,
        linear: linear_3d,
        quadratic: quadratic_3d,
        t_start: tau_start,
        t_end: tau_end,
    };

    // Build NURBS pcurves on cone and plane. Same fidelity choice
    // as the ellipse dispatcher: degree 5, n_samples == n_cps.
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(cone_axis);
    let n_samples = 33;
    let mut cone_pts: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut plane_pts: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let span = tau_end - tau_start;
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let tau = tau_start + frac * span;
        let p3 = vertex_3d + linear_3d * tau + quadratic_3d * (tau * tau);
        let dp = p3 - cone_apex;
        let v_cone = dp.dot(&cone_axis);
        let perp_vec = dp - cone_axis * v_cone;
        let u_cone = perp_vec.dot(&perp2).atan2(perp_vec.dot(&perp1));
        cone_pts.push(Point3::new(u_cone, v_cone, 0.0));
        let dp_plane = p3 - plane_origin;
        let u_plane = dp_plane.dot(&plane_u) / pu2;
        let v_plane_uv = dp_plane.dot(&plane_v) / pv2;
        plane_pts.push(Point3::new(u_plane, v_plane_uv, 0.0));
        params.push(frac);
    }
    // Unwrap atan2 jumps so the NURBS fit doesn't smooth over a
    // 2π discontinuity. Same as the ellipse dispatcher.
    for i in 1..n_samples {
        let du = cone_pts[i].x - cone_pts[i - 1].x;
        if du > std::f64::consts::PI {
            for cp in cone_pts.iter_mut().skip(i) {
                cp.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for cp in cone_pts.iter_mut().skip(i) {
                cp.x += std::f64::consts::TAU;
            }
        }
    }

    let cone_nurbs = match crate::geometry::fitting::fit_curve(&cone_pts, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    let plane_nurbs = match crate::geometry::fitting::fit_curve(&plane_pts, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    let cone_nurbs = reknot_to(cone_nurbs, tau_start, tau_end);
    let plane_nurbs = reknot_to(plane_nurbs, tau_start, tau_end);
    let pcurve_cone = Pcurve::from_nurbs(cone_nurbs);
    let pcurve_plane = Pcurve::from_nurbs(plane_nurbs);

    let plane_surf = SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        u_axis: plane_u,
        v_axis: plane_v,
    };
    let cone_surf = SurfaceGeometry::Cone {
        apex: cone_apex,
        axis: cone_axis,
        half_angle,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if plane_is_b {
        (pcurve_cone, pcurve_plane, &cone_surf, &plane_surf)
    } else {
        (pcurve_plane, pcurve_cone, &plane_surf, &cone_surf)
    };

    let triple = EdgeTriple {
        curve_3d,
        pcurve_a,
        pcurve_b,
        t_start: tau_start,
        t_end: tau_end,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::Identity,
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::PlaneCone,
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

fn reknot_to(
    mut curve: crate::geometry::nurbs_curve::NurbsCurve,
    t0: f64,
    t1: f64,
) -> crate::geometry::nurbs_curve::NurbsCurve {
    if curve.knots.is_empty() {
        return curve;
    }
    let k_lo = *curve.knots.first().unwrap();
    let k_hi = *curve.knots.last().unwrap();
    let span_old = k_hi - k_lo;
    if span_old.abs() < 1.0e-30 {
        return curve;
    }
    let span_new = t1 - t0;
    for k in &mut curve.knots {
        *k = t0 + (*k - k_lo) / span_old * span_new;
    }
    curve
}
