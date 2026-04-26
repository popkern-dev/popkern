//! Plane × cone oblique-slice analytic dispatch (hyperbola sub-case).
//!
//! When the cutting plane is steeper than the cone's generator
//! angle, the discriminant `B² − 4AC` of the conic-equation in
//! plane-uv coordinates is positive, and the intersection is a
//! hyperbola. Each nappe of the cone gives one branch; for the
//! standard half-cone (`v ≥ 0`) we only see one branch through any
//! given anchor pair.
//!
//! Output:
//! - `curve_3d`: exact `CurveGeometry::Hyperbola { center,
//!   semi_major, semi_minor }` derived in closed form from the
//!   conic algebra.
//! - pcurve on cone, pcurve on plane: degree-5 NURBS interpolation
//!   over τ-samples (same fidelity choice as the ellipse and
//!   parabola dispatchers).
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

/// Oblique plane-cone slice — hyperbola sub-case (`B² − 4AC > 0`).
/// Returns `None` for ellipse / parabola / degenerate cases.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_cone_oblique_hyperbola(
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
    let cn = plane_normal.dot(&cone_axis);
    if cn.abs() >= 1.0 - tol.angular() {
        return None;
    }

    // Same conic-coefficients setup as the ellipse / parabola
    // dispatchers. See analytic_trig_plane_cone_oblique.rs for the
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

    let disc = coef_b * coef_b - 4.0 * coef_a * coef_c;
    if disc <= tol.angular() {
        // Ellipse (disc < 0) or parabola (disc ≈ 0). Defer.
        return None;
    }

    // Hyperbola centre (same gradient = 0 system as ellipse).
    let det_centre = 4.0 * coef_a * coef_c - coef_b * coef_b;
    if det_centre.abs() < tol.solver() {
        return None;
    }
    let s_c = (coef_b * coef_e - 2.0 * coef_c * coef_d) / det_centre;
    let r_c = (coef_b * coef_d - 2.0 * coef_a * coef_e) / det_centre;
    let f_prime = coef_f + 0.5 * (coef_d * s_c + coef_e * r_c);

    // Eigendecompose the 2×2 quadratic-form matrix M.
    let m_a = coef_a;
    let m_c = coef_c;
    let m_b = coef_b * 0.5;
    let trace = m_a + m_c;
    let disc_eig = ((m_a - m_c) * (m_a - m_c) + 4.0 * m_b * m_b)
        .max(0.0)
        .sqrt();
    let lam1 = (trace + disc_eig) * 0.5;
    let lam2 = (trace - disc_eig) * 0.5;

    // For a hyperbola, the eigenvalues have opposite signs. Identify
    // which is the "positive" (cosh-arm) eigenvalue and which is the
    // "negative" (sinh-arm) eigenvalue based on the sign of `f_prime`.
    //
    // Shifted equation: lam1 · u² + lam2 · v² = −f_prime.
    // For a hyperbola opening in the eigen-direction with eigenvalue
    // sign opposite to `f_prime`:
    //   semi_axis² = −f_prime / lam_real    (positive)
    //   imag_axis² =  f_prime / lam_imag    (positive)
    if (lam1 > 0.0) == (lam2 > 0.0) {
        return None; // both same sign — not a hyperbola
    }

    // Pick `lam_real` as the eigenvalue with the same sign as `−f_prime`,
    // so `−f_prime / lam_real > 0`. The other eigenvalue is `lam_imag`.
    let (lam_real, lam_imag, real_is_lam1) = if (-f_prime).signum() == lam1.signum() {
        (lam1, lam2, true)
    } else {
        (lam2, lam1, false)
    };
    let semi_real_sq = -f_prime / lam_real;
    let semi_imag_sq = f_prime / lam_imag;
    if semi_real_sq <= 0.0 || semi_imag_sq <= 0.0 {
        return None;
    }
    let semi_real = semi_real_sq.sqrt();
    let semi_imag = semi_imag_sq.sqrt();

    // Eigenvector for `lam_real`. Same robustness pattern as the
    // ellipse dispatcher: try both row formulas, pick the one with
    // larger pre-normalisation magnitude.
    let lam = if real_is_lam1 { lam1 } else { lam2 };
    let w_a = (m_b, lam - m_a);
    let w_b = (lam - m_c, m_b);
    let n_a = (w_a.0 * w_a.0 + w_a.1 * w_a.1).sqrt();
    let n_b = (w_b.0 * w_b.0 + w_b.1 * w_b.1).sqrt();
    let w_real = if n_a >= n_b && n_a > 1.0e-30 {
        (w_a.0 / n_a, w_a.1 / n_a)
    } else if n_b > 1.0e-30 {
        (w_b.0 / n_b, w_b.1 / n_b)
    } else {
        return None;
    };
    let w_imag = (-w_real.1, w_real.0);

    // 3-D hyperbola parameters: centre and axis vectors.
    let center_3d = plane_origin + plane_u * s_c + plane_v * r_c;
    let semi_major = (plane_u * w_real.0 + plane_v * w_real.1) * semi_real;
    let semi_minor = (plane_u * w_imag.0 + plane_v * w_imag.1) * semi_imag;

    // τ from anchor: solve `anchor - center = cosh(τ)·semi_major +
    // sinh(τ)·semi_minor`. Project onto each axis (orthogonal when
    // (pu, pv) are orthogonal):
    //   cosh(τ) = (anchor − center) · semi_major / |semi_major|²
    //   sinh(τ) =                  · semi_minor / |semi_minor|²
    // Then τ = atanh(sinh / cosh) (since asinh works on signed
    // values and cosh > 0 always; atanh of tanh = τ).
    let major_len_sq = semi_major.dot(&semi_major);
    let minor_len_sq = semi_minor.dot(&semi_minor);
    if major_len_sq < tol.solver() || minor_len_sq < tol.solver() {
        return None;
    }

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    // Determine branch sign from the start anchor: if the cosh
    // component is negative, the anchor is on the opposite branch,
    // so we must negate `semi_major` to point at that branch.
    let dp_start = start - center_3d;
    let cosh_start_raw = dp_start.dot(&semi_major) / major_len_sq;
    let (semi_major, branch_flip) = if cosh_start_raw < 0.0 {
        (semi_major * -1.0, true)
    } else {
        (semi_major, false)
    };
    let _ = branch_flip;
    let major_len_sq = semi_major.dot(&semi_major);

    let tau_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let dp = p - center_3d;
        let cosh_tau = dp.dot(&semi_major) / major_len_sq;
        let sinh_tau = dp.dot(&semi_minor) / minor_len_sq;
        if cosh_tau < 1.0 - tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: cosh_tau - 1.0,
            });
        }
        // cosh − sinh validity: cosh² − sinh² should be 1.
        let resid = (cosh_tau * cosh_tau - sinh_tau * sinh_tau - 1.0).abs();
        if resid > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: resid.sqrt(),
            });
        }
        // τ = asinh(sinh_tau) — preserves sign and magnitude exactly.
        Ok(sinh_tau.asinh())
    };
    let tau_start = match tau_of(start, SurfaceSide::A) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let tau_end = match tau_of(end, SurfaceSide::B) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    // Direction handling: hyperbola is monotonic in τ along one
    // branch, so naturally tau_end ≠ tau_start. If tau_end < tau_start,
    // negate semi_minor (which negates sinh, flipping τ → −τ; cosh
    // is even so semi_major is unchanged).
    let (semi_minor, tau_start, tau_end) = if tau_end < tau_start {
        (semi_minor * -1.0, -tau_start, -tau_end)
    } else {
        (semi_minor, tau_start, tau_end)
    };
    if (tau_end - tau_start).abs() < tol.solver() {
        return None;
    }

    let curve_3d = CurveGeometry::Hyperbola {
        center: center_3d,
        semi_major,
        semi_minor,
        t_start: tau_start,
        t_end: tau_end,
    };

    // NURBS pcurves — degree 5, n_samples == n_cps for exact
    // interpolation at the samples.
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(cone_axis);
    let n_samples = 33;
    let mut cone_pts: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut plane_pts: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let span = tau_end - tau_start;
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let tau = tau_start + frac * span;
        let p3 = center_3d + semi_major * tau.cosh() + semi_minor * tau.sinh();
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
            cardinality: None,
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
