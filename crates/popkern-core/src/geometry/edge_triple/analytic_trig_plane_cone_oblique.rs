//! Plane × cone oblique-slice analytic dispatch (ellipse sub-case).
//!
//! When the cutting plane is neither perpendicular to the cone axis
//! (handled by `analytic_trig_cone`) nor parallel to a generator, the
//! intersection is a 3-D ellipse. This dispatcher computes the ellipse
//! parameters in closed form from the cone-plane algebra and emits an
//! exact `CurveGeometry::Ellipse` for `curve_3d`.
//!
//! Pcurves on the cone and plane are sampled and degree-3 NURBS-fit
//! over the ellipse's natural theta-domain, since the cone-surface
//! u-coordinate as a function of theta is `atan2(rational-sinusoid,
//! rational-sinusoid)` — analytically representable but not in any
//! existing pcurve variant. NURBS sample-and-fit gives ~1e-12 residual
//! on a half-arc, well within `tol.linear()`.
//!
//! `ParameterMap::Identity` holds: the NURBS pcurves are built with
//! knot domain `[theta_start, theta_end]` matching the curve_3d.
//!
//! Parabola and hyperbola sub-cases are deferred to stages 2b.13 and
//! 2b.14 respectively, which add `CurveGeometry::Parabola` /
//! `Hyperbola` variants.

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

/// Oblique plane-cone slice. Returns `None` for the perpendicular
/// case (handled elsewhere), the apex-passing case (degenerate
/// 2-line / point), and the parabola / hyperbola sub-cases (deferred
/// to subsequent stages). On the ellipse case, returns an exact
/// `CurveGeometry::Ellipse` for `curve_3d`.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_cone_oblique_ellipse(
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
    // Reject perpendicular (perpendicular dispatcher owns it).
    let cn = plane_normal.dot(&cone_axis);
    if cn.abs() >= 1.0 - tol.angular() {
        return None;
    }

    // Classify the conic. Working in plane-uv coords (s, r):
    //   point on plane: P = po + s · pu + r · pv
    //   axial component: axial(s, r) = delta_a + s mu + r nu
    //     where delta_a = (po - apex)·a, mu = pu·a, nu = pv·a
    //   cone equation: axial² = cos²theta · |P - apex|²
    //
    // Expanding gives a quadratic A·s² + B·s·r + C·r² + D·s + E·r + F = 0.
    // Discriminant B² − 4AC tells the conic type.
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
    if disc >= -tol.angular() {
        // Parabola (disc ≈ 0) or hyperbola (disc > 0). Defer.
        return None;
    }

    // Ellipse centre in plane-uv coords. Solve
    //   2A·s_c +   B·r_c = −D
    //     B·s_c + 2C·r_c = −E
    let det_centre = 4.0 * coef_a * coef_c - coef_b * coef_b;
    if det_centre.abs() < tol.solver() {
        return None; // numerically degenerate
    }
    let s_c = (coef_b * coef_e - 2.0 * coef_c * coef_d) / det_centre;
    let r_c = (coef_b * coef_d - 2.0 * coef_a * coef_e) / det_centre;

    // Constant term after centring: F' = F + (D·s_c + E·r_c) / 2.
    let f_prime = coef_f + 0.5 * (coef_d * s_c + coef_e * r_c);

    // Ellipse principal axes: eigendecompose the 2×2 matrix
    //   M = [[A, B/2], [B/2, C]]
    // — eigenvalues lam1 ≥ lam2, eigenvectors orthonormal w1, w2.
    let m_a = coef_a;
    let m_c = coef_c;
    let m_b = coef_b * 0.5;
    let trace = m_a + m_c;
    let disc_eig = ((m_a - m_c) * (m_a - m_c) + 4.0 * m_b * m_b)
        .max(0.0)
        .sqrt();
    let lam1 = (trace + disc_eig) * 0.5;
    let lam2 = (trace - disc_eig) * 0.5;

    // Eigenvector for lam1: (M − lam1·I) w = 0. Row 1 → w ∝ (m_b, lam1 − m_a);
    // row 2 → w ∝ (lam1 − m_c, m_b). Pick the larger-norm formula to
    // avoid the same FP cancellation we hit in `decompose_2x2`.
    let w_a = (m_b, lam1 - m_a);
    let w_b = (lam1 - m_c, m_b);
    let n_a = (w_a.0 * w_a.0 + w_a.1 * w_a.1).sqrt();
    let n_b = (w_b.0 * w_b.0 + w_b.1 * w_b.1).sqrt();
    let w1 = if n_a >= n_b && n_a > 1.0e-30 {
        (w_a.0 / n_a, w_a.1 / n_a)
    } else if n_b > 1.0e-30 {
        (w_b.0 / n_b, w_b.1 / n_b)
    } else {
        (1.0, 0.0)
    };
    let w2 = (-w1.1, w1.0); // 90° rotation; orthonormal partner.

    // Semi-axis lengths from lam·a² = −F' (assumes ellipse so signs align).
    if f_prime * lam1 >= 0.0 || f_prime * lam2 >= 0.0 {
        return None; // not a real ellipse — defer
    }
    let semi_a = (-f_prime / lam1).sqrt();
    let semi_b = (-f_prime / lam2).sqrt();

    // 3-D ellipse parameters: centre and semi-axis vectors.
    let centre_3d = plane_origin + plane_u * s_c + plane_v * r_c;
    let u_basis = (plane_u * w1.0 + plane_v * w1.1) * semi_a;
    let v_basis = (plane_u * w2.0 + plane_v * w2.1) * semi_b;

    // Anchor → ellipse-theta. Decompose anchor − centre into (U, V) basis.
    // U and V are orthogonal (eigenvectors of a real symmetric matrix
    // are orthogonal), so cos theta · |U| and sin theta · |V| recover the
    // angle. Validity: cos²+sin² should be ~1 if the anchor really is
    // on the ellipse.
    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    let u_len_sq = u_basis.dot(&u_basis);
    let v_len_sq = v_basis.dot(&v_basis);
    if u_len_sq < tol.solver() || v_len_sq < tol.solver() {
        return None;
    }

    let theta_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let dp = p - centre_3d;
        let cos_part = dp.dot(&u_basis) / u_len_sq;
        let sin_part = dp.dot(&v_basis) / v_len_sq;
        let r2 = cos_part * cos_part + sin_part * sin_part;
        if (r2 - 1.0).abs() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: r2.sqrt() - 1.0,
            });
        }
        Ok(sin_part.atan2(cos_part))
    };
    let theta_start = match theta_of(start, SurfaceSide::A) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let theta_end = match theta_of(end, SurfaceSide::B) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let theta_end = if theta_end < theta_start {
        theta_end + std::f64::consts::TAU
    } else {
        theta_end
    };

    let curve_3d = CurveGeometry::Ellipse {
        center: centre_3d,
        u_basis,
        v_basis,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    // Build NURBS pcurves on cone and plane by sampling 17 theta values
    // and fitting. The cone u(theta) = atan2(rational sinusoid,
    // rational sinusoid) — smooth but with significant curvature
    // (factor-of-2 slope swing for shallow cones). Use degree 5
    // with n_samples == n_cps so the LSQ collapses to exact
    // interpolation at the samples; off-sample error is bounded
    // by the 6th derivative of u(theta) over the span and is
    // < 1e-12 for the half-arcs we encounter.
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(cone_axis);
    let n_samples = 33;
    let mut cone_pts: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut plane_pts: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let span = theta_end - theta_start;
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let theta = theta_start + frac * span;
        let p3 = centre_3d + u_basis * theta.cos() + v_basis * theta.sin();
        // Cone uv: u = atan2(perp_proj·perp2, perp_proj·perp1),
        //           v = (p3 − apex) · axis.
        let dp = p3 - cone_apex;
        let v_cone = dp.dot(&cone_axis);
        let perp_vec = dp - cone_axis * v_cone;
        let u_cone = perp_vec.dot(&perp2).atan2(perp_vec.dot(&perp1));
        cone_pts.push(Point3::new(u_cone, v_cone, 0.0));
        // Plane uv: u_plane = (p3 − po) · pu / |pu|², v similarly.
        let dp_plane = p3 - plane_origin;
        let u_plane = dp_plane.dot(&plane_u) / pu2;
        let v_plane = dp_plane.dot(&plane_v) / pv2;
        plane_pts.push(Point3::new(u_plane, v_plane, 0.0));
        params.push(frac);
    }
    // Unwrap the cone u-coordinate: atan2 jumps by ±2π at the seam,
    // which would create a discontinuity that NURBS-fit smooths out
    // (giving large interpolation error near the jump). Detect a
    // jump > π between consecutive samples and add ±2π to keep the
    // sequence monotonic-ish.
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
    // Re-knot the NURBS to live on [theta_start, theta_end] so
    // ParameterMap::Identity holds against curve_3d's theta.
    let cone_nurbs = reknot_to(cone_nurbs, theta_start, theta_end);
    let plane_nurbs = reknot_to(plane_nurbs, theta_start, theta_end);
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
        t_start: theta_start,
        t_end: theta_end,
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

/// Re-knot a NURBS so its parameter range is `[t0, t1]`.
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
