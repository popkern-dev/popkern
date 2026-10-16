//! Rank-3 cone fibre parameterisation for Levin's algorithm.
//!
//! When the pencil's degenerate fibre `D = λA + μB` has rank 3, it
//! is a projective cone with apex at the kernel direction. The
//! intersection `A ∩ D = A ∩ B` is a quartic curve traced by the
//! rulings of the cone — closed-form sampling instead of NURBS
//! marching's iterative Newton steps.
//!
//! ## Pipeline
//!
//! 1. Get apex `V` (kernel of `D`) — affine if `apex_w ≠ 0`,
//!    deferred to NURBS marching otherwise (the cylinder case).
//! 2. Decompose `D`'s 3 non-zero eigenvalues into one **unique-sign**
//!    `λ_axis` and two **same-sign** `λ_perp1, λ_perp2` (otherwise
//!    `D` isn't a real cone — defer).
//! 3. Parameterise the null direction
//!    `u(τ) = e_axis + cos(τ)·√(|λ_axis|/|λ_perp1|)·e_perp1
//!          + sin(τ)·√(|λ_axis|/|λ_perp2|)·e_perp2`
//!    so that `u(τ)ᵀ·Q3·u(τ) = 0` for all `τ`.
//! 4. For each `τ`, the ruling line `P(s, τ) = V + s·u(τ)` lies on
//!    `D`. Substitute into `A`: get quadratic `α(τ)·s² + β(τ)·s +
//!    γ = 0`. Solve for `s`; pick the branch tracking the start
//!    anchor's `σ`.
//! 5. Sample `(τ, s)` pairs into 3-D points and fit a degree-5 NURBS
//!    curve. Pcurves on `surf_a` / `surf_b` are sample-and-fit in UV.
//!
//! Output: `EdgeTriple` with `CurveGeometry::Nurbs` curve_3d and
//! `Pcurve::Nurbs` on both surfaces. `ParameterMap::SharedKnots`.

use super::decompose::{DegenerateShape, classify_degenerate};
use super::pencil::find_degenerate_fiber;
use super::quadric::{Mat4Sym, surface_to_quadric};
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::EdgeAnchors;
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, ParameterMap,
};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

const N_SAMPLES: usize = 129;
const FIT_DEGREE: usize = 5;

/// Try to compute `surf_a ∩ surf_b` via the rank-3 cone fibre of the
/// pencil. Returns `None` when the fibre is rank ≤ 2 (handled by
/// `dispatch::try_levin_pair`), the apex is at infinity (cylinder
/// fibre — defer to NURBS marching), or the cone has degenerate /
/// non-real signature.
pub(crate) fn try_rank3_cone_fiber(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let a_q = surface_to_quadric(surf_a)?;
    let b_q = surface_to_quadric(surf_b)?;
    let fiber = find_degenerate_fiber(&a_q, &b_q)?;
    let scale = matrix_max_abs(&a_q.0).max(matrix_max_abs(&b_q.0)).max(1.0);
    let shape = classify_degenerate(&fiber.quadric, scale);
    let DegenerateShape::Rank3Cone { apex, .. } = shape else {
        return None;
    };

    // Convert apex 4-vector to affine 3-D point. Apex at infinity
    // (w ≈ 0) means D is a cylinder — defer to NURBS marching.
    if apex[3].abs() < 1e-10 {
        return None;
    }
    let inv_w = 1.0 / apex[3];
    let v_affine = Point3::new(apex[0] * inv_w, apex[1] * inv_w, apex[2] * inv_w);

    // Extract the spatial 3×3 sub-block Q3 = D[0..3, 0..3]. The cone's
    // null directions u ∈ ℝ³ satisfy uᵀ·Q3·u = 0; the eigenstructure
    // of Q3 (NOT D's full 4×4 eigenstructure) gives the principal
    // axes for parameterising those directions.
    let mut q3 = [[0.0_f64; 3]; 3];
    for (i, row) in q3.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            *slot = fiber.quadric.0[i][j];
        }
    }
    let (q3_eig, q3_q) = jacobi_eigh_3x3(q3);

    // Identify the unique-sign eigenvalue (the cone's axis direction).
    let signs: [i32; 3] = [
        q3_eig[0].signum() as i32,
        q3_eig[1].signum() as i32,
        q3_eig[2].signum() as i32,
    ];
    let pos = signs.iter().filter(|&&s| s > 0).count();
    let neg = signs.iter().filter(|&&s| s < 0).count();
    let axis_idx = if pos == 1 && neg == 2 {
        signs.iter().position(|&s| s > 0)?
    } else if neg == 1 && pos == 2 {
        signs.iter().position(|&s| s < 0)?
    } else {
        return None; // not a real cone (definite or zero on Q3)
    };
    let perp_indices: Vec<usize> = (0..3).filter(|&i| i != axis_idx).collect();
    let perp1_idx = perp_indices[0];
    let perp2_idx = perp_indices[1];

    let lam_axis = q3_eig[axis_idx].abs();
    let lam_perp1 = q3_eig[perp1_idx].abs();
    let lam_perp2 = q3_eig[perp2_idx].abs();
    if lam_perp1 < 1e-15 || lam_perp2 < 1e-15 {
        return None;
    }

    let e_axis = Vec3::new(q3_q[0][axis_idx], q3_q[1][axis_idx], q3_q[2][axis_idx]);
    let e_perp1 = Vec3::new(q3_q[0][perp1_idx], q3_q[1][perp1_idx], q3_q[2][perp1_idx]);
    let e_perp2 = Vec3::new(q3_q[0][perp2_idx], q3_q[1][perp2_idx], q3_q[2][perp2_idx]);

    let amp_perp1 = (lam_axis / lam_perp1).sqrt();
    let amp_perp2 = (lam_axis / lam_perp2).sqrt();

    let null_dir = |tau: f64| -> Vec3 {
        e_axis + e_perp1 * (amp_perp1 * tau.pcos()) + e_perp2 * (amp_perp2 * tau.psin())
    };

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    // Map each anchor to (τ, σ) in the cone's eigenframe. Anchor in
    // eigenframe coords:
    //   (a, b, c) = σ · (1, amp_perp1·cos τ, amp_perp2·sin τ)
    // ⇒ σ = a, cos τ = b / (σ · amp_perp1), sin τ = c / (σ · amp_perp2).
    // The σ-sign-dependent recovery is essential: opposite axial
    // signs lie on different rulings (τ differs by π), even when
    // (b, c) project to the same direction.
    let to_tau_sigma = |p: Point3| -> Option<(f64, f64)> {
        let d = p - v_affine;
        let sigma = d.dot(&e_axis);
        if sigma.abs() < 1e-14 {
            return None; // anchor at apex — degenerate
        }
        let b_coord = d.dot(&e_perp1);
        let c_coord = d.dot(&e_perp2);
        let cos_t = b_coord / (sigma * amp_perp1);
        let sin_t = c_coord / (sigma * amp_perp2);
        Some((sin_t.patan2(cos_t), sigma))
    };
    let (tau_start, sigma_start) = to_tau_sigma(start)?;
    let (tau_end, _sigma_end) = to_tau_sigma(end)?;
    // Direction-of-travel: pick the shorter sweep (≤ π) — the
    // longer side may cross a non-real region of τ.
    let mut delta_tau = tau_end - tau_start;
    while delta_tau > std::f64::consts::PI {
        delta_tau -= std::f64::consts::TAU;
    }
    while delta_tau < -std::f64::consts::PI {
        delta_tau += std::f64::consts::TAU;
    }
    if delta_tau.abs() < 1e-12 {
        return None; // degenerate sweep
    }

    // Sample 3-D points along the cone fibre, branch-tracked from
    // the start anchor's σ.
    let mut sigma_prev = sigma_start;
    let mut points: Vec<Point3> = Vec::with_capacity(N_SAMPLES);
    for i in 0..N_SAMPLES {
        let frac = i as f64 / (N_SAMPLES as f64 - 1.0);
        let tau = tau_start + frac * delta_tau;
        let u = null_dir(tau);
        // Quadratic α·s² + β·s + γ = 0 from A(V + s·u) = 0.
        let u4 = [u.x, u.y, u.z, 0.0];
        let v4 = [v_affine.x, v_affine.y, v_affine.z, 1.0];
        let alpha = quadratic_form(&a_q, &u4);
        let beta = 2.0 * bilinear_form(&a_q, &v4, &u4);
        let gamma = quadratic_form(&a_q, &v4);
        if alpha.abs() < 1e-14 {
            return None;
        }
        let disc = beta * beta - 4.0 * alpha * gamma;
        if disc < -1e-9 {
            return None; // ruling doesn't intersect surf_a (real solution range broken)
        }
        let sd = disc.max(0.0).sqrt();
        let s_plus = (-beta + sd) / (2.0 * alpha);
        let s_minus = (-beta - sd) / (2.0 * alpha);
        // Pick branch closer to previous σ (continuity tracking).
        let sigma = if (s_plus - sigma_prev).abs() <= (s_minus - sigma_prev).abs() {
            s_plus
        } else {
            s_minus
        };
        let p = v_affine + u * sigma;
        points.push(p);
        sigma_prev = sigma;
    }
    // Chord-length parameterisation — gives substantially lower
    // interpolation error than uniform parameterisation for curves
    // with non-uniform curvature (which spatial quartics generally
    // have). Used identically for curve_3d and pcurve fits so the
    // verify-sample correspondence holds.
    let params = chord_length_params(&points);

    // True NURBS *interpolation* (not LSQ approximation) through the
    // analytic Levin samples. Direct LU solve avoids the
    // condition-number squaring of `fit_curve`'s normal equations,
    // which kept the residual at ~10⁻⁵ even with `num_cps == n`.
    // With degree-5 and 65 samples, inter-sample interpolation error
    // scales as O(h^6) ≈ 10⁻¹¹.
    let curve_nurbs =
        crate::geometry::fitting::interpolate_curve(&points, &params, FIT_DEGREE).ok()?;
    let curve_nurbs = reknot_to(curve_nurbs, 0.0, 1.0);
    let curve_3d = crate::geometry::curve::CurveGeometry::Nurbs(curve_nurbs);

    // Build pcurves by inverting UV at each sampled point.
    let pcurve_a = pcurve_via_uv_inversion(surf_a, &points, &params)?;
    let pcurve_b = pcurve_via_uv_inversion(surf_b, &points, &params)?;

    let triple = EdgeTriple {
        curve_3d,
        pcurve_a,
        pcurve_b,
        t_start: 0.0,
        t_end: 1.0,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::SharedKnots,
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: provenance_kind(surf_a, surf_b),
            branch: BranchSelector::OnlyBranch,
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };
    Some(super::super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

/// Chord-length parameterisation: params[i] = (∑_{k<i} |P_k − P_{k−1}|)
/// / total_length. Gives substantially better NURBS interpolation
/// quality than uniform parameterisation when the data has
/// non-uniform spatial spacing.
fn chord_length_params(points: &[Point3]) -> Vec<f64> {
    let n = points.len();
    if n < 2 {
        return vec![0.0; n];
    }
    let mut cum = vec![0.0_f64; n];
    for i in 1..n {
        cum[i] = cum[i - 1] + points[i].distance_to(&points[i - 1]);
    }
    let total = cum[n - 1];
    if total < 1e-15 {
        return (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
    }
    cum.iter().map(|&c| c / total).collect()
}

/// Cyclic Jacobi eigendecomposition for a 3×3 symmetric matrix.
/// Returns `(eigenvalues, q)` where column `j` of `q` is the unit
/// eigenvector for `eigenvalues[j]`. Same algorithm as the 4×4 in
/// `eigen.rs`, sized down for the spatial-only case.
fn jacobi_eigh_3x3(input: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = input;
    let mut q = [[0.0_f64; 3]; 3];
    for (i, row) in q.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..50 {
        let off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if off < 1e-28 {
            break;
        }
        for p in 0..2 {
            for r in (p + 1)..3 {
                let app = a[p][p];
                let arr = a[r][r];
                let apr = a[p][r];
                if apr.abs() < 1e-300 {
                    continue;
                }
                let theta = (arr - app) / (2.0 * apr);
                let t = if theta.abs() > 1e10 {
                    0.5 / theta
                } else {
                    let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
                    sign / (theta.abs() + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                a[p][p] = app - t * apr;
                a[r][r] = arr + t * apr;
                a[p][r] = 0.0;
                a[r][p] = 0.0;
                for i in (0..3).filter(|i| *i != p && *i != r) {
                    let aip = a[i][p];
                    let air = a[i][r];
                    a[i][p] = c * aip - s * air;
                    a[i][r] = s * aip + c * air;
                    a[p][i] = a[i][p];
                    a[r][i] = a[i][r];
                }
                for row in q.iter_mut() {
                    let qp = row[p];
                    let qr = row[r];
                    row[p] = c * qp - s * qr;
                    row[r] = s * qp + c * qr;
                }
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2]], q)
}

fn quadratic_form(m: &Mat4Sym, v: &[f64; 4]) -> f64 {
    let mut s = 0.0;
    for (i, vi) in v.iter().enumerate() {
        for (j, vj) in v.iter().enumerate() {
            s += vi * m.0[i][j] * vj;
        }
    }
    s
}

fn bilinear_form(m: &Mat4Sym, u: &[f64; 4], v: &[f64; 4]) -> f64 {
    let mut s = 0.0;
    for (i, ui) in u.iter().enumerate() {
        for (j, vj) in v.iter().enumerate() {
            s += ui * m.0[i][j] * vj;
        }
    }
    s
}

fn matrix_max_abs(m: &[[f64; 4]; 4]) -> f64 {
    let mut max = 0.0_f64;
    for row in m.iter() {
        for &v in row {
            max = max.max(v.abs());
        }
    }
    max
}

fn provenance_kind(a: &SurfaceGeometry, b: &SurfaceGeometry) -> AnalyticPair {
    use SurfaceGeometry::*;
    match (a, b) {
        (Cylinder { .. }, Cylinder { .. }) => AnalyticPair::CylinderCylinder,
        (Cone { .. }, Cone { .. }) => AnalyticPair::ConeCone,
        (Cylinder { .. }, Cone { .. }) | (Cone { .. }, Cylinder { .. }) => {
            AnalyticPair::CylinderCone
        }
        (Cylinder { .. }, Sphere { .. }) | (Sphere { .. }, Cylinder { .. }) => {
            AnalyticPair::CylinderSphere
        }
        (Cone { .. }, Sphere { .. }) | (Sphere { .. }, Cone { .. }) => AnalyticPair::SphereCone,
        (Sphere { .. }, Sphere { .. }) => AnalyticPair::SphereSphere,
        _ => AnalyticPair::CylinderCylinder,
    }
}

/// Build a `Pcurve::Nurbs` by inverting the surface's UV at each
/// 3-D sample. Supports `Cylinder`, `Sphere`, `Cone`, `Plane`,
/// `Torus` (rough — only fits if the 3-D curve doesn't wrap a seam).
fn pcurve_via_uv_inversion(
    surf: &SurfaceGeometry,
    points: &[Point3],
    params: &[f64],
) -> Option<Pcurve> {
    let mut uv_pts: Vec<Point3> = Vec::with_capacity(points.len());
    for p in points {
        let (u, v) = invert_uv_on_surface(surf, *p)?;
        uv_pts.push(Point3::new(u, v, 0.0));
    }
    // Unwrap u-jumps so the NURBS fit doesn't smooth over a 2π
    // discontinuity at the seam (cyl/cone/sphere/torus all wrap in u).
    unwrap_periodic_u(&mut uv_pts);
    let nurbs = crate::geometry::fitting::interpolate_curve(&uv_pts, params, FIT_DEGREE).ok()?;
    let nurbs = reknot_to(nurbs, 0.0, 1.0);
    Some(Pcurve::from_nurbs(nurbs))
}

fn invert_uv_on_surface(surf: &SurfaceGeometry, p: Point3) -> Option<(f64, f64)> {
    use SurfaceGeometry::*;
    match surf {
        Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => {
            let d = p - *origin;
            let uu = u_axis.dot(u_axis);
            let vv = v_axis.dot(v_axis);
            if uu < 1e-15 || vv < 1e-15 {
                return None;
            }
            let u = d.dot(u_axis) / uu;
            let v = d.dot(v_axis) / vv;
            Some((u, v))
        }
        Cylinder { origin, axis, .. } => {
            let d = p - *origin;
            let v = d.dot(axis);
            let perp = d - *axis * v;
            let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(*axis);
            let u = perp.dot(&perp2).patan2(perp.dot(&perp1));
            Some((u, v))
        }
        Sphere {
            center,
            radius,
            axis,
        } => {
            let d = crate::geometry::analytic_surfaces::sphere_to_local(axis, p - *center);
            let u = d.y.patan2(d.x);
            let v = (d.z / radius).clamp(-1.0, 1.0).pacos();
            Some((u, v))
        }
        Cone { apex, axis, .. } => {
            let d = p - *apex;
            let v = d.dot(axis);
            let perp = d - *axis * v;
            let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(*axis);
            let u = perp.dot(&perp2).patan2(perp.dot(&perp1));
            Some((u, v))
        }
        Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let d = p - *center;
            let v_axis = d.dot(axis);
            let perp = d - *axis * v_axis;
            let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(*axis);
            let u = perp.dot(&perp2).patan2(perp.dot(&perp1));
            let perp_len = perp.length();
            let cos_v = ((perp_len - major_radius) / minor_radius).clamp(-1.0, 1.0);
            let sin_v = (v_axis / minor_radius).clamp(-1.0, 1.0);
            let v = sin_v.patan2(cos_v);
            Some((u, v))
        }
        _ => None,
    }
}

fn unwrap_periodic_u(uv_pts: &mut [Point3]) {
    for i in 1..uv_pts.len() {
        let du = uv_pts[i].x - uv_pts[i - 1].x;
        if du > std::f64::consts::PI {
            for cp in uv_pts.iter_mut().skip(i) {
                cp.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for cp in uv_pts.iter_mut().skip(i) {
                cp.x += std::f64::consts::TAU;
            }
        }
    }
}

fn reknot_to(
    mut nurbs: crate::geometry::NurbsCurve,
    t_start: f64,
    t_end: f64,
) -> crate::geometry::NurbsCurve {
    let Some(&lo) = nurbs.knots.first() else {
        return nurbs;
    };
    let Some(&hi) = nurbs.knots.last() else {
        return nurbs;
    };
    if (hi - lo).abs() < 1e-15 {
        return nurbs;
    }
    let scale = (t_end - t_start) / (hi - lo);
    for k in nurbs.knots.iter_mut() {
        *k = t_start + (*k - lo) * scale;
    }
    nurbs
}

fn snap_endpoints(
    mut nurbs: crate::geometry::NurbsCurve,
    start: Point3,
    end: Point3,
) -> crate::geometry::NurbsCurve {
    if nurbs.control_points.is_empty() {
        return nurbs;
    }
    let n = nurbs.control_points.len();
    nurbs.control_points[0] = start;
    nurbs.control_points[n - 1] = end;
    nurbs
}
