//! Diagnostic metrics for comparing vertex-blend surface
//! constructions (polar-NURBS vs rational-Gregory). Split from
//! `gregory_diagnostic_tests.rs` to keep both files under the
//! 500-line structural cap.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};

use super::gregory::{extract_triangular_inner_cps, fit_cubic_bezier_through};
use super::surface::BlendSurfaceResult;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod};
use popkern_core::geometry::gregory::TriangularGregory;

#[derive(Debug)]
pub(super) struct DiagnosticMetrics {
    pub method: CornerMethod,
    pub continuity_claimed: ContinuityLevel,
    pub max_boundary_deviation: f64,
    pub max_interior_gaussian_curvature: f64,
    pub mean_interior_gaussian_curvature: f64,
    pub max_boundary_normal_curvature_deviation: f64,
}

impl std::fmt::Display for DiagnosticMetrics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "  method                           = {:?}", self.method)?;
        writeln!(
            f,
            "  continuity_claimed               = {:?}",
            self.continuity_claimed
        )?;
        writeln!(
            f,
            "  max_boundary_deviation           = {:.3e}",
            self.max_boundary_deviation
        )?;
        writeln!(
            f,
            "  max_interior_gaussian_curvature  = {:.3e}",
            self.max_interior_gaussian_curvature
        )?;
        writeln!(
            f,
            "  mean_interior_gaussian_curvature = {:.3e}",
            self.mean_interior_gaussian_curvature
        )?;
        writeln!(
            f,
            "  max_boundary_normal_κ_deviation  = {:.3e}  (§2.4 < 1e-2)",
            self.max_boundary_normal_curvature_deviation
        )
    }
}

/// Gaussian curvature `K = (LN − M²) / (EG − F²)`.
fn gaussian_curvature(surface: &SurfaceGeometry, u: f64, v: f64) -> f64 {
    let pu = surface.partial_u(u, v);
    let pv = surface.partial_v(u, v);
    let e = pu.dot(&pu);
    let f = pu.dot(&pv);
    let g = pv.dot(&pv);
    let det_first = e * g - f * f;
    if det_first.abs() < 1e-15 {
        return 0.0;
    }
    let normal = pu.cross(&pv);
    let n_len = normal.length();
    if n_len < 1e-15 {
        return 0.0;
    }
    let n_hat = normal * (1.0 / n_len);
    let puu = surface.second_partial_uu(u, v);
    let puv = surface.second_partial_uv(u, v);
    let pvv = surface.second_partial_vv(u, v);
    let l = puu.dot(&n_hat);
    let m = puv.dot(&n_hat);
    let nn = pvv.dot(&n_hat);
    (l * nn - m * m) / det_first
}

pub(super) fn compute_metrics(
    result: &BlendSurfaceResult,
    input_boundary_curves: &[NurbsCurve],
    avg_radius: f64,
) -> DiagnosticMetrics {
    let surface = &result.surface;
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let (u0, u1) = if u0.is_finite() && u1.is_finite() {
        (u0, u1)
    } else {
        (0.0, 1.0)
    };
    let (v0, v1) = if v0.is_finite() && v1.is_finite() {
        (v0, v1)
    } else {
        (0.0, 1.0)
    };

    let mut max_boundary_dev: f64 = 0.0;
    for curve in input_boundary_curves {
        let (t0, t1) = curve.parameter_range();
        for k in 0..=4 {
            let t = t0 + (t1 - t0) * (k as f64) / 4.0;
            let p_curve = curve.point_at(t);
            if let Some((u, v)) = surface.closest_point(&p_curve) {
                let p_surf = surface.point_at(u, v);
                let d = p_curve.distance_to(&p_surf);
                max_boundary_dev = max_boundary_dev.max(d);
            }
        }
    }

    let grid_n = 9;
    let mut max_k: f64 = 0.0;
    let mut sum_k = 0.0;
    let mut n_samples = 0usize;
    for i in 1..(grid_n - 1) {
        let u = u0 + (u1 - u0) * (i as f64) / (grid_n - 1) as f64;
        for j in 1..(grid_n - 1) {
            let v = v0 + (v1 - v0) * (j as f64) / (grid_n - 1) as f64;
            if matches!(
                surface,
                SurfaceGeometry::Gregory(GregoryPatch::Triangular(_))
            ) && u + v > 1.0
            {
                continue;
            }
            let k = gaussian_curvature(surface, u, v).abs();
            if !k.is_finite() {
                continue;
            }
            max_k = max_k.max(k);
            sum_k += k;
            n_samples += 1;
        }
    }
    let mean_k = if n_samples > 0 {
        sum_k / (n_samples as f64)
    } else {
        0.0
    };

    let target_kappa = 1.0 / avg_radius;
    let mut max_kappa_dev: f64 = 0.0;
    for curve in input_boundary_curves {
        let (t0, t1) = curve.parameter_range();
        for k in 0..=4 {
            let t = t0 + (t1 - t0) * (k as f64) / 4.0;
            let p = curve.point_at(t);
            let Some((u, v)) = surface.closest_point(&p) else {
                continue;
            };
            let pu = surface.partial_u(u, v);
            let pv = surface.partial_v(u, v);
            let Some(t_dir) = curve.derivative_at(t, 1).normalized() else {
                continue;
            };
            let e = pu.dot(&pu);
            let f = pu.dot(&pv);
            let g = pv.dot(&pv);
            let du = pu.dot(&t_dir);
            let dv = pv.dot(&t_dir);
            let det = e * g - f * f;
            if det.abs() < 1e-15 {
                continue;
            }
            let a = (g * du - f * dv) / det;
            let b = (e * dv - f * du) / det;
            let normal = pu.cross(&pv);
            let n_len = normal.length();
            if n_len < 1e-15 {
                continue;
            }
            let n_hat = normal * (1.0 / n_len);
            let puu = surface.second_partial_uu(u, v);
            let puv = surface.second_partial_uv(u, v);
            let pvv = surface.second_partial_vv(u, v);
            let ii =
                a * a * puu.dot(&n_hat) + 2.0 * a * b * puv.dot(&n_hat) + b * b * pvv.dot(&n_hat);
            let i_first = a * a * e + 2.0 * a * b * f + b * b * g;
            if i_first.abs() < 1e-15 {
                continue;
            }
            let kappa_n = (ii / i_first).abs();
            let denom = kappa_n.max(target_kappa);
            let dev = (kappa_n - target_kappa).abs() / denom;
            if dev.is_finite() {
                max_kappa_dev = max_kappa_dev.max(dev);
            }
        }
    }

    DiagnosticMetrics {
        method: result.method,
        continuity_claimed: result.continuity,
        max_boundary_deviation: max_boundary_dev,
        max_interior_gaussian_curvature: max_k,
        mean_interior_gaussian_curvature: mean_k,
        max_boundary_normal_curvature_deviation: max_kappa_dev,
    }
}

/// Build a rational Gregory patch unconditionally — bypasses the
/// production path's uniform-radii gate so the diagnostic can
/// compare both construction paths on mixed-radius fixtures.
pub(super) fn build_unchecked_rational_gregory(
    boundary_curves: &[NurbsCurve],
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
) -> Option<BlendSurfaceResult> {
    if boundary_curves.len() != 3 || fillet_surfaces.len() != 3 {
        return None;
    }
    let radii: Vec<f64> = fillet_surfaces
        .iter()
        .filter_map(|fr| match &fr.surface {
            SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
            SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    if radii.len() != 3 {
        return None;
    }
    let avg_radius = radii.iter().sum::<f64>() / 3.0;
    let beziers: Vec<[Point3; 4]> = boundary_curves
        .iter()
        .map(fit_cubic_bezier_through)
        .collect();
    let inner_cps = extract_triangular_inner_cps(
        boundary_curves,
        fillet_surfaces,
        avg_radius,
        vertex_position,
    )?;
    let bez_array: [[Point3; 4]; 3] = [beziers[0], beziers[1], beziers[2]];
    let patch = TriangularGregory::from_boundary_beziers(bez_array, inner_cps, 1e-3).ok()?;
    Some(BlendSurfaceResult {
        surface: SurfaceGeometry::Gregory(GregoryPatch::Triangular(Box::new(patch))),
        boundary_curves: boundary_curves.to_vec(),
        method: CornerMethod::Gregory,
        continuity: ContinuityLevel::G1,
    })
}
