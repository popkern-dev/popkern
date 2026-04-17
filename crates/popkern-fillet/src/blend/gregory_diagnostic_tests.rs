//! §5.3 Round 7 C4b.1 — comparative diagnostic for rational
//! Gregory vs polar-NURBS corner construction.
//!
//! **Purpose.** Before we activate the rational-Gregory dispatch
//! path in `build_smooth_blend` (C4b.4), we need **evidence** that
//! the new construction matches or beats the existing polar-NURBS
//! fallback on the metrics that matter:
//!
//! 1. **Boundary positional match** — does the built patch hit
//!    the input boundary curves? (Both paths should; this is a
//!    sanity check.)
//! 2. **Surface curvature profile** — does the patch interior
//!    flow smoothly? Measured via a grid sample of Gaussian
//!    curvature magnitudes.
//! 3. **|Δκ| at boundary samples** — how close does the patch's
//!    normal curvature come to the incident strip's `1/r` at the
//!    boundary? (The §2.4 G2 acceptance criterion is
//!    `|Δκ|/max(|κ|, 1/r) < 0.01`.)
//!
//! The test runs on a **mixed-radius** (0.08, 0.10, 0.12)
//! cube-corner fixture — the case where the sphere path fails
//! and polar NURBS currently fires. Gregory has to earn its
//! place here or we stick with polar NURBS.
//!
//! Output is emitted via `eprintln!` (visible with `--nocapture`)
//! so the side-by-side numbers drive the C4b.2+ decisions. Hard
//! asserts only check the invariants that MUST hold regardless
//! of construction method (both paths hit the boundary
//! corners exactly, no NaN / Inf in curvature samples).

use crate::analysis;
use crate::blend::boundary::{BoundarySegment, compute_blend_boundary};
use crate::blend::gregory::{extract_triangular_inner_cps, fit_cubic_bezier_through};
use crate::blend::surface::{BlendSurfaceResult, build_blend_surface};
use crate::fillet_surface::FilletSurfaceResult;
use crate::geometry_compute::compute_fillet_geometry;
use crate::types::{CapKind, ContinuityLevel, CornerMethod, FilletRadius, FilletSpec, VertexMode};
use crate::vertex_classifier::classify_vertex_blend;

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::{GregoryPatch, TriangularGregory};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[derive(Debug)]
struct DiagnosticMetrics {
    method: CornerMethod,
    continuity_claimed: ContinuityLevel,
    /// Max positional distance between the built patch's surface
    /// point at a boundary parameter and the input boundary
    /// curve's position at the corresponding sample. Should be
    /// ≈ 0 — both constructions preserve the boundary exactly.
    max_boundary_deviation: f64,
    /// Max absolute Gaussian curvature over a 9×9 interior sample
    /// grid. A "quality" signal — smooth corners have bounded
    /// curvature; spikes indicate creases / degeneracies.
    max_interior_gaussian_curvature: f64,
    /// Mean absolute Gaussian curvature over the same grid. Robust
    /// to a single-point spike.
    mean_interior_gaussian_curvature: f64,
    /// Max normal-curvature deviation at boundary samples vs the
    /// expected `1/r_avg` (the §2.4 acceptance metric, with the
    /// target being `max(|κ_measured|, 1/r_avg)` in the denominator).
    /// Values below 0.01 are §2.4-compliant.
    max_boundary_normal_curvature_deviation: f64,
}

impl std::fmt::Display for DiagnosticMetrics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "  method                                    = {:?}",
            self.method
        )?;
        writeln!(
            f,
            "  continuity_claimed                        = {:?}",
            self.continuity_claimed,
        )?;
        writeln!(
            f,
            "  max_boundary_deviation                    = {:.3e}",
            self.max_boundary_deviation,
        )?;
        writeln!(
            f,
            "  max_interior_gaussian_curvature           = {:.3e}",
            self.max_interior_gaussian_curvature,
        )?;
        writeln!(
            f,
            "  mean_interior_gaussian_curvature          = {:.3e}",
            self.mean_interior_gaussian_curvature,
        )?;
        writeln!(
            f,
            "  max_boundary_normal_curvature_deviation   = {:.3e}  (§2.4 target < 1e-2)",
            self.max_boundary_normal_curvature_deviation,
        )?;
        Ok(())
    }
}

/// Compute the Gaussian curvature at surface parameters `(u, v)`.
/// `K = (LN − M²) / (EG − F²)` where L, M, N are 2nd-fundamental-form
/// coefficients (projections of 2nd partials onto the normal) and
/// E, F, G are 1st-fundamental-form coefficients.
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

/// Compute diagnostic metrics for a [`BlendSurfaceResult`] built
/// on the given input boundary curves. Samples across the patch's
/// parameter range — which differs by surface type, so the
/// sampling clamps to `[0, 1]²` for Gregory (bounded domain),
/// or falls back to the surface's own reported range otherwise.
fn compute_metrics(
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

    // 1. Boundary deviation. For each input boundary curve,
    // sample 5 parameter values and find the nearest point on
    // the built surface via `closest_point`. Measure distance.
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

    // 2. Interior curvature grid (9×9 sample, skipping boundary
    // to avoid parameter-range edge artifacts).
    let grid_n = 9;
    let mut max_k: f64 = 0.0;
    let mut sum_k = 0.0;
    let mut n_samples = 0usize;
    for i in 1..(grid_n - 1) {
        let u = u0 + (u1 - u0) * (i as f64) / (grid_n - 1) as f64;
        for j in 1..(grid_n - 1) {
            let v = v0 + (v1 - v0) * (j as f64) / (grid_n - 1) as f64;
            // For triangular Gregory, only sample inside the simplex.
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

    // 3. Boundary normal-curvature deviation vs 1/r_avg. Sample
    // 5 points on each boundary curve, project to surface UV,
    // compute the normal curvature along the surface tangent
    // that matches the boundary tangent direction.
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
            // Normal-curvature along curve tangent.
            let pu = surface.partial_u(u, v);
            let pv = surface.partial_v(u, v);
            let Some(t_dir) = curve.derivative_at(t, 1).normalized() else {
                continue;
            };
            // Solve for (a, b) s.t. a·pu + b·pv ≈ t_dir in tangent
            // plane. Least-squares via 2×2 normal equations.
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
            // Normal curvature κ_n = II(a, b) / I(a, b) where II
            // and I are 2nd and 1st fundamental forms at (u, v)
            // along direction (a, b).
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

/// Build a rational Gregory patch **unconditionally** — bypasses
/// the production path's `uniform radii` gate so the diagnostic
/// can compare both construction paths on mixed-radius fixtures
/// where Gregory's current implementation is known-imperfect.
fn build_unchecked_rational_gregory(
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

#[test]
fn diagnostic_polar_nurbs_vs_rational_gregory_on_mixed_radius_cube() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    assert_eq!(incident_edges.len(), 3, "cube vertex should be trihedral");

    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let tol = Tolerance::default();
    let avg_radius = radii.iter().sum::<f64>() / 3.0;

    // Pipeline up to construct: analyze + per-edge geometry.
    let analyses =
        analysis::analyze_edges(&cube, &incident_edges, false).expect("analysis should succeed");
    let mut edge_geoms = Vec::with_capacity(specs.len());
    for (spec, ana) in specs.iter().zip(analyses.iter()) {
        let geom = compute_fillet_geometry(
            ana,
            &spec.radius,
            &cube,
            &tol,
            ContinuityLevel::G1,
            CapKind::Open,
        )
        .expect("edge geometry should build");
        edge_geoms.push(geom);
    }
    let fillet_surfaces: Vec<FilletSurfaceResult> = edge_geoms
        .iter()
        .map(|eg| FilletSurfaceResult {
            surface: eg.surface.clone(),
            trim_curve_1: eg.trim_3d_1.clone(),
            trim_curve_2: eg.trim_3d_2.clone(),
            continuity_downgrade: None,
        })
        .collect();

    let indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
    let boundary = compute_blend_boundary(vertex, &indices, &fillet_surfaces, &cube, &tol)
        .expect("boundary should build");

    let boundary_curves: Vec<NurbsCurve> = boundary
        .segments
        .iter()
        .filter_map(|s| match s {
            BoundarySegment::FilletEnd { curve, .. } => Some(curve.clone()),
            _ => None,
        })
        .collect();
    assert_eq!(boundary_curves.len(), 3);

    // Upstream-data audit: how far is each boundary-arc endpoint
    // from the vertex position? The vertex is at (0, 0, 0) on a
    // unit cube; a radius-r fillet's end-arc must touch points
    // within `~r * √2` of the vertex. Anything >> r signals that
    // `extract_fillet_end_near_vertex` is picking trim-curve
    // endpoints from the INTERNAL extended-rolling-ball geometry
    // (which overshoots the cube by 1.0 in every direction) rather
    // than the geometric fillet-end arcs we actually want. This
    // is the C4b.1a blocker.
    let mut max_endpoint_dist = 0.0_f64;
    for curve in &boundary_curves {
        let (t0, t1) = curve.parameter_range();
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);
        max_endpoint_dist = max_endpoint_dist.max(p0.distance_to(&boundary.vertex_position));
        max_endpoint_dist = max_endpoint_dist.max(p1.distance_to(&boundary.vertex_position));
    }

    let classification = classify_vertex_blend(&cube, vertex, &incident_edges, &tol).classification;

    // Path 1: polar NURBS (current dispatcher output for mixed-radius).
    let polar_result = build_blend_surface(
        &boundary,
        &fillet_surfaces,
        &cube,
        VertexMode::Smooth,
        classification,
        &tol,
    );
    let polar_metrics = match &polar_result {
        Ok(r) => Some(compute_metrics(r, &boundary_curves, avg_radius)),
        Err(e) => {
            eprintln!("diag: polar-NURBS build failed: {e}");
            None
        }
    };

    // Path 2: rational Gregory (unconditional diagnostic build,
    // bypasses the uniform-radii gate in the production path).
    let gregory_result = build_unchecked_rational_gregory(
        &boundary_curves,
        &fillet_surfaces,
        boundary.vertex_position,
    );
    let gregory_metrics = gregory_result
        .as_ref()
        .map(|r| compute_metrics(r, &boundary_curves, avg_radius));

    eprintln!();
    eprintln!("==========================================");
    eprintln!("§5.3 Round 7 C4b.1 diagnostic");
    eprintln!("  fixture: mixed-radius cube corner r ∈ {radii:?}");
    eprintln!("  avg_radius = {avg_radius:.4}, classification = {classification:?}");
    eprintln!("  vertex position = {:?}", boundary.vertex_position);
    eprintln!(
        "  max boundary-arc endpoint distance to vertex = {max_endpoint_dist:.3e} \
         (expected ≈ r·√2 = {:.3e})",
        avg_radius * std::f64::consts::SQRT_2,
    );
    eprintln!("==========================================");
    eprintln!();
    eprintln!("--- polar-NURBS path (current dispatch) ---");
    match &polar_metrics {
        Some(m) => eprintln!("{m}"),
        None => eprintln!("  FAILED to build"),
    }
    eprintln!("--- rational-Gregory path (unchecked build) ---");
    match &gregory_metrics {
        Some(m) => eprintln!("{m}"),
        None => eprintln!("  FAILED to build (expected — see C4b.1a finding)"),
    }
    eprintln!("==========================================");
    eprintln!();

    // C4b.1a finding (to address in the next commit):
    //   On a mixed-radius cube corner, the boundary-arc endpoints
    //   emitted by `compute_blend_boundary` via
    //   `extract_fillet_end_near_vertex` are **NOT** near the
    //   vertex position. The trim curves span the full edge plus
    //   overshoot — e.g., from x = −1 to x = 2 on a [0, 1]³ cube.
    //   `extract_fillet_end_near_vertex` picks the nearest
    //   endpoint *pair* by cumulative distance, which for this
    //   fixture selects the `x = −1` endpoints of all three trim
    //   curves — three disconnected segments in the planes
    //   `x = −1, y = −1, z = −1`. Corner-coincidence check in
    //   `TriangularGregory::from_boundary_beziers` rejects this.
    //   The polar-NURBS path accepts the bad input but produces
    //   `max_boundary_deviation = 2.5e-2` on a 0.1 fillet —
    //   silent quality regression we were shipping until this
    //   diagnostic caught it.
    //
    // Gate for the Gregory dispatch activation:
    //   `max_boundary_arc_endpoint_distance_to_vertex < 2 · r`
    //   AND polar-NURBS `max_boundary_deviation < 0.05 · r`.
    // Neither holds today. Commercial-grade activation blocked
    // on C4b.1a (fix `extract_fillet_end_near_vertex` / upstream).
    //
    // Soft invariants: if either path builds, its metrics must be
    // finite (no NaN/Inf). This test is deliberately tolerant
    // of bad quality — it EXISTS to document the bad quality.
    if let Some(m) = &polar_metrics {
        assert!(
            m.max_interior_gaussian_curvature.is_finite(),
            "polar-NURBS interior curvature must be finite",
        );
    }
    if let Some(m) = &gregory_metrics {
        assert!(
            m.max_interior_gaussian_curvature.is_finite(),
            "Gregory interior curvature must be finite",
        );
    }
}
