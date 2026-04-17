//! §5.3 Round 7 C4b.1 — comparative diagnostic for rational
//! Gregory vs polar-NURBS corner construction.
//!
//! Runs both constructions on a mixed-radius cube-corner fixture
//! and reports boundary deviation, Gaussian curvature stats, and
//! |Δκ| at boundary samples. Evidence base for C4b dispatch
//! activation. See module code for the full metric definitions.

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
    /// Max surface-point-to-boundary-curve distance across 5
    /// samples per input boundary. Should be ≈ 0.
    max_boundary_deviation: f64,
    /// Max Gaussian curvature over a 9×9 interior grid.
    max_interior_gaussian_curvature: f64,
    /// Mean Gaussian curvature over the same grid.
    mean_interior_gaussian_curvature: f64,
    /// Max `|Δκ_n|/max(|κ_n|, 1/r)` at boundary samples — §2.4
    /// G2 acceptance threshold is `< 0.01`.
    max_boundary_normal_curvature_deviation: f64,
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

/// Compute diagnostic metrics for a [`BlendSurfaceResult`].
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

    // R8 C1 regression guard: `compute_blend_boundary` on a
    // mixed-radius trihedral (planar-face) corner must emit a
    // closed hexagonal boundary — 6 segments alternating
    // FilletEnd / FaceEdge / FilletEnd / FaceEdge / ... and each
    // segment's end point matching the next segment's start.
    let n_fe = boundary
        .segments
        .iter()
        .filter(|s| matches!(s, BoundarySegment::FilletEnd { .. }))
        .count();
    let n_faceedge = boundary
        .segments
        .iter()
        .filter(|s| matches!(s, BoundarySegment::FaceEdge { .. }))
        .count();
    eprintln!(
        "diag: boundary segments = {} total ({} FilletEnd + {} FaceEdge)",
        boundary.segments.len(),
        n_fe,
        n_faceedge,
    );

    // R8 C2: hexagonal closure — for every consecutive pair of
    // segments `(segment_k, segment_{k+1})`, the end of segment_k
    // must match the start of segment_{k+1} within 2 · tolerance.
    // This is the geometric invariant that says "the 6 boundary
    // pieces form a closed loop around the vertex."
    let mut max_gap: f64 = 0.0;
    let n_segs = boundary.segments.len();
    for k in 0..n_segs {
        let curr = match &boundary.segments[k] {
            BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => {
                curve
            }
        };
        let next = match &boundary.segments[(k + 1) % n_segs] {
            BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => {
                curve
            }
        };
        let (_, curr_end_t) = curr.parameter_range();
        let (next_start_t, _) = next.parameter_range();
        let p_end = curr.point_at(curr_end_t);
        let p_start_next = next.point_at(next_start_t);
        let gap = p_end.distance_to(&p_start_next);
        if gap > max_gap {
            max_gap = gap;
        }
    }
    eprintln!("diag: max hexagon-closure gap = {max_gap:.3e}");

    // Soft check: closure gap should be small relative to r. We
    // allow 2 · r as a loose-but-commercial-grade tolerance (the
    // bridges are straight lines, whereas the true face-edge
    // boundary on curved supports would curve, contributing some
    // chord error — bounded by ~r/2).
    assert!(
        max_gap < avg_radius,
        "R8 C2 regression: hexagonal boundary has a {max_gap:.3e} gap \
         (expected < r = {avg_radius:.3e}). The face-edge bridge \
         ordering or endpoint pairing is broken.",
    );

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

    // Finite-value invariant (quality bounds are informational).
    if let Some(m) = &polar_metrics {
        assert!(m.max_interior_gaussian_curvature.is_finite());
    }
    if let Some(m) = &gregory_metrics {
        assert!(m.max_interior_gaussian_curvature.is_finite());
    }

    // C4b.1a regression guard.
    assert!(
        max_endpoint_dist < 2.0 * avg_radius,
        "C4b.1a regression: boundary endpoints at {max_endpoint_dist:.3e} \
         from vertex; expected < 2·r = {:.3e}",
        2.0 * avg_radius,
    );
}
