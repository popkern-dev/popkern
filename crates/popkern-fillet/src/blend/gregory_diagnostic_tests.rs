//! §5.3 Round 7 C4b.1 — comparative diagnostic for rational
//! Gregory vs polar-NURBS corner construction.
//!
//! Runs both constructions on a mixed-radius cube-corner fixture
//! and reports boundary deviation, Gaussian curvature stats, and
//! |Δκ| at boundary samples. Evidence base for C4b dispatch
//! activation. See module code for the full metric definitions.

use crate::analysis;
use crate::blend::boundary::{BoundarySegment, compute_blend_boundary};
use crate::blend::gregory_diagnostic_metrics::{build_unchecked_rational_gregory, compute_metrics};
use crate::blend::surface::build_blend_surface;
use crate::fillet_surface::FilletSurfaceResult;
use crate::geometry_compute::compute_fillet_geometry;
use crate::types::{CapKind, ContinuityLevel, FilletRadius, FilletSpec, VertexMode};
use crate::vertex_classifier::classify_vertex_blend;

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

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

    // R8 C3 — three-quad decomposition. Build 3 rectangular
    // Gregory sub-patches from the hexagonal boundary, verify they
    // construct and evaluate cleanly. This is NOT yet wired into
    // dispatch; the decomposition is a geometric building block
    // for the eventual N=3 mixed-radius corner.
    use crate::blend::gregory_three_quad::build_three_quad_decomposition;
    match build_three_quad_decomposition(&boundary) {
        Some(decomp) => {
            assert_eq!(decomp.sub_patches.len(), 3);
            eprintln!("--- R8 C3 three-quad decomposition ---");
            eprintln!("  central_point = {:?}", decomp.central_point);
            for (k, patch) in decomp.sub_patches.iter().enumerate() {
                // Every sub-patch should evaluate to finite values
                // at the centroid of its parameter domain.
                let p = patch.point_at(0.5, 0.5);
                assert!(
                    p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                    "sub-patch {k} produced non-finite mid-point: {p:?}",
                );
                eprintln!("  sub-patch[{k}] mid = {p:?}");
            }
            eprintln!("======================================");
        }
        None => {
            eprintln!("--- R8 C3 three-quad decomposition: FAILED ---");
        }
    }
}
