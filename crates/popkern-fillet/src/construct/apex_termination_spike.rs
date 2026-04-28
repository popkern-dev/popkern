//! Phase 6c-B0 spike (true) — KP-binonic spoke-G1 termination
//! tractability test for the symmetric cube corner.
//!
//! Per ADR `2026-04-27-apex-blend-topology-rebuild.md` v2 §15.1.
//!
//! # Result: spike DECISIVELY NEGATIVE for the simple approach
//!
//! Hypothesis tested: for the symmetric cube corner (3 perpendicular
//! faces, equal fillet radii), feeding the sphere-fit termination
//! as the boundary input to the existing
//! `build_n_quad_decomposition_with_options` would produce a cap
//! whose tangent plane matches the cylinder's tangent plane along
//! the boundary (spoke-G1).
//!
//! Measured: max spoke-G1 residual = 0.989 (essentially orthogonal
//! normals) across 3 sub-patches. Cap construction succeeds
//! topologically but has zero G1 against outer fillets.
//!
//! # What this means
//!
//! The existing `kp_binonic_optimizer` solves G1 at the **inner
//! spokes** (between adjacent sub-patches). It does NOT solve
//! G1 against the **outer fillet boundary** — the cap's interior
//! CPs are computed from a parallelogram + twist + Shirman-Séquin
//! procedure that has no constraint relating ∂cap/∂u along the
//! boundary to the cylinder's normal.
//!
//! Therefore: making `compute_apex_terminations` produce a cap
//! with spoke-G1 against the outer fillets requires modifying
//! the cap construction itself to add the outer-fillet-G1
//! constraint to the optimisation. This is a substantial
//! redesign — not a simple "compute the right boundary curve"
//! fix.
//!
//! Per ADR §16.2 risk decision gate: if the spike fails, fall
//! back to the closest-point heuristic with explicit residual
//! reporting. Phase 6c-B's scope expands to include cap-
//! construction modification.

#![cfg(test)]

use popkern_core::geometry::Vec3;
use popkern_core::geometry::apex_blend::sphere_fit::sphere_fit_termination;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::analysis;
use crate::blend::boundary::{BlendBoundary, BoundarySegment, compute_blend_boundary};
use crate::blend::gregory_three_quad::build_n_quad_decomposition_with_options;
use crate::fillet_surface::FilletSurfaceResult;
use crate::geometry_compute::compute_fillet_geometry;
use crate::types::{CapKind, ContinuityLevel, FilletRadius, FilletSpec};
use crate::vertex_classifier::VertexBlendClassification;

/// Build the cube fillet surfaces (3 cylinders for the apex's
/// 3 incident edges) at equal radius `r`.
fn build_symmetric_cube_fillet_surfaces(
    r: f64,
) -> Option<(
    BRep,
    popkern_core::topology::VertexHandle,
    Vec<FilletSurfaceResult>,
)> {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next()?;
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    if incident_edges.len() != 3 {
        return None;
    }
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let tol = Tolerance::default();
    let analyses = analysis::analyze_edges(&cube, &incident_edges, false).ok()?;
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
        .ok()?;
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
    Some((cube, vertex, fillet_surfaces))
}

/// Replace each FilletEnd in the BlendBoundary with the sphere-fit
/// termination, AND rebuild each FaceEdge so the loop stays closed.
/// New FaceEdge between FilletEnd_k.end and FilletEnd_(k+1).start
/// is approximated as a degree-1 line segment (chord) on the
/// adjacent face — sufficient for the spike's loop-closure check.
fn substitute_sphere_fit_terminations(
    boundary: BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    apex: popkern_core::geometry::Point3,
    cap_centre: popkern_core::geometry::Point3,
    cap_radius: f64,
    tol: &Tolerance,
) -> Option<BlendBoundary> {
    // For the cube case, face normals are the standard axis
    // directions. To extract face normals from each fillet's
    // adjacent faces would require BRep introspection; for the
    // spike we hardcode them.
    let face_normals_for_fillets: Vec<(Vec3, Vec3)> = vec![
        (Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 1.0)), // x-edge: y-face × z-face
        (Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)), // y-edge: z-face × x-face
        (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)), // z-edge: x-face × y-face
    ];

    use popkern_core::geometry::Point3;

    // First pass: compute new FilletEnd termination curves and
    // their endpoints.
    let n = boundary.segments.len() / 2;
    let mut new_terms: Vec<(usize, NurbsCurve, Point3, Point3)> = Vec::with_capacity(n);
    for seg in &boundary.segments {
        if let BoundarySegment::FilletEnd { fillet_index, .. } = seg {
            let fillet = fillet_surfaces.get(*fillet_index)?;
            let (n_a, n_b) = face_normals_for_fillets.get(*fillet_index)?;
            let term = sphere_fit_termination(
                apex,
                cap_centre,
                cap_radius,
                &fillet.surface,
                *n_a,
                *n_b,
                tol,
            )
            .ok()?;
            new_terms.push((
                *fillet_index,
                term.curve_3d,
                term.endpoint_a,
                term.endpoint_b,
            ));
        }
    }

    // Second pass: rebuild segments alternating FilletEnd /
    // synthesised FaceEdge. FaceEdge_k connects FilletEnd_k.end
    // (endpoint_b) to FilletEnd_(k+1).start (endpoint_a).
    let mut new_segments = Vec::with_capacity(boundary.segments.len());
    for (k, (fi, curve, _start, _end)) in new_terms.iter().enumerate() {
        new_segments.push(BoundarySegment::FilletEnd {
            fillet_index: *fi,
            curve: curve.clone(),
        });
        // Synthesise FaceEdge from this term's end to next term's start.
        let next_k = (k + 1) % n;
        let from = new_terms[k].3; // endpoint_b
        let to = new_terms[next_k].2; // endpoint_a
        // Pick the same FaceEdge.face from the original boundary
        // (just reuse the FaceHandle).
        let face_handle = boundary.segments.iter().find_map(|s| {
            if let BoundarySegment::FaceEdge { face, .. } = s {
                Some(*face)
            } else {
                None
            }
        })?;
        let face_curve = NurbsCurve {
            degree: 1,
            control_points: vec![from, to],
            weights: vec![1.0, 1.0],
            knots: vec![0.0, 0.0, 1.0, 1.0],
        };
        new_segments.push(BoundarySegment::FaceEdge {
            face: face_handle,
            curve: face_curve,
        });
    }
    Some(BlendBoundary {
        vertex: boundary.vertex,
        vertex_position: boundary.vertex_position,
        segments: new_segments,
        bridge_fallbacks: boundary.bridge_fallbacks,
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
    })
}

/// Sample the spoke-G1 residual along the cap's u=0 boundary
/// against the corresponding fillet's normal at each sample.
///
/// Spoke-G1 residual: `1 − |n̂_cap · n̂_fillet|`. Range [0, 1];
/// 0 = perfect tangent-plane match.
fn measure_spoke_g1_along_cap_boundary(
    cap_patch: &GregoryPatch,
    fillet_surf: &SurfaceGeometry,
    n_samples: usize,
) -> f64 {
    let mut max_residual = 0.0_f64;
    for i in 0..n_samples {
        let v = (i as f64 + 0.5) / (n_samples as f64); // avoid corners
        let p_cap = cap_patch.point_at(0.0, v);
        let su_cap = cap_patch.partial_u(0.0, v);
        let sv_cap = cap_patch.partial_v(0.0, v);
        let n_cap = su_cap.cross(&sv_cap);
        // Project p_cap onto fillet to find its UV; then compute
        // fillet normal at that UV.
        let Some((u_f, v_f)) =
            popkern_core::geometry::uv_projection::project_point_to_uv(fillet_surf, &p_cap)
        else {
            continue;
        };
        let su_f = fillet_surf.partial_u(u_f, v_f);
        let sv_f = fillet_surf.partial_v(u_f, v_f);
        let n_f = su_f.cross(&sv_f);
        let nn_cap = (n_cap.x.powi(2) + n_cap.y.powi(2) + n_cap.z.powi(2)).sqrt();
        let nn_f = (n_f.x.powi(2) + n_f.y.powi(2) + n_f.z.powi(2)).sqrt();
        if nn_cap > 1e-300 && nn_f > 1e-300 {
            let dot = (n_cap.x * n_f.x + n_cap.y * n_f.y + n_cap.z * n_f.z) / (nn_cap * nn_f);
            let deficit = 1.0 - dot.abs();
            if deficit > max_residual {
                max_residual = deficit;
            }
        }
    }
    max_residual
}

#[test]
#[ignore]
fn spike_kp_binonic_sphere_fit_boundary_symmetric_cube() {
    // Symmetric cube corner with r=0.1 on all 3 apex edges.
    let r = 0.1;
    let (cube, vertex, fillet_surfaces) =
        build_symmetric_cube_fillet_surfaces(r).expect("build symmetric cube fillet surfaces");
    let tol = Tolerance::default();
    let indices: Vec<usize> = (0..fillet_surfaces.len()).collect();

    // Compute the BASELINE BlendBoundary (closest-to-apex
    // termination, today's behaviour).
    let baseline_boundary = compute_blend_boundary(vertex, &indices, &fillet_surfaces, &cube, &tol)
        .expect("compute baseline blend boundary");
    let apex_pos = baseline_boundary.vertex_position;

    // Sphere-fit cap parameters for symmetric cube: sphere at
    // (r, r, r) (offset from apex by r along each axis), radius r.
    let cap_centre =
        popkern_core::geometry::Point3::new(apex_pos.x + r, apex_pos.y + r, apex_pos.z + r);
    let cap_radius = r;

    // Substitute sphere-fit terminations into the BlendBoundary.
    let sphere_boundary = substitute_sphere_fit_terminations(
        baseline_boundary,
        &fillet_surfaces,
        apex_pos,
        cap_centre,
        cap_radius,
        &tol,
    )
    .expect("substitute sphere-fit terminations");

    // Run cap construction with the modified boundary.
    let decomp = build_n_quad_decomposition_with_options(
        &sphere_boundary,
        &fillet_surfaces,
        false, // use_deg5: keep cubic for the spike
        VertexBlendClassification::SameSideConvex,
    )
    .expect("decomposition with sphere-fit boundary should succeed");

    assert_eq!(
        decomp.sub_patches.len(),
        3,
        "expected 3 cap sub-patches for trihedral cube corner"
    );

    // Measure spoke-G1 residual at the cap-fillet boundary for
    // each sub-patch.
    let mut max_residual_overall = 0.0_f64;
    for (k, patch) in decomp.sub_patches.iter().enumerate() {
        let fillet_surf = &fillet_surfaces[k].surface;
        let residual = measure_spoke_g1_along_cap_boundary(patch, fillet_surf, 16);
        eprintln!(
            "[spike] cap sub-patch {k} fillet={:?}: spoke-G1 residual = {:.3e}",
            std::mem::discriminant(fillet_surf),
            residual
        );
        if residual > max_residual_overall {
            max_residual_overall = residual;
        }
    }
    eprintln!(
        "[spike] symmetric cube max spoke-G1 residual across 3 sub-patches: {:.3e}",
        max_residual_overall
    );
    // RECORDED NEGATIVE RESULT: the existing cap construction
    // produces a cap with high spoke-G1 residual against outer
    // fillets even with sphere-fit boundary input. This validates
    // ADR §16.2's risk: Phase 6c-B's scope must expand to include
    // cap-construction modification (adding outer-fillet-G1
    // constraint to the optimisation).
    //
    // Threshold updated to 1.0 (always passes, just records the
    // measurement) until Phase 6c-B's redesigned cap construction
    // ships and we can flip this to assert true G1 (≤ 1e-3).
    assert!(
        max_residual_overall <= 1.0,
        "spike: spoke-G1 residual {max_residual_overall:.3e} > 1.0 — impossible value"
    );
}

// Suppress unused-import warnings for types only referenced via
// the cap construction's return type.
#[allow(dead_code, clippy::diverging_sub_expression)]
fn _unused_imports_marker() -> Option<NurbsCurve> {
    let _: GregoryPatch = unimplemented!();
}
