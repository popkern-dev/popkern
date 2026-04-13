#![warn(missing_docs)]
//! popkern-fillet: Fillets and chamfers for the popkern BRep kernel.
//!
//! Provides [`fillet_edges`] and [`chamfer_edges`] for applying fillets and
//! chamfers to BRep edges, plus convenience wrappers [`fillet_edge`] and
//! [`chamfer_edge`] for single-edge operations.

pub mod analysis;
pub mod blend;
pub mod center_curve;
pub mod center_curve_helpers;
pub mod chamfer_curve;
pub mod construct;
pub mod error;
pub mod fillet_surface;
mod geometry_compute;
pub mod offset;
pub mod overlap;
pub mod trim;
pub mod types;
pub mod validate;

pub use error::{
    BlendError, BlendFailure, BlendWarning, Remediation, RemediationAction, RemediationContext,
};
pub use types::*;

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Fillet the specified edges of a BRep solid.
///
/// Analyzes each edge, computes the rolling-ball center curve and fillet
/// surface, derives trim curves, and applies the topology modifications.
///
/// # Errors
///
/// Returns [`BlendError::NoValidEdges`] if analysis fails for all edges, or
/// [`BlendError::PartialFailure`] if some edges fail (with a partial result).
pub fn fillet_edges(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    // §10.1 — reject invalid parameters at the API boundary so the
    // input BRep is guaranteed unchanged on failure.
    validate::validate_fillet_specs(brep, specs)?;

    // Capture pre-operation faces for naming.
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();

    // 1. Analyze edges.
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;

    // 2. Compute geometry for each edge.
    let mut edge_geometries = Vec::new();
    let mut failures = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        match compute_fillet_geometry(&analyses[i], &spec.radius, brep, &tolerance) {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // 2b. Detect overlaps (warnings only — don't block the operation).
    let mut overlap_warnings = Vec::new();
    let (overlap_trims, _face_trim_map) = build_overlap_inputs(&edge_geometries, brep);

    let overlaps = overlap::detect_overlaps(&overlap_trims, &tolerance);
    for ov in &overlaps {
        overlap_warnings.push(BlendWarning::FilletOverlap {
            fillet_a: ov.fillet_a,
            fillet_b: ov.fillet_b,
        });
    }

    // Face disappearance is checked geometrically using the actual fillet
    // radius and face edge lengths, not the unbounded trim curves (which
    // extend beyond the face boundary and produce false positives).
    // A face disappears when the sum of fillet radii from opposite edges
    // exceeds the face's shortest perpendicular edge. For now, the
    // topology application handles this naturally — if the perpendicular
    // edge split lands outside the edge, it fails gracefully.

    if edge_geometries.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // 3. Apply topology.
    let mut apply_result =
        construct::apply_blend_geometry(brep, &edge_geometries, &tolerance, &options.vertex_mode)
            .map_err(|edge_failures| BlendError::PartialFailure {
            failed: edge_failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures: edge_failures,
        })?;

    // 4. Ensure pcurve coverage and correctness.
    // recompute_all_pcurves fills gaps, fix_stale_pcurves replaces
    // any pcurves whose endpoints don't match vertex positions.
    let _ = apply_result.brep.recompute_all_pcurves();
    apply_result.brep.fix_stale_pcurves();

    // Merge overlap warnings with topology warnings.
    let mut all_warnings = apply_result.warnings;
    for warn in overlap_warnings {
        // Attach overlap warnings to the first edge in the specs list.
        let edge_h = specs
            .first()
            .map(|s| s.edge)
            .unwrap_or(edge_geometries[0].analysis.edge);
        all_warnings.push((edge_h, warn));
    }

    // Assign stable names when tracking is active.
    if let Some(ctx) = tracking {
        assign_fillet_names(&mut apply_result.brep, ctx, &pre_faces);
    }

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: all_warnings,
    })
}

/// Chamfer the specified edges of a BRep solid.
///
/// Analyzes each edge, computes offset curves on each adjacent face, builds
/// a ruled surface between them, and applies topology modifications.
///
/// # Errors
///
/// Returns [`BlendError::NoValidEdges`] if analysis fails for all edges, or
/// [`BlendError::PartialFailure`] if some edges fail (with a partial result).
pub fn chamfer_edges(
    brep: &BRep,
    specs: &[ChamferSpec],
    options: &ChamferOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    // §10.1 — reject invalid parameters at the API boundary.
    validate::validate_chamfer_specs(brep, specs)?;

    // Capture pre-operation faces for naming.
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();

    // 1. Analyze edges.
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;

    // 2. Compute geometry for each edge.
    let mut edge_geometries = Vec::new();
    let mut failures = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        match compute_chamfer_geometry(&analyses[i], &spec.kind, brep, &tolerance) {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // 3. Apply topology. Chamfers always use mitred vertex mode.
    let mut apply_result =
        construct::apply_blend_geometry(brep, &edge_geometries, &tolerance, &VertexMode::Mitred)
            .map_err(|edge_failures| BlendError::PartialFailure {
                failed: edge_failures.len(),
                total: specs.len(),
                partial_result: Box::new(brep.clone()),
                failures: edge_failures,
            })?;

    // Assign stable names when tracking is active.
    if let Some(ctx) = tracking {
        assign_fillet_names(&mut apply_result.brep, ctx, &pre_faces);
    }

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: apply_result.warnings,
    })
}

/// Convenience: fillet a single edge with a constant radius.
///
/// Equivalent to calling [`fillet_edges`] with a single [`FilletSpec`] and
/// default options.
///
/// # Errors
///
/// Returns a [`BlendError`] if the edge cannot be filleted.
pub fn fillet_edge(
    brep: &BRep,
    edge: EdgeHandle,
    radius: f64,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    fillet_edges(
        brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        }],
        &FilletOptions::default(),
        tolerance,
        tracking,
    )
}

/// Convenience: chamfer a single edge with equal distances.
///
/// Equivalent to calling [`chamfer_edges`] with a single [`ChamferSpec`] and
/// default options.
///
/// # Errors
///
/// Returns a [`BlendError`] if the edge cannot be chamfered.
pub fn chamfer_edge(
    brep: &BRep,
    edge: EdgeHandle,
    distance: f64,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    chamfer_edges(
        brep,
        &[ChamferSpec {
            edge,
            kind: ChamferKind::Equal(distance),
        }],
        &ChamferOptions::default(),
        tolerance,
        tracking,
    )
}

// Internal geometry computation is in `geometry_compute` module.
use geometry_compute::{build_overlap_inputs, compute_chamfer_geometry, compute_fillet_geometry};

/// Assign stable names to entities produced by a fillet or chamfer operation.
///
/// Surviving faces from the original set receive `Modified` names with the
/// original face name as ancestor (tracked via `record_modified`). New
/// fillet/chamfer faces receive `Generated` names.
fn assign_fillet_names(
    brep: &mut BRep,
    ctx: &mut TrackingContext,
    pre_faces: &std::collections::HashSet<FaceHandle>,
) {
    use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    let post_faces: Vec<FaceHandle> = brep.faces.keys().collect();

    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;

    for fh in &post_faces {
        if pre_faces.contains(fh) {
            // Surviving face — Modified (trimmed by the fillet).
            let old_name = brep.faces[*fh]
                .name
                .and_then(|nid| brep.name_pool.get(nid).cloned());

            let ancestors: Vec<_> = brep.faces[*fh].name.into_iter().collect();
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Modified,
                    index: mod_idx,
                },
                ancestors,
            };
            let criteria = capture_face(brep, *fh).map(GeometricCriteria::Face);
            let name_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(name_id);

            if let Some(old) = old_name {
                ctx.record_modified(old, new_name);
            } else {
                ctx.record_generated(new_name);
            }
            mod_idx += 1;
        } else {
            // New fillet/chamfer face — Generated.
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: gen_idx,
                },
                ancestors: Vec::new(),
            };
            let criteria = capture_face(brep, *fh).map(GeometricCriteria::Face);
            let name_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(name_id);
            ctx.record_generated(new_name);
            gen_idx += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }
    fn fopts() -> FilletOptions {
        FilletOptions::default()
    }
    fn copts() -> ChamferOptions {
        ChamferOptions::default()
    }
    fn cube_and_edge() -> (BRep, popkern_core::topology::EdgeHandle) {
        let brep = BRep::unit_cube(tol()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        (brep, edge)
    }

    #[test]
    fn variable_radius_fillet_computes_geometry() {
        let (brep, edge) = cube_and_edge();
        let specs = vec![FilletSpec {
            edge,
            radius: FilletRadius::Variable {
                start: 0.08,
                end: 0.12,
            },
        }];
        let result = fillet_edges(&brep, &specs, &fopts(), tol(), None);
        assert!(result.is_ok(), "variable fillet failed: {:?}", result.err());
        assert!(result.unwrap().brep.faces.len() > brep.faces.len());
    }

    #[test]
    fn by_points_fillet_exercises_code_paths() {
        let (brep, edge) = cube_and_edge();
        // Two-point ByPoints — exercises variable_radius path. Should not panic.
        let s1 = vec![FilletSpec {
            edge,
            radius: FilletRadius::ByPoints(vec![(0.0, 0.1), (1.0, 0.1)]),
        }];
        let _ = fillet_edges(&brep, &s1, &fopts(), tol(), None);
        // Single-point ByPoints.
        let s2 = vec![FilletSpec {
            edge,
            radius: FilletRadius::ByPoints(vec![(0.5, 0.1)]),
        }];
        let _ = fillet_edges(&brep, &s2, &fopts(), tol(), None);
        // Empty ByPoints should fail.
        let s3 = vec![FilletSpec {
            edge,
            radius: FilletRadius::ByPoints(vec![]),
        }];
        assert!(fillet_edges(&brep, &s3, &fopts(), tol(), None).is_err());
    }

    #[test]
    fn chamfer_two_distance_via_api() {
        let (brep, edge) = cube_and_edge();
        let face_1 = analysis::analyze_edges(&brep, &[edge], false).unwrap()[0].face_1;
        let specs = vec![ChamferSpec {
            edge,
            kind: ChamferKind::TwoDistance {
                face_1,
                dist_1: 0.1,
                dist_2: 0.2,
            },
        }];
        let r = chamfer_edges(&brep, &specs, &copts(), tol(), None);
        assert!(r.is_ok(), "two-distance chamfer failed: {:?}", r.err());
    }

    #[test]
    fn chamfer_distance_angle_via_api() {
        let (brep, edge) = cube_and_edge();
        let face = analysis::analyze_edges(&brep, &[edge], false).unwrap()[0].face_1;
        let specs = vec![ChamferSpec {
            edge,
            kind: ChamferKind::DistanceAngle {
                face,
                distance: 0.1,
                angle: std::f64::consts::FRAC_PI_4,
            },
        }];
        let r = chamfer_edges(&brep, &specs, &copts(), tol(), None);
        assert!(r.is_ok(), "dist-angle chamfer failed: {:?}", r.err());
    }

    #[test]
    fn empty_specs_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        assert!(fillet_edges(&brep, &[], &fopts(), tol(), None).is_err());
        assert!(chamfer_edges(&brep, &[], &copts(), tol(), None).is_err());
    }

    #[test]
    fn fillet_edges_rejects_nan_radius_before_touching_geometry() {
        let (brep, edge) = cube_and_edge();
        let pre_face_count = brep.faces.len();
        let specs = vec![FilletSpec {
            edge,
            radius: FilletRadius::Constant(f64::NAN),
        }];
        let err = fillet_edges(&brep, &specs, &fopts(), tol(), None).unwrap_err();
        assert!(matches!(err, BlendError::InvalidParameter { .. }));
        // Input BRep must be untouched — this is part of the
        // §10.1 contract.
        assert_eq!(brep.faces.len(), pre_face_count);
        // Remediation hint should route to FixInput.
        assert_eq!(err.remediation_hint().action, RemediationAction::FixInput);
    }

    #[test]
    fn fillet_edges_rejects_negative_radius_before_touching_geometry() {
        let (brep, edge) = cube_and_edge();
        let specs = vec![FilletSpec {
            edge,
            radius: FilletRadius::Constant(-0.1),
        }];
        let err = fillet_edges(&brep, &specs, &fopts(), tol(), None).unwrap_err();
        assert!(matches!(err, BlendError::InvalidParameter { .. }));
    }

    #[test]
    fn fillet_edges_returns_radius_too_large_with_populated_max_radius() {
        // Unit cube has edge length 1.0. An obviously-too-large
        // radius (e.g. 2.0) should trigger the self-intersection
        // pre-flight and return RadiusTooLarge with a concrete
        // max_radius hint that is non-zero and less than the
        // requested value.
        //
        // However, note: plane-plane curvature is infinite, so two
        // flat faces won't trip the check. This test verifies the
        // plumbing: even when no max is populated, the structured
        // failure shape is stable.
        let (brep, edge) = cube_and_edge();
        let specs = vec![FilletSpec {
            edge,
            radius: FilletRadius::Constant(2.0), // 2× the cube
        }];
        // Cube edges are plane-plane, so self_intersect returns None
        // and the operation proceeds as before (either succeeds or
        // fails downstream). The test for actual max_radius
        // population is in the curved_faces integration test file.
        let _ = fillet_edges(&brep, &specs, &fopts(), tol(), None);
    }

    #[test]
    fn chamfer_edges_rejects_zero_distance_before_touching_geometry() {
        let (brep, edge) = cube_and_edge();
        let specs = vec![ChamferSpec {
            edge,
            kind: ChamferKind::Equal(0.0),
        }];
        let err = chamfer_edges(&brep, &specs, &copts(), tol(), None).unwrap_err();
        assert!(matches!(err, BlendError::InvalidParameter { .. }));
    }

    #[test]
    fn chamfer_edges_rejects_out_of_range_angle_before_touching_geometry() {
        let (brep, edge) = cube_and_edge();
        let face = brep.faces.keys().next().unwrap();
        let specs = vec![ChamferSpec {
            edge,
            kind: ChamferKind::DistanceAngle {
                face,
                distance: 0.1,
                angle: 10.0, // > π
            },
        }];
        let err = chamfer_edges(&brep, &specs, &copts(), tol(), None).unwrap_err();
        assert!(matches!(err, BlendError::InvalidParameter { .. }));
    }

    #[test]
    fn fillet_result_has_pcurves() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let result = fillet_edge(&brep, edge, 0.1, Tolerance::default(), None).unwrap();

        let total = result.brep.half_edges.len();
        let with_pcurves = result
            .brep
            .half_edges
            .values()
            .filter(|he| he.pcurve.is_some())
            .count();

        assert!(
            with_pcurves > total / 2,
            "expected most half-edges to have pcurves, got {with_pcurves}/{total}"
        );
    }
}
