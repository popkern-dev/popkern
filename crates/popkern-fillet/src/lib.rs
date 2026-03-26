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

pub use error::{BlendError, BlendFailure, BlendWarning};
pub use types::*;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

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
) -> Result<BlendResult, BlendError> {
    if specs.is_empty() {
        return Err(BlendError::NoValidEdges("no fillet specs provided".into()));
    }

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

    // 2b. Detect overlaps and face disappearance before topology construction.
    let mut overlap_warnings = Vec::new();
    let (overlap_trims, face_trim_map) = build_overlap_inputs(&edge_geometries, brep);

    let overlaps = overlap::detect_overlaps(&overlap_trims, &tolerance);
    for ov in &overlaps {
        overlap_warnings.push(BlendWarning::FilletOverlap {
            fillet_a: ov.fillet_a,
            fillet_b: ov.fillet_b,
        });
    }

    let face_trims_vec: Vec<(usize, Vec<popkern_core::geometry::NurbsCurve>)> =
        face_trim_map.into_iter().collect();
    let consumed_faces = overlap::detect_face_disappearance(&face_trims_vec, &tolerance);

    // If a face would disappear, report RadiusTooLarge for the responsible edges.
    if !consumed_faces.is_empty() {
        let consumed_face_indices: std::collections::HashSet<usize> =
            consumed_faces.into_iter().collect();
        for (i, eg) in edge_geometries.iter().enumerate() {
            let f1_idx = face_handle_to_index(brep, eg.analysis.face_1);
            let f2_idx = face_handle_to_index(brep, eg.analysis.face_2);
            if consumed_face_indices.contains(&f1_idx) || consumed_face_indices.contains(&f2_idx) {
                failures.push((
                    specs.get(i).map(|s| s.edge).unwrap_or(eg.analysis.edge),
                    BlendFailure::RadiusTooLarge { max_radius: 0.0 },
                ));
            }
        }
    }

    // If all edges now fail (original failures + face-disappearance), bail.
    if edge_geometries.is_empty() || (!failures.is_empty() && failures.len() >= specs.len()) {
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

    // 4. Ensure pcurve coverage. The custom attach_fillet_pcurves handles
    // fillet-specific edges. recompute_all_pcurves fills gaps (half-edges
    // without pcurves). Then replace any stale pcurves whose endpoints
    // no longer match the actual edge curve positions.
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
) -> Result<BlendResult, BlendError> {
    if specs.is_empty() {
        return Err(BlendError::NoValidEdges("no chamfer specs provided".into()));
    }

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
    let apply_result =
        construct::apply_blend_geometry(brep, &edge_geometries, &tolerance, &VertexMode::Mitred)
            .map_err(|edge_failures| BlendError::PartialFailure {
                failed: edge_failures.len(),
                total: specs.len(),
                partial_result: Box::new(brep.clone()),
                failures: edge_failures,
            })?;

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
) -> Result<BlendResult, BlendError> {
    fillet_edges(
        brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        }],
        &FilletOptions::default(),
        tolerance,
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
) -> Result<BlendResult, BlendError> {
    chamfer_edges(
        brep,
        &[ChamferSpec {
            edge,
            kind: ChamferKind::Equal(distance),
        }],
        &ChamferOptions::default(),
        tolerance,
    )
}

// Internal geometry computation is in `geometry_compute` module.
use geometry_compute::{
    build_overlap_inputs, compute_chamfer_geometry, compute_fillet_geometry, face_handle_to_index,
};

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
        let result = fillet_edges(&brep, &specs, &fopts(), tol());
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
        let _ = fillet_edges(&brep, &s1, &fopts(), tol());
        // Single-point ByPoints.
        let s2 = vec![FilletSpec {
            edge,
            radius: FilletRadius::ByPoints(vec![(0.5, 0.1)]),
        }];
        let _ = fillet_edges(&brep, &s2, &fopts(), tol());
        // Empty ByPoints should fail.
        let s3 = vec![FilletSpec {
            edge,
            radius: FilletRadius::ByPoints(vec![]),
        }];
        assert!(fillet_edges(&brep, &s3, &fopts(), tol()).is_err());
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
        let r = chamfer_edges(&brep, &specs, &copts(), tol());
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
        let r = chamfer_edges(&brep, &specs, &copts(), tol());
        assert!(r.is_ok(), "dist-angle chamfer failed: {:?}", r.err());
    }

    #[test]
    fn empty_specs_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        assert!(fillet_edges(&brep, &[], &fopts(), tol()).is_err());
        assert!(chamfer_edges(&brep, &[], &copts(), tol()).is_err());
    }

    #[test]
    fn fillet_result_has_pcurves() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let result = fillet_edge(&brep, edge, 0.1, Tolerance::default()).unwrap();

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
