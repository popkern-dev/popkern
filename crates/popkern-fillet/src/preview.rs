//! Non-committing preview of fillet geometry (§14.6, §17).
//!
//! Preview mode runs the fillet pipeline through the geometry-
//! computation phase but **stops before the topology-modification
//! phase**. The caller gets back the computed fillet / chamfer
//! surfaces and their trim curves — enough to render a
//! highlight-preview in a UI — without paying for the full
//! topology surgery and without mutating the input BRep.
//!
//! Typical UI flow:
//!
//! 1. User selects an edge and drags a radius slider.
//! 2. On each slider move, the UI calls [`preview_fillet`] and
//!    draws the returned [`PreviewGeometry::blend_surfaces`] as
//!    wireframe or lightweight tessellation.
//! 3. On commit, the UI calls the full [`crate::fillet_edges`] to
//!    land the result into the model.
//!
//! The preview path uses exactly the same validation and
//! geometry-computation routines as the real pipeline, so preview
//! results are guaranteed to match the committed result modulo the
//! topology-integration stage.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

use crate::analysis;
use crate::construct::EdgeGeometry;
use crate::error::BlendError;
use crate::geometry_compute::{compute_chamfer_geometry, compute_fillet_geometry};
use crate::types::{ChamferSpec, ContinuityLevel, CornerMethod, FilletSpec};
use crate::validate;
use crate::vertex_classifier::VertexBlendClassification;

/// Per-edge preview geometry — the surfaces and trim curves the
/// full pipeline would assemble into a finished blend, without the
/// topology-integration step.
#[derive(Debug, Clone)]
pub struct PreviewEdgeGeometry {
    /// The edge handle the preview was computed for.
    pub edge: popkern_core::topology::EdgeHandle,
    /// The blend surface (rolling-ball torus section, NURBS sweep,
    /// or chamfer ruled surface).
    pub blend_surface: SurfaceGeometry,
    /// 3D trim curve on the first support face.
    pub trim_curve_1: NurbsCurve,
    /// 3D trim curve on the second support face.
    pub trim_curve_2: NurbsCurve,
    /// The effective radius (or chamfer distance) the preview used.
    pub effective_radius: f64,
}

/// Preview result for a batch of fillet or chamfer specs.
///
/// `per_edge` has one entry per spec in the input list. `failed`
/// reports specs that could not produce preview geometry (e.g.
/// a requested radius that exceeds the self-intersection bound).
#[derive(Debug, Clone)]
pub struct PreviewGeometry {
    /// Preview geometry per successfully-previewed edge.
    pub per_edge: Vec<PreviewEdgeGeometry>,
    /// Edges for which the preview could not be computed, keyed by
    /// spec index in the input list.
    pub failed: Vec<(usize, crate::error::BlendFailure)>,
}

impl PreviewGeometry {
    /// Returns `true` when every spec produced preview geometry.
    pub fn all_succeeded(&self) -> bool {
        self.failed.is_empty()
    }
}

/// Non-destructive preview of fillet geometry.
///
/// Runs parameter validation, edge analysis, and per-edge geometry
/// computation — but **stops before the topology modification
/// phase**. Returns the computed blend surfaces and trim curves
/// without mutating the input BRep.
///
/// The function is a pure query over `brep`. The returned
/// [`PreviewGeometry`] references no part of the BRep and can
/// outlive it.
///
/// # Errors
///
/// Returns the same [`BlendError`] variants as the main pipeline
/// for validation and analysis failures. Per-edge geometry failures
/// are *not* errors at the top level — they go into
/// [`PreviewGeometry::failed`] so the UI can render the remaining
/// successful previews alongside per-edge diagnostics.
pub fn preview_fillet(
    brep: &BRep,
    specs: &[FilletSpec],
    tolerance: Tolerance,
) -> Result<PreviewGeometry, BlendError> {
    validate::validate_fillet_specs(brep, specs)?;

    let edges: Vec<_> = specs.iter().map(|s| s.edge).collect();
    let analyses =
        analysis::analyze_edges(brep, &edges, false).map_err(BlendError::NoValidEdges)?;

    let mut per_edge = Vec::with_capacity(specs.len());
    let mut failed = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        // Preview always uses G1 — preview is for cheap UI feedback
        // and the G2 path is for commit-quality output.
        match compute_fillet_geometry(
            &analyses[i],
            &spec.radius,
            brep,
            &tolerance,
            crate::types::ContinuityLevel::G1,
            crate::types::CapKind::Open,
        ) {
            Ok(geom) => per_edge.push(edge_geometry_to_preview(spec.edge, geom)),
            Err(f) => failed.push((i, f)),
        }
    }

    Ok(PreviewGeometry { per_edge, failed })
}

/// Non-destructive preview of chamfer geometry.
pub fn preview_chamfer(
    brep: &BRep,
    specs: &[ChamferSpec],
    tolerance: Tolerance,
) -> Result<PreviewGeometry, BlendError> {
    validate::validate_chamfer_specs(brep, specs)?;

    let edges: Vec<_> = specs.iter().map(|s| s.edge).collect();
    let analyses =
        analysis::analyze_edges(brep, &edges, false).map_err(BlendError::NoValidEdges)?;

    let mut per_edge = Vec::with_capacity(specs.len());
    let mut failed = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        match compute_chamfer_geometry(&analyses[i], &spec.kind, brep, &tolerance) {
            Ok(geom) => per_edge.push(edge_geometry_to_preview(spec.edge, geom)),
            Err(f) => failed.push((i, f)),
        }
    }

    Ok(PreviewGeometry { per_edge, failed })
}

fn edge_geometry_to_preview(
    edge: popkern_core::topology::EdgeHandle,
    geom: EdgeGeometry,
) -> PreviewEdgeGeometry {
    PreviewEdgeGeometry {
        edge,
        blend_surface: geom.surface,
        trim_curve_1: geom.trim_3d_1,
        trim_curve_2: geom.trim_3d_2,
        effective_radius: geom.effective_radius,
    }
}

/// §5.3 C7 — per-vertex preview of the corner patch the pipeline
/// **would** produce for a given fillet spec list, without running
/// topology surgery.
///
/// Built alongside [`preview_fillet`] so a UI can show corner
/// annotations (method badge, classification, predicted continuity)
/// on a radius-slider drag before the user commits the operation.
#[derive(Debug, Clone)]
pub struct CornerPatchPreview {
    /// The vertex the corner would be resolved at.
    pub vertex: VertexHandle,
    /// Number of incident filleted edges.
    pub n_sides: usize,
    /// §5.3 vertex-blend classification that would drive method
    /// dispatch.
    pub classification: VertexBlendClassification,
    /// Predicted construction method.
    pub predicted_method: CornerMethod,
    /// Predicted continuity with the incident fillet strips.
    pub predicted_continuity: ContinuityLevel,
    /// The corner surface geometry the pipeline would assign.
    pub surface: SurfaceGeometry,
    /// Boundary curves tracing the corner perimeter.
    pub boundary_curves: Vec<NurbsCurve>,
}

/// Non-destructive preview of the N-sided corner patches that would
/// be inserted at each multi-fillet vertex if [`crate::fillet_edges`]
/// were called with the same specs.
///
/// For each vertex touched by two or more specs in `specs`, runs
/// the §5.3 classifier and the corner blend-surface construction
/// (same path as the committed pipeline) but **does not** mutate
/// topology.
///
/// # Errors
///
/// Returns the same validation / analysis errors as
/// [`preview_fillet`]. Per-vertex failures (e.g. degenerate
/// corner geometry) are silently skipped — the return vec holds
/// one entry per successful preview.
pub fn preview_corner_patches(
    brep: &BRep,
    specs: &[FilletSpec],
    tolerance: Tolerance,
) -> Result<Vec<CornerPatchPreview>, BlendError> {
    validate::validate_fillet_specs(brep, specs)?;

    let edges: Vec<_> = specs.iter().map(|s| s.edge).collect();
    let analyses =
        analysis::analyze_edges(brep, &edges, false).map_err(BlendError::NoValidEdges)?;

    let mut edge_geometries: Vec<EdgeGeometry> = Vec::with_capacity(analyses.len());
    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        if let Ok(geom) = compute_fillet_geometry(
            &analyses[i],
            &spec.radius,
            brep,
            &tolerance,
            ContinuityLevel::G1,
            crate::types::CapKind::Open,
        ) {
            edge_geometries.push(geom);
        }
    }

    let mut vertex_edges: std::collections::HashMap<VertexHandle, Vec<usize>> =
        std::collections::HashMap::new();
    for (idx, eg) in edge_geometries.iter().enumerate() {
        vertex_edges
            .entry(eg.analysis.vertex_start)
            .or_default()
            .push(idx);
        vertex_edges
            .entry(eg.analysis.vertex_end)
            .or_default()
            .push(idx);
    }

    let mut out: Vec<CornerPatchPreview> = Vec::new();
    for (vertex, edge_indices) in &vertex_edges {
        if edge_indices.len() < 2 {
            continue;
        }
        let filleted: Vec<_> = edge_indices
            .iter()
            .filter_map(|&i| edge_geometries.get(i).map(|g| g.analysis.edge))
            .collect();
        let classification =
            crate::vertex_classifier::classify_vertex_blend(brep, *vertex, &filleted, &tolerance)
                .classification;

        let fillet_surfaces: Vec<_> = edge_geometries
            .iter()
            .enumerate()
            .filter(|(i, _)| edge_indices.contains(i))
            .map(|(_, eg)| crate::fillet_surface::FilletSurfaceResult {
                surface: eg.surface.clone(),
                trim_curve_1: eg.trim_3d_1.clone(),
                trim_curve_2: eg.trim_3d_2.clone(),
                continuity_downgrade: None,
            })
            .collect();
        let local_indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
        let Ok(boundary) = crate::blend::compute_blend_boundary(
            *vertex,
            &local_indices,
            &fillet_surfaces,
            brep,
            &tolerance,
        ) else {
            continue;
        };
        let Ok(result) = crate::blend::build_blend_surface(
            &boundary,
            &fillet_surfaces,
            brep,
            crate::types::VertexMode::Smooth,
            classification,
            &tolerance,
        ) else {
            continue;
        };

        out.push(CornerPatchPreview {
            vertex: *vertex,
            n_sides: edge_indices.len(),
            classification,
            predicted_method: result.method,
            predicted_continuity: result.continuity,
            surface: result.surface,
            boundary_curves: result.boundary_curves,
        });
    }
    Ok(out)
}

#[cfg(test)]
#[path = "preview_tests.rs"]
mod tests;
