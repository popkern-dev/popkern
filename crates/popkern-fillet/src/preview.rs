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
use popkern_core::topology::BRep;

use crate::analysis;
use crate::construct::EdgeGeometry;
use crate::error::BlendError;
use crate::geometry_compute::{compute_chamfer_geometry, compute_fillet_geometry};
use crate::types::{ChamferSpec, FilletSpec};
use crate::validate;

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
        // Preview always uses the legacy G1 cross-section (g2 = false)
        // — preview is for cheap UI feedback and the G2 path is for
        // commit-quality output.
        match compute_fillet_geometry(&analyses[i], &spec.radius, brep, &tolerance, false) {
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

#[cfg(test)]
#[path = "preview_tests.rs"]
mod tests;
