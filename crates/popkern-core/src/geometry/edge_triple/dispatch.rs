//! Public entry points: single-edge, all-branches, batch, preview,
//! curvature comb, revalidation. Stage 1 returns `Unimplemented`
//! from every analytic / Gregory / marching specialisation.

use serde::{Deserialize, Serialize};

use crate::geometry::edge_triple::analytic::try_dispatch_analytic;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, EdgeHint};
use crate::geometry::edge_triple::options::EdgeTripleOptions;
use crate::geometry::edge_triple::types::{CurvatureSample, EdgeTriple, EdgeTripleQuality};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

// ─── Single-edge entry point ──────────────────────────────────────

/// Construct an [`EdgeTriple`] from a surface pair and edge anchors.
///
/// Stage 1: every dispatch path returns
/// [`EdgeTripleError::Unimplemented`]. Subsequent stages route
/// `(surf_a, surf_b)` through analytic / Gregory-boundary /
/// marching specialisations per the ADR §5 dispatch table.
pub fn compute_edge_geometry(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: EdgeAnchors,
    hint: Option<EdgeHint>,
    tol: &Tolerance,
    opts: &EdgeTripleOptions,
) -> Result<EdgeTriple, EdgeTripleError> {
    if let Some(result) = try_dispatch_analytic(surf_a, surf_b, &anchors, hint.as_ref(), tol, opts)
    {
        return result;
    }
    Err(EdgeTripleError::Unimplemented)
}

/// Return every intersection branch as a separate [`EdgeTriple`].
///
/// Used by `preview_edge_geometry` and by UI surfaces that want to
/// show the user every candidate curve before they pick one.
///
/// Stage 1: returns [`EdgeTripleError::Unimplemented`].
pub fn compute_edge_geometry_all_branches(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: EdgeAnchors,
    tol: &Tolerance,
    opts: &EdgeTripleOptions,
) -> Result<Vec<EdgeTriple>, EdgeTripleError> {
    let _ = (surf_a, surf_b, anchors, tol, opts);
    Err(EdgeTripleError::Unimplemented)
}

// ─── Batch ────────────────────────────────────────────────────────

/// One element of a batch request to [`compute_edge_geometry_batch`].
#[derive(Clone, Debug)]
pub struct EdgeRequest<'a> {
    /// Surface on the A side.
    pub surf_a: &'a SurfaceGeometry,
    /// Surface on the B side.
    pub surf_b: &'a SurfaceGeometry,
    /// Edge endpoint / topology descriptor.
    pub anchors: EdgeAnchors,
    /// Optional disambiguation hint.
    pub hint: Option<EdgeHint>,
}

/// Batch-call options.
#[derive(Clone, Debug)]
pub struct EdgeTripleBatchOptions {
    /// Per-request options (cancellation/progress are shared).
    pub per_request: EdgeTripleOptions,
    /// Run in parallel when batch size exceeds this threshold.
    pub parallelism_threshold: usize,
}

impl Default for EdgeTripleBatchOptions {
    fn default() -> Self {
        Self {
            per_request: EdgeTripleOptions::default(),
            parallelism_threshold: 16,
        }
    }
}

/// Batch entry point. Stage 10 wires up rayon parallelism.
///
/// Determinism: results are returned in request order, not
/// completion order, so identical inputs yield bit-identical output
/// regardless of thread scheduling.
///
/// Stage 1: every element is `Err(Unimplemented)`.
pub fn compute_edge_geometry_batch(
    requests: &[EdgeRequest<'_>],
    tol: &Tolerance,
    opts: &EdgeTripleBatchOptions,
) -> Vec<Result<EdgeTriple, EdgeTripleError>> {
    requests
        .iter()
        .map(|r| {
            compute_edge_geometry(
                r.surf_a,
                r.surf_b,
                r.anchors.clone(),
                r.hint.clone(),
                tol,
                &opts.per_request,
            )
        })
        .collect()
}

// ─── Preview ──────────────────────────────────────────────────────

/// Preview: returns the picked branch plus every candidate. UI
/// renders the candidates as ghost curves so the user can pick.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeTriplePreview {
    /// The branch picked by `BranchPolicy` and any caller hint.
    pub triple: EdgeTriple,
    /// Every candidate branch returned by the dispatch.
    pub all_branches: Vec<EdgeTriple>,
    /// Index into `all_branches` of the picked branch.
    pub picked_branch: usize,
}

/// Non-destructive sibling of [`compute_edge_geometry`]. UI calls
/// this on slider drag.
///
/// Stage 1: returns [`EdgeTripleError::Unimplemented`].
pub fn preview_edge_geometry(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: EdgeAnchors,
    hint: Option<EdgeHint>,
    tol: &Tolerance,
    opts: &EdgeTripleOptions,
) -> Result<EdgeTriplePreview, EdgeTripleError> {
    let _ = (surf_a, surf_b, anchors, hint, tol, opts);
    Err(EdgeTripleError::Unimplemented)
}

// ─── Curvature comb + revalidation ────────────────────────────────

/// Sample the 3-D curve's curvature comb in one call. UI renders
/// the result as a ribbon of normals scaled by curvature magnitude.
///
/// Stage 1: returns an empty vec; stage 10 adds the real
/// implementation.
pub fn edge_triple_curvature_comb(triple: &EdgeTriple, n_samples: u32) -> Vec<CurvatureSample> {
    let _ = (triple, n_samples);
    Vec::new()
}

/// Re-check an existing triple against the latest surface pair.
/// UI calls this after a surface edit to decide whether the edge
/// must be recomputed.
pub fn revalidate_edge_triple(
    triple: &EdgeTriple,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tol: &Tolerance,
) -> EdgeTripleQuality {
    triple
        .verify(surf_a, surf_b, tol, 32)
        .unwrap_or_else(|_| EdgeTripleQuality {
            max_positional_residual: f64::INFINITY,
            max_endpoint_error: f64::INFINITY,
            samples: 32,
            tolerance_used: tol.linear(),
        })
}
