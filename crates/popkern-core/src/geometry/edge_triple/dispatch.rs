//! Public entry points: single-edge, all-branches, batch, preview,
//! curvature comb, revalidation. Stage 1 returns `Unimplemented`
//! from every analytic / Gregory / marching specialisation.

use serde::{Deserialize, Serialize};

use crate::geometry::edge_triple::analytic::{try_dispatch_analytic, try_dispatch_marching};
use crate::geometry::edge_triple::analytic_trig::try_dispatch_trig;
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
    // Stage 2b.39: symbolic short-circuit for coincident surfaces.
    // `A = c·B` for some non-zero rational `c` means the "intersection"
    // is the entire surface — not a 1-D edge. Detected exactly via
    // the Q rational path (no FP tolerance fudge); only fires when
    // both surfaces have rational parameters.
    if let Some(true) = super::levin::surfaces_are_coincident(surf_a, surf_b) {
        return Err(EdgeTripleError::CoincidentSurfaces {
            reason: "surfaces are scalar multiples of one another (A = c · B)",
        });
    }
    // Stage 5: Gregory-boundary dispatch. Runs first when the
    // caller has signalled via `EdgeHint::SurfaceBoundary` that
    // the edge corresponds to one of the Gregory's UV boundaries.
    // Closes the cap-to-fillet seam drift (Task #144) by
    // constructing all three pieces from one source — the
    // Gregory's symbolic boundary curve.
    if let Some(result) = super::gregory_boundary::try_dispatch_gregory_boundary(
        surf_a,
        surf_b,
        &anchors,
        hint.as_ref(),
        tol,
    ) {
        return result;
    }
    if let Some(result) = super::bounded_boundary::try_dispatch_bounded_surface_boundary(
        surf_a,
        surf_b,
        &anchors,
        hint.as_ref(),
        tol,
    ) {
        return result;
    }
    // Try the analytic-trig dispatch first (commercial-style:
    // analytic Circle 3D + Linear cyl pcurve + Elliptical plane
    // pcurve, parameter map = Identity).
    if let Some(result) = try_dispatch_trig(surf_a, surf_b, &anchors, hint.as_ref(), tol) {
        return result;
    }
    // Fall through to the legacy NURBS-via-intersect_surfaces path
    // (handles plane×plane analytically).
    if let Some(result) = try_dispatch_analytic(surf_a, surf_b, &anchors, hint.as_ref(), tol, opts)
    {
        return result;
    }
    // Stage 2b.52: marching catch-all. Runs `intersect_surfaces` on
    // any remaining pair (NURBS×*, torus configurations not in the
    // trig dispatch, etc.) and packages the anchored branch with
    // `Marching` provenance.
    if let Some(result) = try_dispatch_marching(surf_a, surf_b, &anchors, hint.as_ref(), tol, opts)
    {
        return result;
    }
    Err(EdgeTripleError::Unimplemented)
}

/// Return every intersection branch as a separate [`EdgeTriple`].
///
/// Used by `preview_edge_geometry` and by UI surfaces that want to
/// show the user every candidate curve before they pick one. Also
/// the migration target for downstream consumers (`popkern-bool`,
/// `popkern-fillet`) that previously iterated `intersect_surfaces`'s
/// `Vec<SurfaceSurfaceIntersection>`.
///
/// **Multi-branch handling (stage 2b.50):** for **Levin quadric
/// pairs** (cyl×cyl skew, cyl×cone non-coaxial, cone×cone non-coaxial
/// — and the Steinmetz cyl×cyl perpendicular case), the rank-2
/// fibre decomposes into two real planes and this function returns
/// **two** triples, one per plane. Anchors for the plane that doesn't
/// match the user's anchors are synthesised at the **branch points**
/// (intersections of the kernel line with `surf_a`).
///
/// **Single-branch fallback:** for surface pairs without a Levin
/// rank-2 fibre (plane×sphere, sphere×sphere, the various tangent
/// cases, etc.), this function returns a 1-element `Vec` containing
/// the result of [`compute_edge_geometry`].
///
/// **Coverage gaps:**
/// - `cyl × cyl` parallel can have 0/1/2 line branches; the
///   single-branch dispatcher's behaviour is preserved here.
/// - Plane × cylinder oblique with axis-in-plane (2 lines) is
///   not yet enumerated.
///
/// **Error semantics:** If both Levin multi-branch and the
/// single-branch path fail, the single-branch error is returned
/// (e.g., `CoincidentSurfaces`, `NoIntersection`, `Unimplemented`).
/// On success, the returned `Vec` is non-empty.
pub fn compute_edge_geometry_all_branches(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: EdgeAnchors,
    tol: &Tolerance,
    opts: &EdgeTripleOptions,
) -> Result<Vec<EdgeTriple>, EdgeTripleError> {
    // Coincident-surface short-circuit: same exact-symbolic check
    // as `compute_edge_geometry`.
    if let Some(true) = super::levin::surfaces_are_coincident(surf_a, surf_b) {
        return Err(EdgeTripleError::CoincidentSurfaces {
            reason: "surfaces are scalar multiples of one another (A = c · B)",
        });
    }

    // Levin multi-branch enumeration: returns Some(vec) for quadric
    // pairs whose pencil has a rank-≤ 2 degenerate fibre. The vec
    // contains one triple per plane (typically 2 for rank-2
    // mixed-sign).
    if let Some(triples) = super::levin::try_levin_pair_all_branches(surf_a, surf_b, &anchors, tol)
    {
        return Ok(triples);
    }

    // Fall back to single-branch dispatch — pair classes where the
    // intersection has one geometric branch (plane×sphere, etc.) or
    // where Levin doesn't apply.
    let primary = compute_edge_geometry(surf_a, surf_b, anchors, None, tol, opts)?;
    Ok(vec![primary])
}

// ─── Anchorless entry point (stage 2b.58) ────────────────────────

/// Discover every intersection branch between two surfaces and
/// package each as an [`EdgeTriple`], **without** caller-supplied
/// anchors. Use case: face-pair edge discovery in boolean
/// operations / section / imprint, where the curve is being
/// computed from scratch and no endpoints are known *a priori*.
///
/// **API parity:** the UV-range arguments mirror `intersect_surfaces`
/// in `popkern_core::geometry::intersection`. Consumers compute
/// bounded ranges via `uv_projection::bounded_parameter_range`
/// (or similar) from face vertex positions.
///
/// **Implementation:** delegates to `intersect_surfaces` and converts
/// each `CurveOnSurface` branch through `cos_to_triple_unanchored`,
/// using the curve's natural parameter endpoints. Provenance is
/// tagged `Marching` because that's what `intersect_surfaces`
/// reports.
///
/// **Trade-off:** even on analytic pairs (cyl×plane, sphere×plane,
/// etc.), this returns NURBS-fitted curves rather than the analytic
/// primitives the trig dispatcher would produce. Callers that want
/// analytic curves and have anchors should call
/// [`compute_edge_geometry_all_branches`] instead.
///
/// **Coincident-surface short-circuit:** identical to
/// [`compute_edge_geometry`].
#[allow(clippy::too_many_arguments)]
pub fn compute_edge_geometry_unanchored(
    surf_a: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surf_b: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tol: &Tolerance,
    _opts: &EdgeTripleOptions,
) -> Result<Vec<EdgeTriple>, EdgeTripleError> {
    if let Some(true) = super::levin::surfaces_are_coincident(surf_a, surf_b) {
        return Err(EdgeTripleError::CoincidentSurfaces {
            reason: "surfaces are scalar multiples of one another (A = c · B)",
        });
    }
    let result = crate::geometry::intersection::surface_surface::intersect_surfaces(
        surf_a, u1_range, v1_range, surf_b, u2_range, v2_range, tol,
    )
    .map_err(|_| EdgeTripleError::DegenerateSurfacePair {
        reason: "intersect_surfaces returned an internal error",
    })?;
    let mut triples = Vec::new();
    for inter in result.intersections {
        if let crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection::Curve(
            cos,
        ) = inter
        {
            triples.push(super::analytic::cos_to_triple_unanchored(
                &cos, surf_a, surf_b, tol,
            ));
        }
    }
    Ok(triples)
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

#[cfg(test)]
#[path = "dispatch_public_tests.rs"]
mod tests;
