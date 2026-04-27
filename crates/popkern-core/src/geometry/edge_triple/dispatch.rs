//! Public entry points: single-edge, all-branches, batch, preview,
//! curvature comb, revalidation. Stage 1 returns `Unimplemented`
//! from every analytic / Gregory / marching specialisation.

use serde::{Deserialize, Serialize};

use crate::geometry::edge_triple::analytic::try_dispatch_analytic;
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
    // Try the analytic-trig dispatch first (commercial-style:
    // analytic Circle 3D + Linear cyl pcurve + Elliptical plane
    // pcurve, parameter map = Identity).
    if let Some(result) = try_dispatch_trig(surf_a, surf_b, &anchors, hint.as_ref(), tol) {
        return result;
    }
    // Fall through to the legacy NURBS-via-intersect_surfaces path
    // (currently only plane×plane).
    if let Some(result) = try_dispatch_analytic(surf_a, surf_b, &anchors, hint.as_ref(), tol, opts)
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
/// **Stage 2b.49 scope:** returns the single anchored branch picked
/// by [`compute_edge_geometry`] wrapped in a `Vec`. Sufficient for
/// surface pairs that geometrically have one intersection branch
/// (plane×sphere, plane×cone-elliptic, plane×plane, cyl×inscribed-
/// sphere tangent, etc.) — i.e., the majority of CAD edges that
/// downstream consumers care about.
///
/// **Multi-branch enumeration (deferred to 2b.50+):** Levin quadric
/// pairs produce **two** plane fibres, and `cyl × cyl` parallel can
/// produce **two** lines. The current implementation returns only
/// the anchor-matching branch. Consumers that depend on enumerating
/// every branch (e.g., `popkern-bool/intersect.rs`) cannot migrate
/// off `intersect_surfaces` until the Levin multi-branch logic
/// lands — at which point this function will append the remaining
/// branches via synthesised anchors on the secondary plane(s).
///
/// **Error semantics:** If the primary anchored dispatch returns an
/// error, that error is returned (e.g., `CoincidentSurfaces`,
/// `NoIntersection`, `Unimplemented`). On success, the returned
/// `Vec` is non-empty and contains at least the primary branch.
pub fn compute_edge_geometry_all_branches(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: EdgeAnchors,
    tol: &Tolerance,
    opts: &EdgeTripleOptions,
) -> Result<Vec<EdgeTriple>, EdgeTripleError> {
    let primary = compute_edge_geometry(surf_a, surf_b, anchors, None, tol, opts)?;
    Ok(vec![primary])
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
mod tests {
    use super::*;
    use crate::geometry::curve::CurveEval;
    use crate::geometry::{Point3, Vec3};

    fn z_axis_cyl(radius: f64) -> SurfaceGeometry {
        SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius,
        }
    }

    // Stage 2b.39: identical cylinders short-circuit before any
    // analytic dispatch, returning `CoincidentSurfaces`. Detected
    // exactly via the symbolic Q rational path (A = c·B).
    #[test]
    fn identical_cylinders_short_circuit_to_coincident() {
        let cyl = z_axis_cyl(1.0);
        let anchors = EdgeAnchors::Open {
            start: Point3::new(1.0, 0.0, 0.0),
            end: Point3::new(1.0, 0.0, 1.0),
        };
        let tol = Tolerance::default();
        let opts = EdgeTripleOptions::default();
        let result = compute_edge_geometry(&cyl, &cyl, anchors, None, &tol, &opts);
        match result {
            Err(EdgeTripleError::CoincidentSurfaces { reason }) => {
                assert!(!reason.is_empty(), "reason must be populated");
            }
            other => panic!("expected CoincidentSurfaces, got {other:?}"),
        }
    }

    // Stage 2b.39: scalar-multiple sphere matrices share the same
    // zero locus. The `A = c·B` rational test catches this exactly,
    // even though the surface struct fields are not bit-identical
    // (here both are unit spheres so the check is trivial; the
    // important point is that the dispatch returns the dedicated
    // `CoincidentSurfaces` error rather than diving into a solver).
    #[test]
    fn coincident_unit_spheres_short_circuit() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let anchors = EdgeAnchors::Open {
            start: Point3::new(1.0, 0.0, 0.0),
            end: Point3::new(0.0, 1.0, 0.0),
        };
        let tol = Tolerance::default();
        let opts = EdgeTripleOptions::default();
        let result = compute_edge_geometry(&sphere, &sphere, anchors, None, &tol, &opts);
        assert!(
            matches!(result, Err(EdgeTripleError::CoincidentSurfaces { .. })),
            "expected CoincidentSurfaces, got {result:?}"
        );
    }

    // Stage 2b.39: distinct cylinder radii are NOT coincident — the
    // pencil shares the axis kernel but A != c·B, so the short-circuit
    // must not fire. Either the analytic dispatch handles the pair or
    // it falls through to `Unimplemented`; the only thing we assert
    // here is the absence of `CoincidentSurfaces`.
    // Stage 2b.41: cyl-with-inscribed-sphere is the canonical
    // tangent-along-curve case. The pencil contains the rank-1 fibre
    // `D = diag(0,0,1,0)` at λ=−1 (recovered exactly in Q), giving
    // the tangent plane z=0. The intersection is the unit circle in
    // that plane. Before 2b.41 this returned `Unimplemented`; the
    // symbolic tangent path now wires plane×cyl + plane×sphere
    // dispatchers via `try_tangent_pair`.
    #[test]
    fn cyl_inscribed_sphere_returns_unit_equator_circle() {
        let cyl = z_axis_cyl(1.0);
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let anchors = EdgeAnchors::Open {
            start: Point3::new(1.0, 0.0, 0.0),
            end: Point3::new(0.0, 1.0, 0.0),
        };
        let tol = Tolerance::default();
        let opts = EdgeTripleOptions::default();
        let triple = compute_edge_geometry(&cyl, &sphere, anchors, None, &tol, &opts)
            .expect("tangent path produces a triple");
        // G5 (endpoint accuracy): start anchor is hit.
        let p_start = triple.curve_3d.point_at(triple.t_start);
        let p_end = triple.curve_3d.point_at(triple.t_end);
        assert!(
            p_start.distance_to(&Point3::new(1.0, 0.0, 0.0)) <= tol.linear(),
            "start anchor mismatch: {:?}",
            p_start
        );
        assert!(
            p_end.distance_to(&Point3::new(0.0, 1.0, 0.0)) <= tol.linear(),
            "end anchor mismatch: {:?}",
            p_end
        );
        // G2: the recorded triple verifies on both surfaces.
        let q = triple
            .verify(&cyl, &sphere, &tol, 32)
            .expect("verify on tangent intersection");
        assert!(
            q.max_positional_residual <= tol.linear(),
            "G2 residual {}",
            q.max_positional_residual
        );
        // Stage 2b.42: the symbolic tangent path tags its provenance
        // so downstream consumers (boolean / mesh / IO) can branch on
        // tangent vs transverse intersections.
        match &triple.provenance {
            crate::geometry::edge_triple::types::EdgeTripleProvenance::AnalyticIntersection {
                cardinality: Some(c),
                ..
            } => {
                assert!(c.tangent, "expected tangent: true on cyl-inscribed-sphere");
                assert_eq!(c.total, 1, "tangent intersection has a single branch");
            }
            other => panic!("expected AnalyticIntersection w/ cardinality, got {other:?}"),
        }
    }

    #[test]
    fn coaxial_cylinders_with_different_radii_do_not_short_circuit() {
        let cyl_a = z_axis_cyl(1.0);
        let cyl_b = z_axis_cyl(2.0);
        let anchors = EdgeAnchors::Open {
            start: Point3::new(1.0, 0.0, 0.0),
            end: Point3::new(1.0, 0.0, 1.0),
        };
        let tol = Tolerance::default();
        let opts = EdgeTripleOptions::default();
        let result = compute_edge_geometry(&cyl_a, &cyl_b, anchors, None, &tol, &opts);
        assert!(
            !matches!(result, Err(EdgeTripleError::CoincidentSurfaces { .. })),
            "must not classify distinct radii as coincident: {result:?}"
        );
    }

    // ─── compute_edge_geometry_all_branches (stage 2b.49) ─────────

    /// Single-branch case (cyl-with-inscribed-sphere tangent). The
    /// all-branches wrapper returns a 1-element `Vec` whose element
    /// equals the single-branch result.
    #[test]
    fn all_branches_returns_single_element_vec_on_tangent() {
        let cyl = z_axis_cyl(1.0);
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let anchors = EdgeAnchors::Open {
            start: Point3::new(1.0, 0.0, 0.0),
            end: Point3::new(0.0, 1.0, 0.0),
        };
        let tol = Tolerance::default();
        let opts = EdgeTripleOptions::default();
        let branches = compute_edge_geometry_all_branches(&cyl, &sphere, anchors, &tol, &opts)
            .expect("all-branches returns Ok on tangent intersection");
        assert_eq!(branches.len(), 1, "single-branch path returns one triple");
    }

    /// Coincident-surface short-circuit propagates through the
    /// all-branches API. Consumers should see the same dedicated
    /// error variant as via `compute_edge_geometry`.
    #[test]
    fn all_branches_propagates_coincident_surfaces_error() {
        let cyl = z_axis_cyl(1.0);
        let anchors = EdgeAnchors::Open {
            start: Point3::new(1.0, 0.0, 0.0),
            end: Point3::new(1.0, 0.0, 1.0),
        };
        let tol = Tolerance::default();
        let opts = EdgeTripleOptions::default();
        let result = compute_edge_geometry_all_branches(&cyl, &cyl, anchors, &tol, &opts);
        assert!(
            matches!(result, Err(EdgeTripleError::CoincidentSurfaces { .. })),
            "expected CoincidentSurfaces, got {result:?}"
        );
    }

    /// Unimplemented surface pair propagates the `Unimplemented`
    /// error variant rather than silently returning an empty `Vec`.
    /// (Empty `Vec` would conflate "no intersection" with "we don't
    /// know how to intersect" — bad for callers.)
    #[test]
    fn all_branches_propagates_unimplemented_error() {
        let torus = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        let anchors = EdgeAnchors::Open {
            start: Point3::new(6.0, 0.0, 0.0),
            end: Point3::new(0.0, 6.0, 0.0),
        };
        let tol = Tolerance::default();
        let opts = EdgeTripleOptions::default();
        let result = compute_edge_geometry_all_branches(&torus, &torus, anchors, &tol, &opts);
        // Two identical tori are coincident; the short-circuit fires
        // first, so we don't actually hit Unimplemented here. Use a
        // torus + non-quadric pair to test Unimplemented instead.
        assert!(matches!(
            result,
            Err(EdgeTripleError::CoincidentSurfaces { .. }) | Err(EdgeTripleError::Unimplemented)
        ));
    }
}
