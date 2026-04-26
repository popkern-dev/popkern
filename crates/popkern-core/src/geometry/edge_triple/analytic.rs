//! Analytic dispatch for the EdgeTriple factory.
//!
//! Stage 2a: plane-pair specialisations (plane×plane, plane×cylinder,
//! plane×sphere, plane×cone). These reuse the existing analytic
//! fast paths in `geometry::intersection::surface_surface`, then
//! adapt the resulting [`CurveOnSurface`] into an [`EdgeTriple`]
//! that satisfies hard guarantees G1–G6.
//!
//! Stage 2b will add closed-form dispatches for sphere×sphere,
//! cylinder×sphere, cylinder×cylinder (parallel/coplanar), and
//! cylinder×cone (coplanar). Stages 4+ add torus and Gregory.

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::edge_triple::analytic_helpers::{
    derive_search_box, pick_branch, rebuild_consistent_cos, reverse_nurbs, slice_curve,
    snap_endpoints_3d, snap_endpoints_pcurve,
};
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, EdgeHint, SurfaceSide};
use crate::geometry::edge_triple::options::EdgeTripleOptions;
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance, EdgeTripleQuality,
};
use crate::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, intersect_surfaces,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// Try to dispatch via the analytic whitelist.
///
/// Returns:
/// - `Some(Ok(triple))` — analytic path succeeded.
/// - `Some(Err(e))` — analytic path was tried and failed (e.g.
///   anchors off surface, no intersection, ambiguous branch with no
///   hint). Caller should NOT fall through to other dispatch paths.
/// - `None` — the pair is not in the analytic whitelist for this
///   stage. Caller should continue to other dispatch arms.
pub(super) fn try_dispatch_analytic(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    tol: &Tolerance,
    opts: &EdgeTripleOptions,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let pair = classify_pair(surf_a, surf_b)?;
    Some(dispatch_analytic_inner(
        surf_a, surf_b, anchors, hint, tol, opts, pair,
    ))
}

/// Identify the surface pair if it has an analytic specialisation
/// available at the current stage.
///
/// Stage 2a ships plane×plane only. The trig pairs (plane×cylinder,
/// plane×sphere, plane×cone) require a special pcurve handling —
/// rational NURBS circles are not angle-linear except at cardinal
/// midpoints, so a same-knots/weights pcurve does not satisfy
/// `surf(pcurve(t)) ≡ curve_3d(t)` strictly. Stage 2b will
/// introduce a high-degree NURBS pcurve fit (or a non-NURBS
/// "trigonometric pcurve" type) to close that gap.
fn classify_pair(surf_a: &SurfaceGeometry, surf_b: &SurfaceGeometry) -> Option<AnalyticPair> {
    use SurfaceGeometry::*;
    match (surf_a, surf_b) {
        (Plane { .. }, Plane { .. }) => Some(AnalyticPair::PlanePlane),
        // Stage 2b adds: PlaneCylinder, PlaneSphere, PlaneCone,
        // SphereSphere, CylinderSphere, CylinderCylinder
        // (parallel/coplanar), CylinderCone (coplanar). Until then,
        // these fall through to Unimplemented at the dispatch entry.
        _ => None,
    }
}

#[allow(clippy::too_many_arguments)]
fn dispatch_analytic_inner(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    tol: &Tolerance,
    opts: &EdgeTripleOptions,
    pair: AnalyticPair,
) -> Result<EdgeTriple, EdgeTripleError> {
    if opts.cancel.as_ref().is_some_and(|c| c.is_cancelled()) {
        return Err(EdgeTripleError::Cancelled);
    }

    let (u1, v1) = derive_search_box(surf_a, anchors);
    let (u2, v2) = derive_search_box(surf_b, anchors);

    let result = intersect_surfaces(surf_a, u1, v1, surf_b, u2, v2, tol).map_err(|_| {
        EdgeTripleError::DegenerateSurfacePair {
            reason: "intersect_surfaces returned an internal error",
        }
    })?;

    let curves: Vec<&CurveOnSurface> = result
        .intersections
        .iter()
        .filter_map(|i| match i {
            SurfaceSurfaceIntersection::Curve(c) => Some(c),
            _ => None,
        })
        .collect();

    if curves.is_empty() {
        return Err(EdgeTripleError::NoIntersection);
    }

    let picked = pick_branch(&curves, anchors, hint, opts)?;
    let chosen = curves[picked];
    let branch = if curves.len() == 1 {
        BranchSelector::OnlyBranch
    } else {
        BranchSelector::Indexed(picked as u32)
    };

    // Rebuild pcurves to share curve_3d's parameterisation exactly.
    // `intersect_surfaces` sometimes returns mismatched degree/knot
    // pcurves (e.g. degree-1 hline UV against degree-2 rational
    // 3-D circle in `case_perpendicular`); the factory's contract
    // (G2/G3) requires shared parameterisation.
    let chosen_repaired = rebuild_consistent_cos(chosen, surf_a, surf_b);

    cos_to_triple(
        &chosen_repaired,
        anchors,
        EdgeTripleProvenance::AnalyticIntersection {
            kind: pair,
            branch,
            cardinality: None,
        },
        surf_a,
        surf_b,
        tol,
    )
}

/// Convert a [`CurveOnSurface`] from the intersection solver into
/// an [`EdgeTriple`] satisfying G1, G3, G5: trim to the anchor
/// parameter range, reverse if direction is wrong, snap endpoint
/// CPs to anchor positions exactly.
fn cos_to_triple(
    cos: &CurveOnSurface,
    anchors: &EdgeAnchors,
    provenance: EdgeTripleProvenance,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tol: &Tolerance,
) -> Result<EdgeTriple, EdgeTripleError> {
    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        EdgeAnchors::Closed { sample, .. } => {
            return Ok(closed_loop_triple(cos, *sample, provenance));
        }
        EdgeAnchors::Seam { .. } => {
            return Err(EdgeTripleError::DegenerateSurfacePair {
                reason: "seam anchors not supported by analytic dispatch (stage 8)",
            });
        }
    };

    // 1. Find anchor parameters on the 3-D curve.
    let t_start = cos
        .curve_3d
        .closest_point(&start)
        .ok_or(EdgeTripleError::AnchorOffSurface {
            which: SurfaceSide::A,
            distance: f64::INFINITY,
        })?;
    let t_end = cos
        .curve_3d
        .closest_point(&end)
        .ok_or(EdgeTripleError::AnchorOffSurface {
            which: SurfaceSide::B,
            distance: f64::INFINITY,
        })?;

    // 2. Verify anchors lie on the curve within tolerance.
    let d_start = cos.curve_3d.point_at(t_start).distance_to(&start);
    let d_end = cos.curve_3d.point_at(t_end).distance_to(&end);
    if d_start > tol.linear() {
        return Err(EdgeTripleError::AnchorOffSurface {
            which: SurfaceSide::A,
            distance: d_start,
        });
    }
    if d_end > tol.linear() {
        return Err(EdgeTripleError::AnchorOffSurface {
            which: SurfaceSide::B,
            distance: d_end,
        });
    }

    // 3. Order parameters; if reversed, flip all three curves.
    let (t_lo, t_hi, reversed) = if t_start <= t_end {
        (t_start, t_end, false)
    } else {
        (t_end, t_start, true)
    };

    let (mut curve_3d, mut pcurve_a, mut pcurve_b) = (
        slice_curve(&cos.curve_3d, t_lo, t_hi).ok_or(EdgeTripleError::DegenerateSurfacePair {
            reason: "curve_3d slice failed",
        })?,
        slice_curve(&cos.curve_uv1, t_lo, t_hi).ok_or(EdgeTripleError::DegenerateSurfacePair {
            reason: "pcurve_a slice failed",
        })?,
        slice_curve(&cos.curve_uv2, t_lo, t_hi).ok_or(EdgeTripleError::DegenerateSurfacePair {
            reason: "pcurve_b slice failed",
        })?,
    );

    if reversed {
        curve_3d = reverse_nurbs(&curve_3d);
        pcurve_a = reverse_nurbs(&pcurve_a);
        pcurve_b = reverse_nurbs(&pcurve_b);
    }

    // 4. Snap endpoint CPs to anchors exactly (G1).
    snap_endpoints_3d(&mut curve_3d, start, end);
    snap_endpoints_pcurve(&mut pcurve_a, surf_a, start, end);
    snap_endpoints_pcurve(&mut pcurve_b, surf_b, start, end);

    let (t_start_new, t_end_new) = curve_3d.parameter_range();

    use crate::geometry::curve::CurveGeometry;
    use crate::geometry::edge_triple::types::ParameterMap;
    let triple = EdgeTriple {
        curve_3d: CurveGeometry::Nurbs(curve_3d),
        pcurve_a: pcurve_a.into(),
        pcurve_b: pcurve_b.into(),
        t_start: t_start_new,
        t_end: t_end_new,
        topology: EdgeTopology::Open,
        // Stage 2a's plane×plane and other intersect_surfaces-based
        // dispatches share knots across the three NURBS — the
        // legacy invariant. Stage 2b.3+ trig dispatches use Identity
        // with analytic curves directly.
        parameter_map: ParameterMap::SharedKnots,
        provenance,
        quality: EdgeTripleQuality::default(),
    };

    finalize_quality(triple, surf_a, surf_b, tol)
}

fn closed_loop_triple(
    cos: &CurveOnSurface,
    _sample: Point3,
    provenance: EdgeTripleProvenance,
) -> EdgeTriple {
    use crate::geometry::curve::CurveGeometry;
    use crate::geometry::edge_triple::types::ParameterMap;
    let (t_start, t_end) = cos.curve_3d.parameter_range();
    EdgeTriple {
        curve_3d: CurveGeometry::Nurbs(cos.curve_3d.clone()),
        pcurve_a: cos.curve_uv1.clone().into(),
        pcurve_b: cos.curve_uv2.clone().into(),
        t_start,
        t_end,
        topology: EdgeTopology::Closed {
            period: t_end - t_start,
        },
        parameter_map: ParameterMap::SharedKnots,
        provenance,
        quality: EdgeTripleQuality::default(),
    }
}

/// Crate-internal alias used by `analytic_trig.rs` so it can run
/// the same finalize-and-verify pass without duplicating logic.
pub(super) fn finalize_quality_pub(
    triple: EdgeTriple,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tol: &Tolerance,
) -> Result<EdgeTriple, EdgeTripleError> {
    finalize_quality(triple, surf_a, surf_b, tol)
}

fn finalize_quality(
    mut triple: EdgeTriple,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tol: &Tolerance,
) -> Result<EdgeTriple, EdgeTripleError> {
    let q = triple.verify(surf_a, surf_b, tol, 32).map_err(|_| {
        EdgeTripleError::DegenerateSurfacePair {
            reason: "post-construction verify failed",
        }
    })?;
    let endpoint_err = endpoint_err_3d(&triple);
    triple.quality = EdgeTripleQuality {
        max_endpoint_error: endpoint_err,
        ..q
    };
    Ok(triple)
}

fn endpoint_err_3d(triple: &EdgeTriple) -> f64 {
    use crate::geometry::curve::CurveGeometry;
    let p_start = triple.curve_3d.point_at(triple.t_start);
    let p_end = triple.curve_3d.point_at(triple.t_end);
    // For analytic curve_3d (Line, Circle, Arc) the endpoint is
    // exact-by-construction; for NURBS we check first/last CPs.
    match &triple.curve_3d {
        CurveGeometry::Nurbs(n) => {
            if let (Some(first), Some(last)) = (n.control_points.first(), n.control_points.last()) {
                let d_start = first.distance_to(&p_start);
                let d_end = last.distance_to(&p_end);
                d_start.max(d_end)
            } else {
                f64::INFINITY
            }
        }
        _ => 0.0,
    }
}
