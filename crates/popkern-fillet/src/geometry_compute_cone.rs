//! Walker-based center-curve routing. Pulls the edge curve from the
//! BRep, builds a `WalkRequest` against it, and routes through
//! `walker_based_center_curve` which handles the full radius-type
//! matrix (Constant, Variable, ByPoints).
//!
//! **When this path runs:**
//! - Phase 1b: either support surface is a Cone (offset-intersect
//!   path can't run because cones have no analytic offset).
//! - Phase 2 step 4: the caller requested G2 cross-sections via
//!   `FilletOptions::g2_cross_section`. The walker populates
//!   `CenterCurveResult::contact_solutions`, which lets the downstream
//!   `build_fillet_surface` route through `g2_sweep_from_contact_solutions`
//!   on any surface pair (not just cones).
//!
//! Split out of `geometry_compute.rs` to keep the parent under the
//! 500-line structural cap.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::analysis;
use crate::center_curve;
use crate::center_curve::CenterCurveResult;
use crate::error::BlendFailure;
use crate::types::FilletRadius;

#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_center_via_walker(
    analysis: &analysis::EdgeAnalysis,
    radius: &FilletRadius,
    brep: &BRep,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    edge_midpoint: Point3,
) -> Result<center_curve::CenterCurveResult, BlendFailure> {
    let edge = brep
        .edges
        .get(analysis.edge)
        .ok_or_else(|| BlendFailure::InvalidTopology {
            reason: format!("edge {:?} not found", analysis.edge),
            location: None,
        })?;
    let curve_on_edge = edge
        .curve
        .as_ref()
        .ok_or_else(|| BlendFailure::InvalidTopology {
            reason: "cone fillet requires the edge to carry curve geometry".to_string(),
            location: Some(edge_midpoint),
        })?;

    // Map the radius profile to a closure over the [0, 1] parameter
    // that `walker_based_center_curve` expects.
    let radius_at_normalized: Box<dyn Fn(f64) -> f64> = match radius {
        FilletRadius::Constant(r) => {
            let r = *r;
            Box::new(move |_t| r)
        }
        FilletRadius::Variable { start, end } => {
            let start = *start;
            let end = *end;
            Box::new(move |t| start + (end - start) * t)
        }
        FilletRadius::ByPoints(pts) => {
            if pts.is_empty() {
                return Err(BlendFailure::InvalidTopology {
                    reason: "ByPoints radius with empty point list".to_string(),
                    location: Some(edge_midpoint),
                });
            }
            let pts = pts.clone();
            Box::new(move |t| center_curve::bypoints_local_radius(&pts, t))
        }
    };

    // Project a query point onto both surfaces to get an initial
    // UV seed. Use the spine midpoint, nudged toward each surface
    // normal to break atan2 degeneracies on periodic surfaces.
    //
    // For blend-on-blend (§12.3), the spine midpoint may NOT lie
    // on both surfaces (e.g. edge from cube corner to cylinder
    // rim). Fall back to the edge endpoint that IS on the non-
    // planar surface by picking the vertex closest to the face
    // centroid.
    let (initial_seed, spine_t_range) = {
        use popkern_core::geometry::SurfaceEval;
        use popkern_core::geometry::curve::CurveEval;
        let ts = curve_on_edge.t_start;
        let te = curve_on_edge.t_end;
        let t_mid = 0.5 * (ts + te);
        let p_mid = curve_on_edge.geometry.point_at(t_mid);
        let p_start = curve_on_edge.geometry.point_at(ts);
        let p_end = curve_on_edge.geometry.point_at(te);

        let dist_to_surfs = |p: &Point3| -> f64 {
            let d = |s: &popkern_core::geometry::surface::SurfaceGeometry| {
                s.closest_point(p)
                    .map(|(u, v)| (s.point_at(u, v) - *p).length())
                    .unwrap_or(f64::INFINITY)
            };
            d(&analysis.surface_1).max(d(&analysis.surface_2))
        };
        let mut best_pt = p_mid;
        for candidate in [p_start, p_end] {
            if dist_to_surfs(&candidate) < dist_to_surfs(&best_pt) {
                best_pt = candidate;
            }
        }

        let nudge = 1.0e-4;
        let proj = |s: &popkern_core::geometry::surface::SurfaceGeometry,
                    ur: (f64, f64),
                    vr: (f64, f64)|
         -> Option<(f64, f64)> {
            let n = s.normal_at(0.5 * (ur.0 + ur.1), 0.5 * (vr.0 + vr.1));
            let (u, v) = s.closest_point(&(best_pt + n * nudge))?;
            Some((u.clamp(ur.0, ur.1), v.clamp(vr.0, vr.1)))
        };
        let uv0 = proj(&analysis.surface_1, u1_range, v1_range);
        let uv1 = proj(&analysis.surface_2, u2_range, v2_range);
        let seed = match (uv0, uv1) {
            (Some(a), Some(b)) => Some((a, b)),
            _ => None,
        };
        // Walk from the end nearest to the seed point so the
        // initial contact has the best chance of converging.
        let ds = (best_pt - p_start).length();
        let de = (best_pt - p_end).length();
        let t_range = if de < ds { (te, ts) } else { (ts, te) };
        (seed, t_range)
    };

    center_curve::walker_based_center_curve(
        &analysis.surface_1,
        u1_range,
        v1_range,
        &analysis.surface_2,
        u2_range,
        v2_range,
        radius_at_normalized,
        &curve_on_edge.geometry,
        spine_t_range,
        tolerance,
        analysis.orientation_1 == analysis.shell_outward,
        analysis.orientation_2 == analysis.shell_outward,
        analysis.convexity,
        initial_seed,
    )
    .map_err(|e| BlendFailure::IntersectionFailed {
        reason: e,
        location: Some(edge_midpoint),
    })
}

/// Walker-or-legacy center curve dispatch with §7 continuity
/// downgrade tracking. Returns `(center_result, downgrade_reason)` —
/// the reason string is `Some` when the walker failed on a non-cone
/// G2 request and the pipeline degraded to the legacy path.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_center_with_fallback(
    analysis: &analysis::EdgeAnalysis,
    radius: &FilletRadius,
    brep: &BRep,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    edge_midpoint: Point3,
    either_needs_walker: bool,
    wants_g2: bool,
) -> Result<(CenterCurveResult, Option<String>), BlendFailure> {
    if either_needs_walker || wants_g2 {
        let w = compute_center_via_walker(
            analysis,
            radius,
            brep,
            u1_range,
            v1_range,
            u2_range,
            v2_range,
            tolerance,
            edge_midpoint,
        );
        match w {
            Ok(r) => Ok((r, None)),
            Err(walker_err) if !either_needs_walker && wants_g2 => {
                let reason = format!("walker failed: {walker_err:?}");
                let legacy = compute_legacy_center_curve(
                    analysis,
                    radius,
                    u1_range,
                    v1_range,
                    u2_range,
                    v2_range,
                    tolerance,
                    edge_midpoint,
                )?;
                Ok((legacy, Some(reason)))
            }
            Err(e) => Err(e),
        }
    } else {
        compute_legacy_center_curve(
            analysis,
            radius,
            u1_range,
            v1_range,
            u2_range,
            v2_range,
            tolerance,
            edge_midpoint,
        )
        .map(|r| (r, None))
    }
}

/// Legacy offset-intersect center curve: the pre-phase-1 path that
/// analytically offsets both surfaces, intersects the offsets, and
/// picks the curve closest to the edge midpoint. Used as the primary
/// route when neither the cone trigger nor the G2 opt-in fires, and
/// as a fallback when the walker fails on a G2 opt-in request
/// (silently degrading the caller's G2 request to G1 rather than
/// propagating a walker-specific error).
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_legacy_center_curve(
    analysis: &analysis::EdgeAnalysis,
    radius: &FilletRadius,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    edge_midpoint: Point3,
) -> Result<CenterCurveResult, BlendFailure> {
    let wrap_err = |e: String| BlendFailure::IntersectionFailed {
        reason: e,
        location: Some(edge_midpoint),
    };
    match radius {
        FilletRadius::Constant(r) => center_curve::constant_radius_center_curve(
            &analysis.surface_1,
            u1_range,
            v1_range,
            &analysis.surface_2,
            u2_range,
            v2_range,
            *r,
            &edge_midpoint,
            tolerance,
            analysis.orientation_1 == analysis.shell_outward,
            analysis.orientation_2 == analysis.shell_outward,
            analysis.convexity,
        )
        .map_err(wrap_err),
        FilletRadius::Variable { start, end } => center_curve::variable_radius_center_curve(
            &analysis.surface_1,
            u1_range,
            v1_range,
            &analysis.surface_2,
            u2_range,
            v2_range,
            *start,
            *end,
            &edge_midpoint,
            tolerance,
            analysis.orientation_1 == analysis.shell_outward,
            analysis.orientation_2 == analysis.shell_outward,
            analysis.convexity,
        )
        .map_err(wrap_err),
        FilletRadius::ByPoints(pts) => {
            if pts.is_empty() {
                return Err(BlendFailure::InvalidTopology {
                    reason: "ByPoints radius with empty point list".to_string(),
                    location: Some(edge_midpoint),
                });
            }
            center_curve::parametric_radius_center_curve(
                &analysis.surface_1,
                u1_range,
                v1_range,
                &analysis.surface_2,
                u2_range,
                v2_range,
                |t| center_curve::bypoints_local_radius(pts, t),
                &edge_midpoint,
                tolerance,
                analysis.orientation_1 == analysis.shell_outward,
                analysis.orientation_2 == analysis.shell_outward,
                analysis.convexity,
            )
            .map_err(wrap_err)
        }
    }
}
