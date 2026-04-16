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

    // §12.3: project the spine MIDPOINT onto both surfaces to
    // get an initial UV seed. Clamped to the face UV ranges so
    // periodic surfaces (torus) don't wrap to the wrong sheet.
    let initial_seed = {
        use popkern_core::geometry::SurfaceEval;
        use popkern_core::geometry::curve::CurveEval;
        let t_mid = 0.5 * (curve_on_edge.t_start + curve_on_edge.t_end);
        let p_mid = curve_on_edge.geometry.point_at(t_mid);
        // The spine midpoint sits exactly on the torus major
        // circle, where project_point_to_uv degenerates (atan2(0,0)
        // → v=0 instead of v=π/2). Nudge the query point slightly
        // toward each surface's inward normal to break the symmetry
        // and get a well-defined minor angle.
        let nudge = 1.0e-4;
        let n1_at_mid = analysis.surface_1.normal_at(
            0.5 * (u1_range.0 + u1_range.1),
            0.5 * (v1_range.0 + v1_range.1),
        );
        let p_nudged_0 = p_mid + n1_at_mid * nudge;
        let n2_at_mid = analysis.surface_2.normal_at(
            0.5 * (u2_range.0 + u2_range.1),
            0.5 * (v2_range.0 + v2_range.1),
        );
        let p_nudged_1 = p_mid + n2_at_mid * nudge;
        let uv0 = analysis.surface_1.closest_point(&p_nudged_0);
        let uv1 = analysis.surface_2.closest_point(&p_nudged_1);
        match (uv0, uv1) {
            (Some((u0, v0)), Some((u1, v1))) => {
                let cu0: f64 = u0.clamp(u1_range.0, u1_range.1);
                let cv0: f64 = v0.clamp(v1_range.0, v1_range.1);
                let cu1: f64 = u1.clamp(u2_range.0, u2_range.1);
                let cv1: f64 = v1.clamp(v2_range.0, v2_range.1);
                Some(((cu0, cv0), (cu1, cv1)))
            }
            _ => None,
        }
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
        (curve_on_edge.t_start, curve_on_edge.t_end),
        tolerance,
        analysis.orientation_1,
        analysis.orientation_2,
        initial_seed,
    )
    .map_err(|e| BlendFailure::IntersectionFailed {
        reason: e,
        location: Some(edge_midpoint),
    })
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
            analysis.orientation_1,
            analysis.orientation_2,
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
            analysis.orientation_1,
            analysis.orientation_2,
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
                analysis.orientation_1,
                analysis.orientation_2,
            )
            .map_err(wrap_err)
        }
    }
}
