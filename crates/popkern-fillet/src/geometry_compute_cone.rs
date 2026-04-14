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
        })?;
    let curve_on_edge = edge
        .curve
        .as_ref()
        .ok_or_else(|| BlendFailure::InvalidTopology {
            reason: "cone fillet requires the edge to carry curve geometry".to_string(),
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
                });
            }
            let pts = pts.clone();
            Box::new(move |t| center_curve::bypoints_local_radius(&pts, t))
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
    )
    .map_err(|e| BlendFailure::IntersectionFailed {
        reason: e,
        location: Some(edge_midpoint),
    })
}
