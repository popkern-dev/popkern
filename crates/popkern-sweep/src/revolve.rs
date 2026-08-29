//! Revolve operation — creates a solid of revolution by rotating a profile
//! around an axis.
//!
//! The profile is analyzed for smooth spans, and each span is converted to a
//! NURBS surface of revolution using [`popkern_core::geometry::revolution::revolution_surface`].
//! The resulting surfaces are assembled into a BRep using the existing sweep
//! topology builder.

use popkern_core::geometry::revolution::{revolution_surface, rotate_point};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;

use crate::error::SweepError;
use crate::profile::analyze_profile;
use crate::topology::side_faces::build_sweep_topology;
use crate::types::{SweepProfile, SweepResult};

/// Revolve a profile around an axis by a given angle.
///
/// Creates a solid of revolution by rotating the profile curve(s) around
/// the axis. A full revolution (angle = 2pi) creates a closed solid.
/// Partial revolutions create an open-sided solid with start/end cap faces.
///
/// # Arguments
/// * `profile` - The profile to revolve (wire or face).
/// * `axis_origin` - A point on the revolution axis.
/// * `axis_direction` - Direction of the revolution axis.
/// * `angle` - Revolution angle in radians, must be in (0, 2pi].
/// * `tolerance` - Geometric tolerance.
///
/// # Errors
///
/// Returns [`SweepError`] if the profile is invalid, the angle is out of
/// range, or topology construction fails.
pub fn revolve(
    profile: &SweepProfile,
    axis_origin: &Point3,
    axis_direction: &Vec3,
    angle: f64,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    revolve_impl(
        profile,
        axis_origin,
        axis_direction,
        angle,
        tolerance,
        tracking,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// Honest per-run stage count of [`revolve_impl`] for progress totals
/// (PK-276).
pub(crate) const REVOLVE_STAGE_COUNT: usize = 3;

/// Shared implementation for [`revolve`]. `tracker` threads PK-276
/// interactive control: a cancellation poll before each stage
/// (`SweepError::Cancelled`, nothing produced) and one progress event
/// per stage; the surface stage polls per profile span.
#[allow(clippy::too_many_arguments)]
pub(crate) fn revolve_impl(
    profile: &SweepProfile,
    axis_origin: &Point3,
    axis_direction: &Vec3,
    angle: f64,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<SweepResult, SweepError> {
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(SweepError::Cancelled);
            }
        };
    }
    // 1. Analyze profile into smooth spans.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let analyzed = analyze_profile(profile, &tolerance)?;

    if analyzed.spans.is_empty() {
        return Err(SweepError::InvalidProfile("no profile spans".into()));
    }
    tracker.finish_stage("revolve.analyze_profile", t0, analyzed.spans.len());

    // 2. Build a revolution surface for each span (per-span poll).
    let t0 = popkern_core::time::Stopwatch::start();
    let mut surfaces = Vec::with_capacity(analyzed.spans.len());
    for span in &analyzed.spans {
        stage_gate!();
        let surface = revolution_surface(span, axis_origin, axis_direction, angle)
            .map_err(SweepError::Core)?;
        surfaces.push(surface);
    }
    tracker.finish_stage("revolve.surfaces", t0, surfaces.len());

    // 3. Build cap curves (rotated profile endpoints) for partial revolutions.
    let is_full_revolution = (angle - std::f64::consts::TAU).abs() < 1e-10;

    let (start_caps, end_caps) = if analyzed.is_closed && !is_full_revolution {
        // For closed profiles with partial revolution, we need cap faces.
        let start_caps = build_cap_curves(&analyzed.spans);
        let end_caps =
            build_rotated_cap_curves(&analyzed.spans, axis_origin, axis_direction, angle);
        (start_caps, end_caps)
    } else if analyzed.is_closed && is_full_revolution {
        // Full revolution with closed profile: no caps needed, but the
        // topology builder expects cap curves for closed profiles.
        // We pass empty caps and treat it as an open profile to avoid
        // adding cap faces — the surface wraps around.
        (vec![], vec![])
    } else {
        (vec![], vec![])
    };

    // 4. Build topology.
    // For full revolution with closed profile, the body is already closed
    // so we build with is_closed_profile = false to avoid cap faces.
    let use_closed_topology = analyzed.is_closed && !is_full_revolution;

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let brep = build_sweep_topology(
        &surfaces,
        &start_caps,
        &end_caps,
        use_closed_topology,
        tolerance,
        tracking,
    )?;
    tracker.finish_stage("revolve.topology", t0, brep.faces.len());

    Ok(SweepResult {
        brep,
        warnings: vec![],
    })
}

/// Build cap curves from the profile spans at the start position (angle = 0).
///
/// Each span contributes one curve segment to the cap boundary.
fn build_cap_curves(spans: &[NurbsCurve]) -> Vec<NurbsCurve> {
    spans.to_vec()
}

/// Build cap curves from the profile spans rotated to the end position.
///
/// Each span's control points are rotated around the axis by the given angle.
fn build_rotated_cap_curves(
    spans: &[NurbsCurve],
    axis_origin: &Point3,
    axis_direction: &Vec3,
    angle: f64,
) -> Vec<NurbsCurve> {
    spans
        .iter()
        .map(|span| {
            let rotated_cps: Vec<Point3> = span
                .control_points
                .iter()
                .map(|cp| rotate_point(cp, axis_origin, axis_direction, angle))
                .collect();
            NurbsCurve {
                degree: span.degree,
                knots: span.knots.clone(),
                control_points: rotated_cps,
                weights: span.weights.clone(),
            }
        })
        .collect()
}

#[cfg(test)]
#[path = "revolve_tests.rs"]
mod tests;
