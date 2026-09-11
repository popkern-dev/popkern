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
use crate::topology::revolution_faces::build_full_revolution_topology;
use crate::types::{SweepProfile, SweepResult};

/// Translate a curve by `shift`.
fn shift_curve(curve: &NurbsCurve, shift: Vec3) -> NurbsCurve {
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: curve.control_points.iter().map(|p| *p + shift).collect(),
        weights: curve.weights.clone(),
    }
}

/// The same curve traversed the other way (parameter range kept).
fn reverse_curve(curve: &NurbsCurve) -> NurbsCurve {
    let (lo, hi) = (
        curve.knots.first().copied().unwrap_or(0.0),
        curve.knots.last().copied().unwrap_or(1.0),
    );
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.iter().rev().map(|k| lo + hi - k).collect(),
        control_points: curve.control_points.iter().rev().copied().collect(),
        weights: curve.weights.iter().rev().copied().collect(),
    }
}

/// Signed area of the closed profile in the (radial, axial) half-plane
/// frame of the revolution: positive when the profile runs
/// counter-clockwise with the radial direction as the abscissa and the
/// axis as the ordinate. `du × dθ` of a revolution surface points
/// INWARD for such a profile (at the outermost point the profile
/// tangent is +axis and the sweep tangent +θ, whose cross product is
/// −radial), so a positive area means the spans must be traversed the
/// other way for outward normals.
fn profile_signed_area(spans: &[NurbsCurve], axis_origin: &Point3, axis_dir: &Vec3) -> f64 {
    use popkern_core::geometry::CurveEval;
    let Some(axis) = axis_dir.normalized() else {
        return 0.0;
    };
    // Radial reference: the profile centroid's offset from the axis.
    let mut acc = Vec3::new(0.0, 0.0, 0.0);
    let mut count: f64 = 0.0;
    for span in spans {
        for i in 0..8 {
            let (t0, t1) = span.parameter_range();
            let p = span.point_at(t0 + (t1 - t0) * (i as f64 + 0.5) / 8.0);
            acc = acc + (p - *axis_origin);
            count += 1.0;
        }
    }
    let mean = acc * (1.0 / count.max(1.0));
    let radial = mean - axis * mean.dot(&axis);
    let Some(e_r) = radial.normalized() else {
        return 0.0;
    };
    let mut points: Vec<(f64, f64)> = Vec::new();
    for span in spans {
        let (t0, t1) = span.parameter_range();
        for i in 0..16 {
            let p = span.point_at(t0 + (t1 - t0) * (i as f64) / 16.0) - *axis_origin;
            points.push((p.dot(&e_r), p.dot(&axis)));
        }
    }
    let mut area = 0.0;
    for i in 0..points.len() {
        let (r0, z0) = points[i];
        let (r1, z1) = points[(i + 1) % points.len()];
        area += r0 * z1 - r1 * z0;
    }
    area * 0.5
}

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
    // The analysis centres the profile at its centroid (sweeps move it
    // along the path); a revolution keeps the profile where it is, so
    // put the spans back in world coordinates.
    let shift = analyzed.centroid - Point3::new(0.0, 0.0, 0.0);
    let mut spans: Vec<NurbsCurve> = analyzed
        .spans
        .iter()
        .map(|span| shift_curve(span, shift))
        .collect();
    // Outward normals: traverse a closed profile clockwise in the
    // (radial, axial) frame (see `profile_signed_area`).
    if analyzed.is_closed && profile_signed_area(&spans, axis_origin, axis_direction) > 0.0 {
        spans = spans.iter().rev().map(reverse_curve).collect();
    }
    // A span lying on the axis sweeps nothing: drop it. The remaining
    // chain of a closed profile then runs pole to pole (a disc, a
    // cone, a sphere); the welded full-revolution builder handles
    // both shapes, the capped partial builder cannot.
    let is_full_revolution = (angle - std::f64::consts::TAU).abs() < 1e-10;
    let on_axis = |span: &NurbsCurve| {
        use popkern_core::geometry::CurveEval;
        let (t0, t1) = span.parameter_range();
        (0..5).all(|i| {
            let p = span.point_at(t0 + (t1 - t0) * (i as f64) / 4.0);
            let d = p - *axis_origin;
            let axis = axis_direction.normalized().unwrap_or(*axis_direction);
            (d - axis * d.dot(&axis)).length() < tolerance.linear() * 10.0
        })
    };
    let dropped = spans.iter().filter(|s| on_axis(s)).count();
    if dropped > 0 {
        if !is_full_revolution {
            return Err(SweepError::InvalidProfile(
                "a profile span on the axis is only supported for full revolutions".into(),
            ));
        }
        // Rotate the chain so it starts right after a dropped span,
        // then keep the off-axis spans in order (pole to pole).
        let n = spans.len();
        let first_off = (0..n)
            .find(|&i| on_axis(&spans[(i + n - 1) % n]) && !on_axis(&spans[i]))
            .unwrap_or(0);
        spans = (0..n)
            .map(|k| spans[(first_off + k) % n].clone())
            .filter(|s| !on_axis(s))
            .collect();
        if spans.is_empty() {
            return Err(SweepError::InvalidProfile(
                "the whole profile lies on the revolution axis".into(),
            ));
        }
    }
    let wires = vec![spans.len()];
    let analyzed = crate::profile::AnalyzedProfile {
        spans,
        wires,
        is_closed: analyzed.is_closed,
        centroid: analyzed.centroid,
    };
    tracker.finish_stage("revolve.analyze_profile", t0, analyzed.spans.len());

    // 2. Build a revolution surface for each span (per-span poll).
    let t0 = popkern_core::time::Stopwatch::start();
    let mut surfaces = Vec::with_capacity(analyzed.spans.len());
    for span in &analyzed.spans {
        stage_gate!();
        let mut surface = revolution_surface(span, axis_origin, axis_direction, angle)
            .map_err(SweepError::Core)?;
        // The core surface parameterises the angle as v ∈ [0, angle]
        // radians; every topology builder evaluates corners, iso-curves
        // and pcurves in a [0, 1] chart, so normalise v here (PK-408:
        // v = 1 used to land at 57.3°, opening every loop).
        for k in &mut surface.knots_v {
            *k /= angle;
        }
        surfaces.push(surface);
    }
    tracker.finish_stage("revolve.surfaces", t0, surfaces.len());

    // 3. Build cap curves (rotated profile endpoints) for partial revolutions.

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
    let brep = if analyzed.is_closed && is_full_revolution {
        // Welded periodic faces: one seam per span, one circle per
        // profile vertex (PK-408).
        let mut brep = build_full_revolution_topology(
            &analyzed.spans,
            &surfaces,
            axis_origin,
            axis_direction,
            tolerance,
        )?;
        if let Some(ctx) = tracking {
            crate::topology::side_faces::assign_sweep_names(&mut brep, ctx);
        }
        brep
    } else {
        crate::topology::side_faces::build_sweep_topology_with(
            &surfaces,
            &start_caps,
            &end_caps,
            use_closed_topology,
            crate::topology::side_faces::SideFaceGeometry::Nurbs,
            tolerance,
            tracking,
        )?
    };
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
