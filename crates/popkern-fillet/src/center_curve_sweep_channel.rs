//! PK-141 edge-span sweep channel for legacy offset-intersect
//! center curves.
//!
//! Offset intersections computed over padded parameter ranges return
//! the whole carrier (measured: a boolean body's plane×cylinder pair
//! produced the full z∈[−5,7] line for a 0.5-long edge; sections
//! outside the edge got clamped garbage footprints and the sweep
//! twisted — SelfIntersectingBlend). The LEGACY center/footprint
//! artifacts keep the full carrier range — the face-trim splitter
//! relies on trims extending past face boundaries (measured:
//! trimming them broke a pre-filleted lip union split, and bare-span
//! trimming regressed the Gregory mixed-radius setbacks). Only the
//! SWEEP consumes this channel, and only for GROSS overhangs
//! (carrier > 3× span + one radius of headroom per side).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, NurbsCurve, Point3};

use crate::center_curve_helpers::{fit_nurbs_through_points, project_onto_surface};

/// Sweep-channel artifacts: `(trimmed center curve, footprints on
/// surface 1, footprints on surface 2)` — all `None` when no trim
/// applied.
pub(crate) type SweepChannelArtifacts =
    (Option<NurbsCurve>, Option<Vec<Point3>>, Option<Vec<Point3>>);

/// Trimmed sweep span when the carrier grossly overhangs the edge;
/// `None` keeps the legacy full-range behavior.
pub(crate) fn gross_overhang_sweep_span(
    center_curve: &NurbsCurve,
    edge_anchor_at: &dyn Fn(f64) -> Point3,
    radius: f64,
) -> Option<(f64, f64)> {
    let (full_min, full_max) = center_curve.parameter_range();
    let anchor_start = edge_anchor_at(0.0);
    let anchor_end = edge_anchor_at(1.0);
    let t_a = crate::center_curve_station::closest_t_on_nurbs(center_curve, &anchor_start);
    let t_b = crate::center_curve_station::closest_t_on_nurbs(center_curve, &anchor_end);
    let (lo, hi) = if t_a <= t_b { (t_a, t_b) } else { (t_b, t_a) };
    let span_valid = hi - lo > (full_max - full_min) * 1.0e-6;
    let chord = anchor_start.distance_to(&anchor_end).max(1.0e-12);
    let pad = (hi - lo) * (radius / chord);
    let padded = ((lo - pad).max(full_min), (hi + pad).min(full_max));
    (span_valid && (full_max - full_min) > 3.0 * (padded.1 - padded.0)).then_some(padded)
}

/// Trimmed-center refit + aligned footprints for the sweep channel.
#[allow(clippy::too_many_arguments)]
pub(crate) fn sweep_channel_artifacts(
    center_curve: &NurbsCurve,
    sweep_span: Option<(f64, f64)>,
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    n_samples: usize,
) -> Result<SweepChannelArtifacts, String> {
    let Some((lo, hi)) = sweep_span else {
        return Ok((None, None, None));
    };
    let mut centers = Vec::with_capacity(n_samples);
    let mut f1 = Vec::with_capacity(n_samples);
    let mut f2 = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = lo + (hi - lo) * (i as f64) / ((n_samples - 1) as f64);
        let center_pt = center_curve.point_at(t);
        f1.push(project_onto_surface(
            surface1, u1_range, v1_range, &center_pt,
        )?);
        f2.push(project_onto_surface(
            surface2, u2_range, v2_range, &center_pt,
        )?);
        centers.push(center_pt);
    }
    Ok((
        Some(fit_nurbs_through_points(&centers)?),
        Some(f1),
        Some(f2),
    ))
}
