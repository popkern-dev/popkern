//! A pipe along a circular arc or circle is exactly a torus (or a
//! sector of one): the section circle at the path's start, square to
//! it, revolved about the arc's axis through its angle. The general
//! pipe skins its surface through sampled frames, which read a tube
//! round a circle of radius 4 0.06% short; this path is exact (the
//! revolve recognises the torus), and every other path keeps the
//! general pipe.

use std::f64::consts::TAU;

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::ReplayError;

/// The tube of `radius` along `path` when the path is an arc or circle;
/// `None` for any other curve.
pub(super) fn pipe_along_arc(
    path: &CurveOnEdge,
    radius: f64,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
) -> Result<Option<BRep>, ReplayError> {
    let (center, normal, bend) = match path.geometry {
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        }
        | CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => (center, normal, radius),
        _ => return Ok(None),
    };
    let fail = |detail: String| ReplayError::FeatureFailed {
        label: "pipe".into(),
        error: crate::failure::FeatureFailure::sweep("pipe", &detail).into(),
    };
    if radius >= bend {
        return Err(fail(format!(
            "the tube radius {radius} reaches the path's centre of curvature (bend radius {bend}): the tube would cross itself"
        )));
    }
    let (t0, t1) = (path.t_start, path.t_end);
    let sweep = (t1 - t0).abs().min(TAU);
    let start = path.geometry.point_at(t0);
    let tangent = path.geometry.tangent_at(t0);
    let Some(tangent) = tangent.normalized() else {
        return Ok(None);
    };
    // Turning about the normal in the direction the path runs.
    let Some(axis) = normal.normalized() else {
        return Ok(None);
    };
    let axis: Vec3 = if t1 >= t0 { axis } else { -axis };
    let section = CurveGeometry::Circle {
        center: start,
        normal: tangent,
        radius,
    }
    .to_nurbs();
    popkern_sweep::revolve(
        &popkern_sweep::SweepProfile::Wire(section),
        &center,
        &axis,
        sweep,
        *tolerance,
        Some(tracking),
    )
    .map(|r| Some(r.brep))
    .map_err(|e| fail(e.to_string()))
}
