//! Torus strips ending on a MERIDIAN end plane (PK-769): an arc rim
//! filleted alone, its strip an analytic torus, running out square on a
//! flat face that contains the torus axis (the D profile's flat side).
//! There the terminal section is the exact minor circle, centred on the
//! spine circle at that station, in the plane normal to the spine
//! tangent. Locally that is the section of the cylinder with that
//! centre, axis (the spine tangent) and radius (the minor radius), so
//! the cylinder terminal-arc repair settles it unchanged.
//!
//! Any other end plane cuts the torus in a quartic, not a circle: that
//! end declines as `NonCylinderStrip` (fixture-less; the stale stubs
//! stay and the pipeline's curve-surface gate reports them).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::terminal_arc_trim::{CornerCandidate, StripCylinder};
use crate::types::TerminalArcTrimDeclineReason;

/// Distance from `p` to the torus `(center, unit axis, major, minor)`.
pub(super) fn torus_distance(
    p: &Point3,
    center: &Point3,
    axis: &Vec3,
    major: f64,
    minor: f64,
) -> f64 {
    let w = *p - *center;
    let axial = w.dot(axis);
    let radial = (w - *axis * axial).length();
    (((radial - major).powi(2) + axial * axial).sqrt() - minor).abs()
}

/// The local section cylinder of the torus strip `c.strip_face` at its
/// stale corner, when the strip's end plane (the stubs' common cap)
/// contains the torus axis.
///
/// # Errors
///
/// `NonCylinderStrip` when the end plane is not a meridian plane (or
/// the strip is not a torus); `NoCommonPlanarCap` when the stubs do not
/// bound one planar face.
pub(super) fn meridian_section_cylinder(
    brep: &BRep,
    c: &CornerCandidate,
    tolerance: &Tolerance,
) -> Result<StripCylinder, TerminalArcTrimDeclineReason> {
    let geometry = |face: FaceHandle| {
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
    };
    let Some(SurfaceGeometry::Torus {
        center,
        axis,
        major_radius,
        minor_radius,
    }) = geometry(c.strip_face)
    else {
        return Err(TerminalArcTrimDeclineReason::NonCylinderStrip);
    };
    let twin_face = |he: HalfEdgeHandle| -> Option<FaceHandle> {
        let tw = brep.half_edges.get(brep.half_edges.get(he)?.twin)?;
        Some(brep.loops.get(tw.loop_)?.face)
    };
    let cap = twin_face(c.he_in)
        .filter(|f| twin_face(c.he_out) == Some(*f))
        .ok_or(TerminalArcTrimDeclineReason::NoCommonPlanarCap)?;
    let Some(SurfaceGeometry::Plane { origin, normal, .. }) = geometry(cap) else {
        return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap);
    };
    let (Some(a), Some(n)) = (axis.normalized(), normal.normalized()) else {
        return Err(TerminalArcTrimDeclineReason::NonCylinderStrip);
    };
    // Meridian: the plane holds the axis (its normal ⟂ the axis, the
    // torus centre on it).
    let meridian = n.dot(&a).abs() <= 1e-9 && (center - origin).dot(&n).abs() <= tolerance.linear();
    if !meridian {
        return Err(TerminalArcTrimDeclineReason::NonCylinderStrip);
    }
    // The spine point in this plane on the corner's side of the axis.
    let corner = brep
        .vertices
        .get(c.corner)
        .map(|v| v.position)
        .ok_or(TerminalArcTrimDeclineReason::NonCylinderStrip)?;
    let w = corner - center;
    let radial = (w - a * w.dot(&a))
        .normalized()
        .ok_or(TerminalArcTrimDeclineReason::NonCylinderStrip)?;
    Ok(StripCylinder {
        origin: center + radial * major_radius,
        axis: n,
        radius: minor_radius,
    })
}
