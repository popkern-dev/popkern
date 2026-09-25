//! Seams on a dual cross-jet patch boundary certify against the
//! patch's polynomial boundary surrogate (PK-541).
//!
//! The seam quantities (position, normal, second form) read only the
//! side's 2-jet along its pcurve. On a patch boundary, a
//! [`RectangularGregoryBiNonicDualJet`] has the 2-jet of its plain
//! [`boundary_surrogate`] exactly (see that module). Along a pcurve
//! that lies on the boundary, every Taylor coefficient of every seam
//! quantity is therefore the surrogate's.
//!
//! The pcurve must lie on the boundary exactly, in real arithmetic:
//! - a linear pcurve with both ends on it; or
//! - a NURBS pcurve whose every control point has that coordinate, so
//!   that the curve is a positive-weight convex combination of them.
//!
//! A wider interval enclosure of such a pcurve is still sound: the
//! surrogate's enclosure over a superset contains the surrogate along
//! the real pcurve, and that equals the patch there.
//!
//! The rational exact form has weight `Π (wu² + wv²)`, which vanishes
//! at every activated corner. Near a corner its enclosures decline
//! (`DenominatorContainsZero`) or blow up (second-form bounds of 1e8
//! and more on the KP20 setback fan), even though the patch is
//! polynomial along its boundaries.
//!
//! [`boundary_surrogate`]: RectangularGregoryBiNonicDualJet::boundary_surrogate

use crate::geometry::gregory::{GregoryPatch, PatchBoundary, RectangularGregoryBiNonicDualJet};
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use crate::geometry::surface::SurfaceGeometry;

/// The patch boundary `pcurve` lies on exactly, if any.
pub(super) fn pcurve_patch_boundary(pcurve: &Pcurve) -> Option<PatchBoundary> {
    match &pcurve.geometry {
        PcurveGeometry::Linear { start, end } => {
            if start.0 == end.0 {
                PatchBoundary::of_chart_line(true, start.0)
            } else if start.1 == end.1 {
                PatchBoundary::of_chart_line(false, start.1)
            } else {
                None
            }
        }
        PcurveGeometry::Nurbs(curve) => {
            let (first, rest) = curve.control_points.split_first()?;
            if curve.weights.iter().any(|&w| w.is_nan() || w <= 0.0) {
                return None;
            }
            if rest.iter().all(|p| p.x == first.x) {
                PatchBoundary::of_chart_line(true, first.x)
            } else if rest.iter().all(|p| p.y == first.y) {
                PatchBoundary::of_chart_line(false, first.y)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// The surface a side certifies against: the dual-jet patch's boundary
/// surrogate when the pcurve lies on a boundary of a patch with an
/// activated corner, else `None` (the surface itself is used).
pub(super) fn boundary_surrogate_surface(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
) -> Option<SurfaceGeometry> {
    let SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonicDualJet(patch)) = surface else {
        return None;
    };
    let patch: &RectangularGregoryBiNonicDualJet = patch;
    if patch.corners.iter().all(Option::is_none) {
        return None;
    }
    let boundary = pcurve_patch_boundary(pcurve)?;
    Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(
        Box::new(patch.boundary_surrogate(boundary)),
    )))
}

#[cfg(test)]
#[path = "boundary_surrogate_tests.rs"]
mod tests;
