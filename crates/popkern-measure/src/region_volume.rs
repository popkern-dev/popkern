//! The volume of each solid region of a body (PK-1107).
//!
//! A region is bounded by the sides of faces facing it: a face whose
//! back lies in the region bounds it with its folded normal pointing
//! out, a face whose front lies in it with its normal pointing in. Each
//! face's divergence term (`⅓ ∬ r · n dA`, exact on planes, quadrature
//! on curved surfaces) therefore adds to the region behind it and
//! subtracts from the region in front, so an internal face between two
//! solids counts once for each, and a lamina inside a region (both
//! sides in it) cancels.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind, RegionHandle};

use crate::error::MeasureError;

/// Every solid region's volume, in region order.
///
/// # Errors
///
/// [`MeasureError`] for a body that is not a solid, general or empty
/// (the kind gate), or a tessellation failure.
pub fn region_volumes(
    brep: &BRep,
    _tol: &Tolerance,
) -> Result<Vec<(RegionHandle, f64)>, MeasureError> {
    brep.require_kind(
        "region volumes",
        &[BodyKind::Solid, BodyKind::General, BodyKind::Empty],
    )?;
    let solids = brep.solid_regions();
    if solids.is_empty() {
        return Ok(Vec::new());
    }
    let mesh = crate::body::tessellate_fine(brep)?;
    let gate = 100.0 * brep.tolerance.linear().max(1e-9);
    let terms = crate::body_quadrature::quadrature_face_volumes(brep, &mesh, gate);
    let mut volumes: Vec<(RegionHandle, f64)> = solids.iter().map(|&r| (r, 0.0)).collect();
    for (face, term) in terms {
        let Some(sides) = brep.face_sides(face) else {
            continue;
        };
        for (region, volume) in &mut volumes {
            if brep.shell_region(sides.back) == Some(*region) {
                *volume += term;
            }
            if brep.shell_region(sides.front) == Some(*region) {
                *volume -= term;
            }
        }
    }
    Ok(volumes)
}

#[cfg(test)]
#[path = "region_volume_tests.rs"]
mod tests;
