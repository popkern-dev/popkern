//! Half-space split of a closed solid by a plane (the planar arm of §2.3).
//!
//! [`split_solid_by_plane`] produces up to two closed sub-solids from a
//! single input solid and an infinite cutting plane, using the two
//! plane half-spaces and the boolean subtraction path. Curved
//! sheets require the connected-component approach in
//! [`crate::sheet::split_curved`].

use popkern_core::topology::BRep;

/// Split a solid by an infinite plane into closed half-solids.
///
/// Each half is the solid minus the opposite half-space
/// ([`crate::ops::half_space`]: a prism whose cap lies exactly on the
/// plane, oriented to it and sized from the solid's conservative face
/// boxes). Returns the pieces; a single clone of the input when the
/// plane misses the body. The cutting boxes used to be axis-aligned and
/// sized from the vertices: a tilted plane left no box face on it, and
/// a cylinder's single rim vertex under-sized the box.
pub fn split_solid_by_plane(
    solid: &BRep,
    plane_origin: popkern_core::geometry::Point3,
    plane_normal: popkern_core::geometry::Vec3,
    tol: popkern_core::tolerance::Tolerance,
) -> Result<Vec<BRep>, crate::BoolError> {
    if solid.faces.is_empty() {
        return Err(crate::BoolError::DegenerateInput(
            "solid has no faces".to_string(),
        ));
    }
    let normal = plane_normal * -1.0;
    let below = crate::ops::half_space::half_space(solid, plane_origin, plane_normal, tol)?;
    let above = crate::ops::half_space::half_space(solid, plane_origin, normal, tol)?;
    let (Some(below), Some(above)) = (below, above) else {
        // The solid lies wholly on one side: no cut.
        return Ok(vec![solid.clone()]);
    };
    // Each half by SUBTRACTING the opposite half-space: `solid - below`
    // leaves the positive half.
    let mut pieces = Vec::new();
    if let Ok(positive) = crate::subtraction_ref(solid, &below, tol, None)
        && !positive.faces.is_empty()
    {
        pieces.push(positive);
    }
    if let Ok(negative) = crate::subtraction_ref(solid, &above, tol, None)
        && !negative.faces.is_empty()
    {
        pieces.push(negative);
    }
    Ok(pieces)
}
