//! A closed sweep or loft faced out of the volume it bounds.
//!
//! The wall's faces take their sense from the sections' winding and the
//! caps from the wall, so a loft whose sections run against the normal
//! of their winding (a square at z = 0 lofted to one at z = -12, or the
//! sections of a profile on the xz plane, whose normal is -y, lofted
//! toward +y) came out with every face facing into it: a solid inside
//! out, its signed volume negative. The whole body is reversed together
//! when the sampled enclosed volume reads negative; a set too small or
//! degenerate to call is left as built.

use popkern_core::topology::BRep;

/// Reverse `brep` (every loop's traversal and every face's flag
/// together) when its faces face into the volume they bound.
pub(crate) fn orient_outward(brep: &mut BRep) {
    let faces: Vec<_> = brep.faces.keys().collect();
    if brep.folded_wedge_volume(&faces).is_some_and(|v| v < 0.0) {
        brep.reverse_body_orientation();
    }
}
