//! The pieces of a body (design: operations on the cell complex, 5.4):
//! its solids, each kept or dropped whole with the faces, edges and
//! vertices that bound it, their names and lineage unchanged.
//!
//! A body of several solids apart from each other (lumps, as
//! [`crate::concatenate_solids`] joins them) holds each with its own
//! topology. Taking pieces out therefore removes whole components;
//! nothing is cut. A
//! body whose solids share an entity (a vertex or edge used by two
//! solids) is not pieces of this kind and declines.

use popkern_core::topology::{BRep, RegionHandle};

use crate::error::BoolError;

/// The body's pieces, its solids, in slot order (the order a split
/// served them in).
pub fn piece_handles(body: &BRep) -> Vec<RegionHandle> {
    body.solid_regions()
}

/// `body` with only the solids in `keep`: every entity bounding a
/// dropped solid is removed, every kept entity keeps its handle, name
/// and geometry (`BRep::retain_pieces`).
///
/// ```
/// use popkern_bool::{concatenate_solids, pieces};
/// use popkern_core::{geometry::Transform3, tolerance::Tolerance, topology::BRep};
///
/// let tol = Tolerance::default();
/// let a = BRep::unit_cube(tol).unwrap();
/// let b = BRep::unit_cube(tol)
///     .unwrap()
///     .transform(&Transform3::translation(3.0, 0.0, 0.0))
///     .unwrap();
/// let both = concatenate_solids(a, b).unwrap();
/// let first = pieces::piece_handles(&both)[0];
/// let one = pieces::retain_pieces(&both, &[first]).unwrap();
/// assert_eq!((one.solid_regions().len(), one.faces.len()), (1, 6));
/// ```
///
/// # Errors
///
/// [`BoolError::DegenerateInput`] when a kept and a dropped solid share
/// a half-edge, edge or vertex (the body is not separable pieces).
pub fn retain_pieces(body: &BRep, keep: &[RegionHandle]) -> Result<BRep, BoolError> {
    body.retain_pieces(keep)
        .map_err(|e| BoolError::DegenerateInput(format!("{e}: its pieces cannot be taken apart")))
}

/// Each piece of `body` as its own body (see [`retain_pieces`]), in
/// [`piece_handles`] order.
///
/// # Errors
///
/// As [`retain_pieces`].
pub fn piece_bodies(body: &BRep) -> Result<Vec<BRep>, BoolError> {
    piece_handles(body)
        .into_iter()
        .map(|s| retain_pieces(body, &[s]))
        .collect()
}

#[cfg(test)]
#[path = "pieces_tests.rs"]
mod tests;
