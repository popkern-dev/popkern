//! The merged cell complex (PK-964).
//!
//! Both operands of a boolean are imprinted into ONE structure in
//! which every edge keeps a radially ordered cycle of **fins** (one fin
//! per face use of the edge), every face has two **sides**, and the
//! face sides connected through the wedges between radially adjacent
//! fins form **shells**: the boundary components of the regions of
//! space the operands divide. Region membership (which operands contain
//! a region) is decided by propagation over owned face sides, not by a
//! ray per fragment, and the result is the boundary of the selected
//! regions — extracted with an explicit pairing at every edge, so a
//! non-manifold edge (four or more fins) is resolved by radial order
//! instead of by the order fragments happen to arrive in.
//!
//! References: Weiler 1986 (radial edge structure); Lee & Lee 2001
//! (partial entities: a fin is a p-edge, a face side a p-face);
//! Masuda 1993 (merge/extract split; owner attributes); Rossignac &
//! O'Connor 1989 (subdivide / select / simplify).
//!
//! Module boundary: `model`, `ids`, `operand_set` and `validate` use
//! only `popkern_core` and each other, so the data-structure half can
//! move to `popkern_core` unchanged when a second crate needs it
//! (plan §4.13; pinned by `structural_complex_module_boundary`).

pub(crate) mod artifact;
#[cfg(test)]
pub(crate) mod casting_scenes;
pub(crate) mod fans;
pub(crate) mod hole_phantoms;
pub(crate) mod ids;
pub(crate) mod lift;
pub(crate) mod model;
pub(crate) mod operand_set;
pub(crate) mod pipeline;
pub(crate) mod radial;
pub(crate) mod radial_ladder;
mod radial_stations;
pub(crate) mod record_use;
pub(crate) mod regions;
pub(crate) mod same_domain;
pub(crate) mod same_domain_check;
pub(crate) mod select;
pub(crate) mod serve;
pub(crate) mod shadow_compare;
pub(crate) mod spans;
pub(crate) mod split_plan;
pub(crate) mod validate;

#[cfg(test)]
#[path = "complex_tests.rs"]
pub(crate) mod tests;
