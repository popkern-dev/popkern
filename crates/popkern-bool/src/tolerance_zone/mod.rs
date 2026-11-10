//! Tolerance zones and candidate-vertex resolution (plan §4.3, §4.3a).
//!
//! Deciding vertices one interference at a time is how a boolean ends
//! up with order-dependent snaps, T-junction repair passes and open
//! intersection loops: each interference sees part of the picture.
//! This module defers the decision: every interference contributes a
//! candidate vertex (`candidates`), and resolution (`resolve`) decides
//! how many real vertices each group of nearby candidates becomes,
//! with the rule behind every separation recorded.
//!
//! **Stage.** The committed table decides the intersection graph's
//! record-end vertices and paves on every cellular boolean
//! (`table::table_enabled`). The graph's own identity
//! (`IntersectionGraph::find_node`: first match within 2 × linear,
//! record ends adopted within 20 × linear) serves only where the
//! resolution declines, a recorded fallback, or under the kill switch
//! `POPKERN_VERTEX_TABLE=0`. `shadow` compares the two on every
//! cellular boolean when `POPKERN_COMPLEX_SHADOW=1`.
//!
//! **Tolerance.** Every radius is an entity's own tolerance, the model
//! linear tolerance when it has none; the growth cap is an explicit
//! parameter. No threshold here is global state.

pub(crate) mod branches;
pub(crate) mod candidates;
pub(crate) mod closure;
pub(crate) mod common_blocks;
pub(crate) mod consequences;
pub(crate) mod face_info;
pub(crate) mod face_path;
pub(crate) mod face_side;
pub(crate) mod incidence;
pub(crate) mod meetings;
mod meetings_foot;
pub(crate) mod partition;
pub(crate) mod pave_blocks;
pub(crate) mod record_order;
pub(crate) mod region;
pub(crate) mod resolution_record;
pub(crate) mod resolve;
mod resolve_record;
pub(crate) mod shadow;
pub(crate) mod table;

#[cfg(test)]
mod branch_tests;
#[cfg(test)]
mod consequence_tests;
#[cfg(test)]
mod face_info_tests;
#[cfg(test)]
mod face_path_tests;
#[cfg(test)]
mod face_side_tests;
#[cfg(test)]
mod incidence_tests;
#[cfg(test)]
mod partition_tests;
#[cfg(test)]
mod pave_block_tests;
#[cfg(test)]
mod record_order_tests;
#[cfg(test)]
mod tests;
