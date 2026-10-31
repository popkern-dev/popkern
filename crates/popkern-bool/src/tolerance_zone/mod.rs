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
//! **Stage.** The resolution runs in shadow beside the intersection
//! graph's own node identity (`IntersectionGraph::find_node`: first
//! match within 2 × linear, record ends adopted within 20 × linear).
//! `shadow` compares the two partitions on every cellular boolean when
//! `POPKERN_COMPLEX_SHADOW=1`, so the table can replace graph node
//! identity once the disagreements are explained (plan Stage 2 gate).
//!
//! **Tolerance.** Every radius is an entity's own tolerance, the model
//! linear tolerance when it has none; the growth cap is an explicit
//! parameter. No threshold here is global state.

pub(crate) mod branches;
pub(crate) mod candidates;
pub(crate) mod common_blocks;
pub(crate) mod consequences;
pub(crate) mod face_path;
pub(crate) mod face_side;
pub(crate) mod incidence;
pub(crate) mod partition;
pub(crate) mod pave_blocks;
pub(crate) mod record_order;
pub(crate) mod region;
pub(crate) mod resolve;
pub(crate) mod shadow;
pub(crate) mod table;

#[cfg(test)]
mod branch_tests;
#[cfg(test)]
mod consequence_tests;
#[cfg(test)]
mod face_path_tests;
#[cfg(test)]
mod face_side_tests;
#[cfg(test)]
mod partition_tests;
#[cfg(test)]
mod pave_block_tests;
#[cfg(test)]
mod record_order_tests;
#[cfg(test)]
mod tests;
