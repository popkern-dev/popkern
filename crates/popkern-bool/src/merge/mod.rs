//! Stage 5: Merge selected face fragments into a single result BRep.
//!
//! Uses the intersection graph to walk face boundaries, classify fragments,
//! and assemble a valid closed BRep from the selected face fragments of
//! two split operands.

pub(crate) mod assemble;
pub(crate) mod curve_split;
pub(crate) mod fragment_classify;
mod fragment_classify_overlap;
pub(crate) mod fragment_sample;
pub(crate) mod fragment_sample_interior;
pub(crate) mod graph;
pub(crate) mod graph_build;
pub(crate) mod graph_intersect;
pub(crate) mod graph_merge;
pub(crate) mod graph_split;
pub(crate) mod graph_split_node;
pub(crate) mod graph_split_pcurve;
pub(crate) mod graph_split_shared;
pub(crate) mod graph_weld;
pub(crate) mod manifold_check;
pub(crate) mod radial;
pub(crate) mod walk;
