//! Post-processing for boolean results: redundant edge removal and face simplification.

pub mod collinear_vertex;
pub mod exact_rims;
pub mod face_simplify;
pub mod loop_repair;
pub mod parallel_weld;
pub mod point_loop_retire;
pub mod quality;
pub(crate) mod redundant_chain;
pub mod redundant_edge;
pub mod self_twin_stub;
pub mod slit_retire;
pub mod slit_seam;
pub mod sliver_removal;
pub mod spur_edge;

#[cfg(test)]
#[path = "postprocess_tests.rs"]
mod tests;
