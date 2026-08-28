//! Apex-region topology inspection.
//!
//! Historical context: this module once carried the full route-B
//! apex-region rebuild stack (plan → mef-execute → lens-collapse,
//! see `docs/decisions/2026-04-28-phase-6c-c-topology-design.md`).
//! Route-B corner reconstruction replaced the anchor-hexagon cap
//! anatomy that machinery was built to clean up, leaving it
//! unreachable from any production path; it was retired with
//! PK-181. `detect` stays: its tests pin the reconstructed
//! spherical-triangle cap anatomy.
//!
//! All entry points are `pub(crate)`. Public-API exposure is
//! through `popkern_fillet::vertex_blend::compute_vertex_blend_terminations`
//! and the existing `fillet_edges` pipeline.

pub(crate) mod detect;

#[cfg(test)]
#[path = "detect_tests.rs"]
mod detect_tests;
