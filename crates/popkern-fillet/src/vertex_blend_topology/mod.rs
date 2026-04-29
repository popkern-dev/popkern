//! Phase 6c-C — apex-region topology rebuild for sphere-fit caps.
//!
//! See `docs/decisions/2026-04-28-phase-6c-c-topology-design.md`
//! for the architectural rationale and implementation phases.
//!
//! Sub-phases:
//! * **6c-C.1 — detection** (`detect`): given an apex vertex, find
//!   the BRep sub-graph that needs to be replaced.
//! * **6c-C.2 — deletion** (`delete`, future): remove the detected
//!   sub-graph via Euler ops.
//! * **6c-C.3 — stitching** (`stitch`, future): insert the
//!   sphere-fit topology in its place.
//!
//! All entry points are `pub(crate)`. Public-API exposure is
//! through `popkern_fillet::vertex_blend::compute_vertex_blend_terminations`
//! and the existing `fillet_edges` pipeline (Phase 6c-C.4).

pub(crate) mod collapse;
pub(crate) mod detect;
pub(crate) mod execute;
pub(crate) mod plan;
pub(crate) mod rebuild;
pub(crate) mod var_radius;

#[cfg(test)]
#[path = "detect_tests.rs"]
mod detect_tests;

#[cfg(test)]
#[path = "plan_tests.rs"]
mod plan_tests;

#[cfg(test)]
#[path = "execute_tests.rs"]
mod execute_tests;

#[cfg(test)]
#[path = "collapse_tests.rs"]
mod collapse_tests;

#[cfg(test)]
#[path = "rebuild_tests.rs"]
mod rebuild_tests;
