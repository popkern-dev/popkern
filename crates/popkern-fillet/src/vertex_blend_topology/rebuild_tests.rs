//! End-to-end Phase 6c-C orchestration tests.
//!
//! `cube_uniform_apex_rebuild_e2e` exercises the full
//! `fillet_edges` pipeline with `POPKERN_APEX_TOPOLOGY_REBUILD=1`,
//! verifies the cap face's outer loop drops from 6 edges to 3
//! edges (the canonical sphere-fit triangle), and measures the
//! cap-vs-cylinder G0 / G1 deviation post-rebuild.
//!
//! Counterpart `cube_uniform_apex_no_rebuild_baseline` measures
//! the same metrics WITHOUT the env-var, so the deltas can be
//! reported as a regression-guard.
//!
//! Marked `#[ignore]` because they set process-level env vars;
//! run with `--ignored`.

#[path = "rebuild_tests_cube_mixed.rs"]
mod tests_cube_mixed;
#[path = "rebuild_tests_cube_uniform.rs"]
mod tests_cube_uniform;
#[path = "rebuild_tests_pyramid_diagnostics.rs"]
mod tests_pyramid_diagnostics;
#[path = "rebuild_tests_reroute_gates.rs"]
mod tests_reroute_gates;
#[path = "rebuild_tests_support.rs"]
mod tests_support;
