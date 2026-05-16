//! Phase 6c-A2 — Direct sphere-fit apex construction.
//!
//! The fillet pipeline's apex resolution path: plan apex topology
//! (with sphere-fit terminations) first, then build the BRep apex
//! region directly via Euler operators.  Replaces the legacy
//! substitute-then-patch overlay (Phases 6c-A through 6c-X).
//!
//! See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-construction.md`.
//!
//! # Public API
//!
//! Three entry points (per ADR §5.1):
//!
//! * [`plan_apex`] — pure planner, returns an [`ApexPlan`] without
//!   mutating the BRep.  UI-friendly: callers can preview, inspect
//!   telemetry, or cancel.
//! * [`construct_apex`] — applies a precomputed [`ApexPlan`] to a
//!   BRep via Euler ops.  Snapshot-rollback on failure.
//! * [`resolve_apex`] — combined `plan_apex` + `construct_apex` for
//!   callers who don't need to inspect the plan first.
//!
//! All three accept a [`Tolerance`](popkern_core::tolerance::Tolerance)
//! and an optional [`CancelToken`](crate::cancellation::CancelToken).
//!
//! # Status
//!
//! A2.1 in progress: type API committed; planner under
//! implementation.  A2.2 (builder), A2.3 (pipeline integration),
//! A2.4 (cap surface), A2.5 (validation), A2.6 (cleanup) follow.

pub(crate) mod apex_vertex_authoring;
pub(crate) mod apex_vertex_removal;
pub(crate) mod cap_kind;
pub(crate) mod cap_loop_construction;
pub(crate) mod cap_loop_construction_options;
pub(crate) mod cap_surface_assignment;
pub(crate) mod construct;
pub(crate) mod construct_cleanup;
pub(crate) mod construct_edges;
pub(crate) mod construct_result;
pub(crate) mod construct_surface;
pub(crate) mod construct_vertices;
pub(crate) mod constructive_builder;
pub(crate) mod constructive_builder_options;
pub(crate) mod constructive_geometry;
pub(crate) mod constructive_topology;
pub(crate) mod cylinder_face_coherence;
pub(crate) mod error;
pub(crate) mod face_edge_construction;
pub(crate) mod gregory_cap;
pub(crate) mod gregory_cap_boundary;
pub(crate) mod gregory_cap_topology;
pub(crate) mod perp_anchor_curve_refit;
pub(crate) mod perp_anchor_relocate;
pub(crate) mod plan;
pub(crate) mod plan_rung_fallback;
pub(crate) mod plan_terminations;
pub(crate) mod plan_topology;
pub(crate) mod plan_validate;
pub(crate) mod planner_prepass;
pub(crate) mod repeated_apex_loop_plan;
pub(crate) mod repeated_apex_partition;
pub(crate) mod repeated_apex_partition_topology;
pub(crate) mod repeated_apex_region_candidates;
pub(crate) mod repeated_apex_region_graph;
pub(crate) mod repeated_apex_serial;
pub(crate) mod repeated_apex_side_support;
pub(crate) mod repeated_apex_side_support_boundary;
pub(crate) mod repeated_apex_side_support_bridge;
pub(crate) mod repeated_apex_side_support_bridge_operator;
pub(crate) mod repeated_apex_side_support_execution;
pub(crate) mod repeated_apex_side_support_facts;
pub(crate) mod repeated_apex_side_support_live_preflight;
pub(crate) mod repeated_apex_side_support_merge;
pub(crate) mod repeated_apex_side_support_merge_geometry;
pub(crate) mod repeated_apex_side_support_outside;
pub(crate) mod repeated_apex_side_support_outside_core;
pub(crate) mod repeated_apex_side_support_outside_operator;
pub(crate) mod repeated_apex_side_support_plan;
pub(crate) mod repeated_apex_side_support_plan_access;
pub(crate) mod repeated_apex_side_support_result;
pub(crate) mod repeated_apex_topology_knowledge;
pub(crate) mod repeated_apex_topology_knowledge_collect;
pub(crate) mod repeated_apex_topology_knowledge_trim;
pub(crate) mod repeated_apex_trim_subcomplex;
pub(crate) mod target_topology;
pub(crate) mod telemetry;

#[cfg(test)]
#[path = "plan_tests.rs"]
mod plan_tests;

#[cfg(test)]
#[path = "plan_topology_tests.rs"]
mod plan_topology_tests;

#[cfg(test)]
#[path = "plan_e2e_test_helpers.rs"]
mod plan_e2e_test_helpers;

#[cfg(test)]
#[path = "plan_e2e_tests.rs"]
mod plan_e2e_tests;

#[cfg(test)]
#[path = "plan_e2e_extras_tests.rs"]
mod plan_e2e_extras_tests;

#[cfg(test)]
#[path = "construct_tests.rs"]
mod construct_tests;

#[cfg(test)]
#[path = "construct_test_helpers.rs"]
mod construct_test_helpers;

#[cfg(test)]
#[path = "repeated_apex_side_support_exact_tests.rs"]
mod repeated_apex_side_support_exact_tests;

#[cfg(test)]
#[path = "repeated_apex_migration_parity_tests.rs"]
mod repeated_apex_migration_parity_tests;

pub use cap_kind::CapKind;
pub use construct::{construct_apex, construct_apex_with_surfaces};
pub use construct_result::ApexBuildResult;
pub use error::{ApexBuildError, ApexError, ApexPlanError};
pub use plan::{ApexPlan, plan_apex_resolution};
pub use target_topology::{TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology};
pub use telemetry::{ConstructionRung, PlanTelemetry};
