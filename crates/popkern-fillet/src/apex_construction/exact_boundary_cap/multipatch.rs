//! Candidate-only multi-patch exact-boundary cap construction.
//!
//! This module starts as a dormant implementation surface. Runtime routing
//! remains unchanged until the orchestration and transaction slices are built.

#[path = "multipatch_budget.rs"]
pub(crate) mod budget;

#[path = "multipatch_commit.rs"]
mod commit;

#[path = "multipatch_inner.rs"]
mod inner;

#[path = "multipatch_inventory.rs"]
mod inventory;

#[path = "multipatch_report.rs"]
mod report;

#[path = "multipatch_strip.rs"]
mod strip;

#[path = "multipatch_validation.rs"]
mod validation;
