//! Repair operations for BRep healing.
//!
//! Provides gap repair, tolerance management, surface simplification,
//! a multi-step repair pipeline, and a preview-before-apply workflow.
//!
//! # Submodules
//!
//! - [`gap_repair`] -- close geometric gaps between adjacent faces
//! - [`tolerance_mgmt`] -- adjust and enforce entity tolerance hierarchy
//! - [`simplification`] -- detect and replace over-parameterized NURBS
//! - [`pipeline`] -- orchestrate repair passes with transactional rollback
//! - [`preview`] -- preview proposed repairs before committing

/// Gap repair strategies and execution.
pub mod gap_repair;
/// Multi-step repair pipeline with transactional rollback.
pub mod pipeline;
/// Preview proposed repairs before applying.
pub mod preview;
/// Surface simplification (NURBS to analytical).
pub mod simplification;
/// Tolerance hierarchy enforcement and auto-adjustment.
pub mod tolerance_mgmt;

pub use gap_repair::{GapRepairResult, GapRepairStrategy};
pub use pipeline::HealResult;
pub use preview::RepairPreview;
pub use simplification::{AnalyticalKind, SimplificationResult};
pub use tolerance_mgmt::ToleranceAdjustment;
