//! Repair operations for BRep healing.
//!
//! Provides gap repair, tolerance management, surface simplification,
//! degeneracy repair, coplanar face merging, boundary stitching,
//! self-intersection repair, hole filling, a multi-step repair pipeline,
//! and a preview-before-apply workflow.
//!
//! # Submodules
//!
//! - [`gap_repair`] -- close geometric gaps between adjacent faces
//! - [`tolerance_mgmt`] -- adjust and enforce entity tolerance hierarchy
//! - [`simplification`] -- detect and replace over-parameterized NURBS
//! - [`degeneracy`] -- remove zero-length edges and collapse small faces
//! - [`merge_faces`] -- detect and merge adjacent coplanar faces
//! - [`stitch`] -- identify and merge matching boundary edges
//! - [`self_intersection`] -- detect and repair self-intersecting surfaces
//! - [`hole_fill`] -- detect and fill holes in shell boundaries
//! - [`complexity`] -- reduce unnecessary NURBS knots
//! - [`rederive`] -- recompute edge curves from adjacent surface intersections
//! - [`pipeline`] -- orchestrate repair passes with transactional rollback
//! - [`preview`] -- preview proposed repairs before committing

/// NURBS complexity reduction (knot removal).
pub mod complexity;
/// Coons patch surface construction for hole filling.
pub mod coons_patch;
/// Degeneracy repair (zero-length edges, small faces).
pub mod degeneracy;
/// Analytical surface detection (plane, cylinder, sphere, cone, torus).
pub mod detect_analytical;
/// G1-continuous Coons patch with cross-boundary derivative matching.
pub mod g1_coons_patch;
/// Gap repair strategies and execution.
pub mod gap_repair;
/// Hole filling (wound healing) for shell boundaries.
pub mod hole_fill;
/// Import repair pipeline for models from external formats.
pub mod import_repair;
/// Coplanar face merging.
pub mod merge_faces;
/// N-sided hole filling by decomposition into 4-sided patches.
pub mod n_sided_fill;
pub mod n_sided_fill_g1;
/// Multi-step repair pipeline with transactional rollback.
pub mod pipeline;
/// Preview proposed repairs before applying.
pub mod preview;
/// Progressive healing with checkpoints and severity-based phases.
pub mod progressive;
/// Re-derivation gap repair (recompute edge curves from adjacent surfaces).
pub mod rederive;
/// Arc-length reparameterization for NURBS UV uniformity.
pub mod reparameterize;
/// Self-intersection detection and repair.
pub mod self_intersection;
/// Surface simplification (NURBS to analytical).
pub mod simplification;
/// Spatial grid for proximity acceleration.
pub(crate) mod spatial_grid;
/// Boundary assembly (matching and merging free boundary edges).
pub mod stitch;
/// Tolerance hierarchy enforcement and auto-adjustment.
pub mod tolerance_mgmt;
/// Geometry-type-specific healing for fillet, offset, and swept surfaces.
pub mod type_specific;

pub use complexity::{ComplexityReductionResult, CurveReductionResult, SurfaceReductionResult};
pub use coons_patch::{CoonsBoundary, CoonsPatchResult};
pub use degeneracy::{
    DegeneracyRepairResult, detect_cusps, detect_self_tangent_faces, detect_zero_area_faces,
    repair_degeneracies_full,
};
pub use g1_coons_patch::{CrossDerivatives, SecondDerivatives};
pub use gap_repair::{GapRepairResult, GapRepairStrategy};
pub use hole_fill::{DetectedHole, HoleFillResult, HoleFillStrategy, fill_by_extension};
pub use import_repair::{
    CustomRepairOp, ImportFormat, ImportRepairConfig, ImportRepairOp, ImportRepairReport,
    RepairOperation, StepLabel,
};
pub use merge_faces::FaceMergeResult;
pub use n_sided_fill::fill_n_sided_hole;
pub use n_sided_fill_g1::{G1HoleFillError, fill_n_sided_hole_g1};
pub use pipeline::{ConvergenceStatus, HealResult};
pub use preview::RepairPreview;
pub use progressive::{HealingCheckpoint, ProgressiveHealResult, heal_progressive};
pub use rederive::RederiveResult;
pub use reparameterize::ReparamResult;
pub use self_intersection::{SelfIntersectionRepairResult, repair_by_splitting};
pub use simplification::{AnalyticalKind, SimplificationResult};
pub use stitch::StitchResult;
pub use tolerance_mgmt::ToleranceAdjustment;
pub use type_specific::{SurfaceType, TypeDetectionResult, TypeSpecificRepairResult};
