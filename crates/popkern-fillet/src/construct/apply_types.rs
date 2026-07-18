//! Public result types for blend topology construction.
//!
//! Extracted from `construct::mod` (file-size cap): the per-edge
//! computed geometry bundle and the apply-phase result envelope.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle};

use crate::analysis::EdgeAnalysis;
use crate::error::BlendWarning;

/// Geometry computed for one fillet/chamfer edge.
#[derive(Debug)]
pub struct EdgeGeometry {
    /// The edge analysis.
    pub analysis: EdgeAnalysis,
    /// Fillet or chamfer surface.
    pub surface: SurfaceGeometry,
    /// Trim curve on face 1 (3D).
    pub trim_3d_1: NurbsCurve,
    /// Trim curve on face 2 (3D).
    pub trim_3d_2: NurbsCurve,
    /// Exact trim circle on face 1 `(center, normal, radius)` when
    /// the contact is circular (see `trim::TrimCurves::circle_1`).
    pub trim_circle_1: Option<(
        popkern_core::geometry::Point3,
        popkern_core::geometry::Vec3,
        f64,
    )>,
    /// Exact trim circle on face 2; see `trim_circle_1`.
    pub trim_circle_2: Option<(
        popkern_core::geometry::Point3,
        popkern_core::geometry::Vec3,
        f64,
    )>,
    /// Pcurve on face 1.
    pub pcurve_1: NurbsCurve,
    /// Pcurve on face 2.
    pub pcurve_2: NurbsCurve,
    /// Average fillet/chamfer radius (kept for legacy callers).
    pub effective_radius: f64,
    /// Radius at the edge's start vertex (constant ⇒ `effective_radius`; variable ⇒ `start`).
    pub radius_at_start: f64,
    /// Radius at the edge's end vertex.
    pub radius_at_end: f64,
    /// §18 telemetry: walker-path Newton iterations.
    pub iterations_used: Option<usize>,
    /// §18 telemetry: line-search rejections.
    pub step_rejections_used: Option<usize>,
    /// §8.3 start-endpoint cap (`None` when `cap_kind == Open`).
    pub start_cap: Option<crate::cap::CapGeometry>,
    /// §8.3 end-endpoint cap (`None` when `cap_kind == Open`).
    pub end_cap: Option<crate::cap::CapGeometry>,
    /// §7 continuity downgrade record — `Some((achieved, reason))` when the engine fell back; caller emits [`BlendWarning::ContinuityDowngraded`](crate::error::BlendWarning::ContinuityDowngraded).
    pub continuity_downgrade: Option<(crate::types::ContinuityLevel, String)>,
}

/// Result of applying blend geometry, including warnings.
pub struct ApplyBlendResult {
    /// The modified BRep.
    pub brep: BRep,
    /// Warnings generated during vertex mode handling.
    pub warnings: Vec<(EdgeHandle, BlendWarning)>,
    /// §5.3 per-corner patch reports — one entry per vertex where a blend patch was inserted.
    pub corner_patches: Vec<crate::types::CornerPatchInfo>,
    /// Structured diagnostics for repair and hardening actions.
    pub repair_report: crate::types::BlendRepairReport,
}
