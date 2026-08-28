//! Public result types for blend topology construction.
//!
//! Extracted from `construct::mod` (file-size cap): the per-edge
//! computed geometry bundle and the apply-phase result envelope.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle};

use crate::analysis::EdgeAnalysis;
use crate::error::{BlendFailure, BlendWarning};

/// Failure modes of [`crate::construct::apply_blend_geometry`].
///
/// Both variants are transactional: the caller's input BRep is
/// untouched and the working clone is discarded.
#[derive(Debug, thiserror::Error)]
pub enum ApplyBlendError {
    /// Per-edge construction failures — the whole apply declines.
    #[error("{} edge(s) failed during blend apply", .0.len())]
    Failures(Vec<(EdgeHandle, BlendFailure)>),
    /// A cancellation token fired between apply stages (§17 clean
    /// termination — polls sit between phases and between corner
    /// constructions, bounding the post-cancel stretch to one corner).
    #[error("blend apply cancelled by token")]
    Cancelled,
}

impl ApplyBlendError {
    /// Adapt to the pipeline-level error: `Cancelled` maps straight
    /// through; `Failures` becomes the transactional `PartialFailure`
    /// carrying an unmodified copy of the input BRep.
    pub fn into_blend_error(self, total: usize, brep: &BRep) -> crate::error::BlendError {
        match self {
            Self::Cancelled => crate::error::BlendError::Cancelled,
            Self::Failures(failures) => crate::error::BlendError::PartialFailure {
                failed: failures.len(),
                total,
                partial_result: Box::new(brep.clone()),
                failures,
            },
        }
    }
}

/// Emit one PK-276 per-corner progress event and advance the corner
/// cursor. Factored out of the vertex-blend loop (two exits per
/// corner: direct-apex success `continue` and the vertex-blend tail).
pub(crate) fn emit_corner_progress(
    reporter: &mut Option<&mut (dyn crate::progress::ProgressReporter + '_)>,
    vertex: popkern_core::topology::VertexHandle,
    corner_index: &mut usize,
    total_corners: usize,
    started: std::time::Instant,
) {
    if let Some(r) = reporter.as_deref_mut() {
        r.on_corner(crate::progress::CornerProgress {
            vertex,
            corner_index: *corner_index,
            total_corners,
            duration: started.elapsed(),
        });
    }
    *corner_index += 1;
}

/// Poll a cancellation token between apply stages: `Err(Cancelled)`
/// the moment any holder has fired it, `Ok(())` otherwise (including
/// when no token was supplied).
pub fn apply_cancel_poll(
    cancel: Option<&crate::cancellation::CancellationToken>,
) -> Result<(), ApplyBlendError> {
    match cancel {
        Some(token) if token.is_cancelled() => Err(ApplyBlendError::Cancelled),
        _ => Ok(()),
    }
}

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
    /// Average fillet/chamfer radius (kept for existing callers).
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
