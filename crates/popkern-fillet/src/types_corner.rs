use popkern_core::topology::{FaceHandle, VertexHandle};

use crate::types::ContinuityLevel;

/// §5.3 corner-resolution method used to build a single vertex-
/// blend patch. Reported on [`CornerPatchInfo::method`] so a UI
/// can tell at a glance which algorithm fired and colour the
/// corners accordingly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CornerMethod {
    /// Rolling-ball sphere — equal-radius cube-corner case. G1.
    RollingBallSphere,
    /// Polar NURBS surface fit through the boundary curves. G0.
    PolarNurbs,
    /// Flat (mitred) plane fit — legacy G0 fallback.
    Mitred,
    /// §5.3 Coons patch for N=4 corners. G1. [Pending — C3.]
    Coons,
    /// §5.3 Gregory rational patch for general N. G1; G2 for N=3, 4.
    Gregory,
    /// Polyhedral-Net Splines (KK-Peters 2015) bi-4 cap. G1 with
    /// the supporting fillet/face surfaces at the boundary; valid
    /// for N in {3, 5, 6, 7, 8, 9, 10, 11, 12}.
    PolyhedralNetSpline,
    /// §5.3 item #7 — planar cap at a tangent corner.
    TangentPlane,
    /// §5.3 Plate-minimisation (thin-plate energy). [BEYOND.]
    PlateMin,
}

impl CornerMethod {
    /// Display string for diagnostics and UI.
    pub fn as_str(&self) -> &'static str {
        match self {
            CornerMethod::RollingBallSphere => "RollingBallSphere",
            CornerMethod::PolarNurbs => "PolarNurbs",
            CornerMethod::Mitred => "Mitred",
            CornerMethod::Coons => "Coons",
            CornerMethod::Gregory => "Gregory",
            CornerMethod::PolyhedralNetSpline => "PolyhedralNetSpline",
            CornerMethod::TangentPlane => "TangentPlane",
            CornerMethod::PlateMin => "PlateMin",
        }
    }
}

impl std::fmt::Display for CornerMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// §5.3 per-corner patch report emitted on
/// [`crate::types::BlendResult::corner_patches`].
#[derive(Debug, Clone)]
pub struct CornerPatchInfo {
    /// The original BRep vertex the corner was resolved at.
    pub vertex: VertexHandle,
    /// Number of incoming blend boundaries.
    pub n_sides: usize,
    /// Which construction method produced the patch.
    pub method: CornerMethod,
    /// The continuity level actually achieved at the corner boundary.
    pub achieved_continuity: ContinuityLevel,
    /// §5.3 C6 — the vertex-blend classification that drove dispatch.
    pub classification: crate::vertex_classifier::VertexBlendClassification,
    /// The corner face handle in the output BRep.
    pub face: FaceHandle,
    /// Surface-variant names for blend-bridge fallbacks.
    pub bridge_fallbacks: Vec<&'static str>,
    /// True when a self-intersecting corner was remediated.
    pub overlap_remediated: bool,
    /// Worst coplanarity when K-P binonic relaxed the continuity gate.
    pub relaxed_continuity_coplanarity: Option<f64>,
    /// Per-fillet sphere fit `(centre, radius)` for this sub-face.
    pub cap_sphere: Option<(popkern_core::geometry::Point3, f64)>,
}
