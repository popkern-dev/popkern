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

/// Semantic ownership for a sub-face that represents one piece of a
/// larger multi-sided macro-cap construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CornerSubPatchProvenance {
    /// Multi-sided macro-cap construction family that produced the sub-patch.
    pub family: CornerSubPatchFamily,
    /// Original polar sector index before macro subdivision.
    pub sector_index: usize,
    /// Number of original polar sectors at the corner.
    pub sector_count: usize,
    /// Piece index within the sector-local macro layout.
    pub piece_index: usize,
    /// Number of pieces per sector.
    pub piece_count: usize,
    /// CAD-domain role of this piece in the macro layout.
    pub role: CornerSubPatchRole,
}

/// Multi-sided macro-cap family that produced a sector sub-patch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CornerSubPatchFamily {
    /// KP13 bi-5 macro cap.
    Kp13Bi5,
    /// KP20 bi-6 macro cap.
    Kp20Bi6,
}

/// CAD-domain role of one macro-cap sub-patch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CornerSubPatchRole {
    /// KP20 piece closer to the extraordinary point, left of the sector median.
    NearOriginLeft,
    /// KP20 piece closer to the extraordinary point, right of the sector median.
    NearOriginRight,
    /// KP20 piece closer to the outer midpoint, right of the sector median.
    NearMidpointRight,
    /// KP20 piece closer to the outer midpoint, left of the sector median.
    NearMidpointLeft,
}

impl CornerSubPatchProvenance {
    /// Build provenance for KP20's 2×2 bi-6 macro-sector layout.
    pub fn macro_cap(
        family: CornerSubPatchFamily,
        sector_index: usize,
        sector_count: usize,
        piece_index: usize,
    ) -> Self {
        let role = match piece_index {
            0 => CornerSubPatchRole::NearOriginLeft,
            1 => CornerSubPatchRole::NearOriginRight,
            2 => CornerSubPatchRole::NearMidpointRight,
            3 => CornerSubPatchRole::NearMidpointLeft,
            _ => panic!("KP20 macro piece index must be in 0..4, got {piece_index}"),
        };
        Self {
            family,
            sector_index,
            sector_count,
            piece_index,
            piece_count: 4,
            role,
        }
    }
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
    /// Macro-cap ownership metadata for multi-piece sector constructions.
    pub sub_patch_provenance: Option<CornerSubPatchProvenance>,
}
