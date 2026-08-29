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
    /// Flat (mitred) plane fit — the mitred G0 fallback.
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
    /// Boundary-side lane emitted from a split KP20 support-owned piece.
    SupportBoundaryLane,
    /// Interior-side lane emitted from a split KP20 support-owned piece.
    SupportInteriorLane,
    /// Rectangular support-cell emitted from a KP20 support-owned cell primitive
    /// whose original near-midpoint owner is not encoded.
    SupportCell,
    /// Rectangular support-cell replacing the near-midpoint-right KP20 owner.
    SupportCellNearMidpointRight,
    /// Rectangular support-cell replacing the near-midpoint-left KP20 owner.
    SupportCellNearMidpointLeft,
}

impl CornerSubPatchRole {
    /// True for every support-cell role, including owner-specific cells.
    pub fn is_support_cell(self) -> bool {
        matches!(
            self,
            CornerSubPatchRole::SupportCell
                | CornerSubPatchRole::SupportCellNearMidpointRight
                | CornerSubPatchRole::SupportCellNearMidpointLeft
        )
    }

    /// Original macro owner replaced by an owner-specific support cell.
    pub fn support_cell_owner_role(self) -> Option<Self> {
        match self {
            CornerSubPatchRole::SupportCellNearMidpointRight => {
                Some(CornerSubPatchRole::NearMidpointRight)
            }
            CornerSubPatchRole::SupportCellNearMidpointLeft => {
                Some(CornerSubPatchRole::NearMidpointLeft)
            }
            _ => None,
        }
    }
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

    /// Build provenance for a KP20 expanded support-lane patchwork.
    pub fn kp20_expanded_support_lane(
        sector_index: usize,
        sector_count: usize,
        patch_index: usize,
        patch_count: usize,
        boundary_lane: bool,
    ) -> Self {
        Self {
            family: CornerSubPatchFamily::Kp20Bi6,
            sector_index,
            sector_count,
            piece_index: patch_index,
            piece_count: patch_count,
            role: if boundary_lane {
                CornerSubPatchRole::SupportBoundaryLane
            } else {
                CornerSubPatchRole::SupportInteriorLane
            },
        }
    }

    /// Build provenance for a KP20 support-cell primitive surface.
    pub fn kp20_support_cell(
        sector_index: usize,
        sector_count: usize,
        patch_index: usize,
        patch_count: usize,
    ) -> Self {
        Self {
            family: CornerSubPatchFamily::Kp20Bi6,
            sector_index,
            sector_count,
            piece_index: patch_index,
            piece_count: patch_count,
            role: CornerSubPatchRole::SupportCell,
        }
    }

    /// Build owner-aware provenance for a KP20 support-cell primitive surface.
    pub fn kp20_support_cell_for_owner_role(
        sector_index: usize,
        sector_count: usize,
        patch_index: usize,
        patch_count: usize,
        owner_role: CornerSubPatchRole,
    ) -> Self {
        let role = match owner_role {
            CornerSubPatchRole::NearMidpointRight => {
                CornerSubPatchRole::SupportCellNearMidpointRight
            }
            CornerSubPatchRole::NearMidpointLeft => CornerSubPatchRole::SupportCellNearMidpointLeft,
            _ => {
                panic!(
                    "KP20 support-cell owner role must be a near-midpoint support owner, got {owner_role:?}"
                )
            }
        };
        Self {
            family: CornerSubPatchFamily::Kp20Bi6,
            sector_index,
            sector_count,
            piece_index: patch_index,
            piece_count: patch_count,
            role,
        }
    }

    /// Build provenance for a macro piece fragmented by expanded KP20 support topology.
    pub fn kp20_expanded_macro_fragment(
        sector_index: usize,
        sector_count: usize,
        patch_index: usize,
        patch_count: usize,
        role: CornerSubPatchRole,
    ) -> Self {
        Self {
            family: CornerSubPatchFamily::Kp20Bi6,
            sector_index,
            sector_count,
            piece_index: patch_index,
            piece_count: patch_count,
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
    /// True when this patch was produced by the SETBACK vertex-blend
    /// surgery (PK-094/PK-097 cap families). The setback surgery
    /// owns its apex topology: the vertex-blend post passes (sliver absorb,
    /// cylinder off-surface topology repair) and vertex-blend runoff
    /// scoping key off this to keep the two worlds apart.
    pub setback_cap: bool,
    /// True when this patch is a committed sector face of an ACCEPTED
    /// exact-boundary cap candidate (one-patch Coons wedge or its
    /// multipatch partition fallback). The exact-boundary route owns
    /// its topology and pcurves — acceptance gates every cap/support
    /// pcurve pair transactionally — so vertex-blend corner-patch pcurve
    /// completion must not re-author them (same quarantine contract
    /// as `setback_cap`, PK-087/PK-119).
    pub exact_boundary_cap: bool,
    /// Macro-cap ownership metadata for multi-piece sector constructions.
    pub sub_patch_provenance: Option<CornerSubPatchProvenance>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kp20_provenance_records_sector_piece_and_role() {
        let provenance =
            CornerSubPatchProvenance::macro_cap(CornerSubPatchFamily::Kp20Bi6, 2, 7, 3);

        assert_eq!(provenance.family, CornerSubPatchFamily::Kp20Bi6);
        assert_eq!(provenance.sector_index, 2);
        assert_eq!(provenance.sector_count, 7);
        assert_eq!(provenance.piece_index, 3);
        assert_eq!(provenance.piece_count, 4);
        assert_eq!(provenance.role, CornerSubPatchRole::NearMidpointLeft);
    }

    #[test]
    fn kp13_and_kp20_macro_families_are_distinct() {
        let kp13 = CornerSubPatchProvenance::macro_cap(CornerSubPatchFamily::Kp13Bi5, 0, 5, 0);
        let kp20 = CornerSubPatchProvenance::macro_cap(CornerSubPatchFamily::Kp20Bi6, 0, 5, 0);

        assert_ne!(kp13.family, kp20.family);
        assert_eq!(kp13.role, kp20.role);
    }

    #[test]
    fn kp20_expanded_support_lane_provenance_records_lane_role() {
        let boundary = CornerSubPatchProvenance::kp20_expanded_support_lane(1, 6, 7, 30, true);
        let interior = CornerSubPatchProvenance::kp20_expanded_support_lane(1, 6, 8, 30, false);

        assert_eq!(boundary.family, CornerSubPatchFamily::Kp20Bi6);
        assert_eq!(boundary.sector_index, 1);
        assert_eq!(boundary.sector_count, 6);
        assert_eq!(boundary.piece_index, 7);
        assert_eq!(boundary.piece_count, 30);
        assert_eq!(boundary.role, CornerSubPatchRole::SupportBoundaryLane);
        assert_eq!(interior.role, CornerSubPatchRole::SupportInteriorLane);
    }

    #[test]
    fn kp20_support_cell_provenance_records_cell_role() {
        let cell = CornerSubPatchProvenance::kp20_support_cell(2, 6, 9, 24);

        assert_eq!(cell.family, CornerSubPatchFamily::Kp20Bi6);
        assert_eq!(cell.sector_index, 2);
        assert_eq!(cell.sector_count, 6);
        assert_eq!(cell.piece_index, 9);
        assert_eq!(cell.piece_count, 24);
        assert_eq!(cell.role, CornerSubPatchRole::SupportCell);
    }

    #[test]
    fn kp20_owner_support_cell_provenance_records_replaced_owner() {
        let cell = CornerSubPatchProvenance::kp20_support_cell_for_owner_role(
            3,
            6,
            10,
            24,
            CornerSubPatchRole::NearMidpointRight,
        );

        assert_eq!(cell.family, CornerSubPatchFamily::Kp20Bi6);
        assert_eq!(cell.sector_index, 3);
        assert_eq!(cell.sector_count, 6);
        assert_eq!(cell.piece_index, 10);
        assert_eq!(cell.piece_count, 24);
        assert_eq!(cell.role, CornerSubPatchRole::SupportCellNearMidpointRight);
        assert!(cell.role.is_support_cell());
        assert_eq!(
            cell.role.support_cell_owner_role(),
            Some(CornerSubPatchRole::NearMidpointRight)
        );
    }

    #[test]
    fn kp20_expanded_macro_fragment_provenance_preserves_macro_role() {
        let fragment = CornerSubPatchProvenance::kp20_expanded_macro_fragment(
            2,
            6,
            12,
            36,
            CornerSubPatchRole::NearOriginLeft,
        );

        assert_eq!(fragment.family, CornerSubPatchFamily::Kp20Bi6);
        assert_eq!(fragment.sector_index, 2);
        assert_eq!(fragment.sector_count, 6);
        assert_eq!(fragment.piece_index, 12);
        assert_eq!(fragment.piece_count, 36);
        assert_eq!(fragment.role, CornerSubPatchRole::NearOriginLeft);
    }
}
