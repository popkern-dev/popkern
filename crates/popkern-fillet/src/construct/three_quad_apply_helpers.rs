//! Small helpers for three-quad cap application.

use crate::construct::promotion_ladder::PromotionResult;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{
    ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchFamily, CornerSubPatchProvenance,
};

pub(crate) struct ThreeQuadApplyResult {
    pub patches: Vec<CornerPatchInfo>,
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

/// §5.3 item #10 — why the three-quad dispatch failed. Callers
/// use this to distinguish an honest topology / geometry failure
/// (caller should record the reason + fall through to polar) from
/// a self-intersection remediation (caller must ALSO emit the
/// `CornerRemediationOverlap` warning).
#[derive(Debug)]
pub(crate) enum ThreeQuadFailure {
    /// The Gregory sub-patches self-intersect; remediate to polar.
    #[allow(dead_code)]
    Overlap,
    /// Anything else — boundary shape, decomposition rejection,
    /// poke topology failure. Error string included for logging.
    Other(#[allow(dead_code)] String),
}

impl From<String> for ThreeQuadFailure {
    fn from(s: String) -> Self {
        ThreeQuadFailure::Other(s)
    }
}

pub(super) fn record_support_cell_nonselection(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    reason_code: &str,
    n_sides: usize,
    patch_count: usize,
) {
    diagnostics.push(
        crate::construct::three_quad_kp_diagnostics::support_cell_candidate_nonselection(
            reason_code,
            n_sides,
            patch_count,
        ),
    );
}

pub(crate) fn author_cap_face_boundary_pcurves(
    brep: &mut BRep,
    sub_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    subface_index: usize,
    old_surface: Option<&SurfaceGeometry>,
) -> Result<(), String> {
    super::three_quad_reproject::reproject_loop_pcurves_onto_gregory(
        brep,
        sub_face,
        old_surface,
        false,
    );
    let _diagnostics = super::three_quad_seam_pcurves::author_boundary_owned_seam_pcurves(
        brep,
        sub_face,
        boundary,
        fillet_surfaces,
        subface_index,
        None,
    );
    super::three_quad_reproject::author_missing_outer_loop_pcurves_from_vertices(brep, sub_face);
    super::three_quad_reproject::verify_outer_loop_pcurves(brep, sub_face)
}

pub(super) fn sub_patch_provenance(
    macro_family: Option<CornerSubPatchFamily>,
    n_sides: usize,
    sub_patch_index: usize,
) -> Option<CornerSubPatchProvenance> {
    macro_family.map(|family| {
        CornerSubPatchProvenance::macro_cap(
            family,
            sub_patch_index / 4,
            n_sides,
            sub_patch_index % 4,
        )
    })
}

pub(super) fn macro_family_from_promotion(
    promotion_result: PromotionResult,
) -> Option<CornerSubPatchFamily> {
    match promotion_result {
        PromotionResult::Standard
        | PromotionResult::ExpandedKp20Support
        | PromotionResult::Kp20SupportCell => None,
        PromotionResult::Macro4n(family) => Some(family),
    }
}

pub(super) fn decomp_self_intersects(
    decomp: &crate::blend::gregory_three_quad::NQuadDecomposition,
    tolerance: f64,
) -> bool {
    decomp
        .sub_patches
        .iter()
        .any(|patch| patch.self_intersects_with_tolerance(tolerance))
}

pub(super) fn max_outer_boundary_residual(
    decomp: &crate::blend::gregory_three_quad::NQuadDecomposition,
    boundary: &BlendBoundary,
) -> Option<f64> {
    use popkern_core::geometry::curve::CurveEval;
    if boundary.segments.len() != 2 * decomp.sub_patches.len() {
        return None;
    }
    let mut worst = 0.0_f64;
    for (k, patch) in decomp.sub_patches.iter().enumerate() {
        for (segment_index, side) in [(2 * k, 0usize), (2 * k + 1, 1usize)] {
            let curve = match boundary.segments.get(segment_index)? {
                BoundarySegment::FilletEnd { curve, .. }
                | BoundarySegment::FaceEdge { curve, .. } => curve,
            };
            let (t0, t1) = curve.parameter_range();
            for sample in 0..=10 {
                let s = sample as f64 / 10.0;
                let reference = curve.point_at(t0 + s * (t1 - t0));
                let surface_point = match side {
                    0 => patch.point_at(0.0, s),
                    _ => patch.point_at(s, 1.0),
                };
                worst = worst.max(reference.distance_to(&surface_point));
            }
        }
    }
    Some(worst)
}
