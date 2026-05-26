//! Small helpers for three-quad cap application.

use crate::construct::promotion_ladder::PromotionResult;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use crate::blend::boundary::BlendBoundary;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{
    ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchFamily, CornerSubPatchProvenance,
};

pub(crate) struct ThreeQuadApplyResult {
    pub patches: Vec<CornerPatchInfo>,
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
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
        PromotionResult::Standard => None,
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
