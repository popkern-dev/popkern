//! Disabled-by-default PnS cap candidate for constructive DirectApex Gregory caps.

use popkern_core::error::ValidationError;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend::boundary::BlendBoundary;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo};
use crate::vertex_classifier::VertexBlendClassification;

use super::gregory_cap::GregoryCapSector;

pub(super) struct BoundaryPreservingPnsCapInput<'a> {
    pub brep: &'a mut BRep,
    pub cap_face: FaceHandle,
    pub boundary: &'a BlendBoundary,
    pub fillet_surfaces: &'a [FilletSurfaceResult],
    pub vertex: VertexHandle,
    pub n_sides: usize,
    pub classification: VertexBlendClassification,
    pub sectors: &'a [GregoryCapSector],
}

pub(super) fn try_boundary_preserving_pns_cap(
    input: BoundaryPreservingPnsCapInput<'_>,
) -> Option<Vec<CornerPatchInfo>> {
    let BoundaryPreservingPnsCapInput {
        brep,
        cap_face,
        boundary,
        fillet_surfaces,
        vertex,
        n_sides,
        classification,
        sectors,
    } = input;
    let enabled = std::env::var("POPKERN_DIRECT_APEX_PNS_CANDIDATE")
        .map(|value| value == "1")
        .unwrap_or(false);
    if !enabled {
        return None;
    }

    let snapshot = brep.clone();
    let pns = match crate::blend::pns_promote::insert_pns_cap_into_blend_face(
        brep,
        cap_face,
        boundary,
        fillet_surfaces,
    ) {
        Ok(pns) => pns,
        Err(_) => {
            *brep = snapshot;
            return None;
        }
    };

    let face_edge_handles = sectors
        .iter()
        .map(|sector| sector.face_edge)
        .collect::<Vec<_>>();
    super::cap_surface_assignment::author_plane_side_face_edge_pcurves(brep, &face_edge_handles);
    for (index, sub_face) in pns.sub_faces.iter().copied().enumerate() {
        if crate::construct::three_quad_apply::author_cap_face_boundary_pcurves(
            brep,
            sub_face,
            boundary,
            fillet_surfaces,
            index,
            None,
        )
        .is_err()
        {
            *brep = snapshot;
            return None;
        }
    }

    let patches = pns
        .sub_faces
        .iter()
        .copied()
        .map(|face| CornerPatchInfo {
            vertex,
            n_sides,
            method: CornerMethod::PolyhedralNetSpline,
            achieved_continuity: ContinuityLevel::G1,
            classification,
            face,
            bridge_fallbacks: boundary.bridge_fallbacks.clone(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: None,
            cap_sphere: None,
            sub_patch_provenance: None,
            setback_cap: false,
            exact_boundary_cap: false,
        })
        .collect::<Vec<_>>();

    let topology_error = brep.validate_collect_fast().into_iter().find(|err| {
        matches!(
            err,
            ValidationError::BrokenCycle { .. }
                | ValidationError::ShellGap { .. }
                | ValidationError::OpenLoop { .. }
        )
    });
    let cap_faces = patches.iter().map(|patch| patch.face).collect::<Vec<_>>();
    let continuity_error = crate::construct::cap_outer_acceptance::cap_support_acceptance_error(
        brep, &cap_faces, 5.0e-3, 5.0e-2,
    );
    if topology_error.is_some() || continuity_error.is_some() {
        *brep = snapshot;
        return None;
    }

    Some(patches)
}
