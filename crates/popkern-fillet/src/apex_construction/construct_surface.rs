//! Cap-surface assignment for the legacy DirectApex mutating builder.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;
use crate::types::CornerPatchInfo;

use super::cap_kind::CapKind;
use super::error::ApexBuildError;
use super::gregory_cap_boundary::gregory_sectors_from_plan_and_apex_edge_handles;
use super::plan::ApexPlan;

#[allow(clippy::too_many_arguments)]
pub(super) fn assign_cap_surface(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    apex_position: Point3,
    apex_vertex_handles: &[VertexHandle],
    apex_edge_handles: &[Option<EdgeHandle>],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Vec<CornerPatchInfo>, ApexBuildError> {
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    match &plan.cap_kind {
        CapKind::Sphere { center, radius, .. } => {
            brep.faces[cap_face].surface = Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Sphere {
                    center: *center,
                    radius: *radius,
                },
                orientation: true,
            });
            Ok(vec![CornerPatchInfo {
                vertex: plan.apex_vertex,
                n_sides: plan.cap_kind.n_fillets(),
                method: crate::types::CornerMethod::RollingBallSphere,
                achieved_continuity: crate::types::ContinuityLevel::G1,
                classification: crate::vertex_classifier::VertexBlendClassification::SameSideConvex,
                face: cap_face,
                bridge_fallbacks: Vec::new(),
                overlap_remediated: false,
                relaxed_continuity_coplanarity: None,
                cap_sphere: Some((*center, *radius)),
                sub_patch_provenance: None,
                setback_cap: false,
                exact_boundary_cap: false,
            }])
        }
        CapKind::Gregory { .. } => {
            let sectors = gregory_sectors_from_plan_and_apex_edge_handles(plan, apex_edge_handles)?;
            super::gregory_cap::assemble_gregory_cap_surface(
                brep,
                cap_face,
                plan,
                apex_position,
                apex_vertex_handles,
                apex_edge_handles,
                &sectors,
                fillet_surfaces,
                tolerance,
            )
            .map(|assembly| assembly.corner_patches)
        }
        CapKind::Open { .. } => unreachable!("Open handled in construct_apex_inner"),
    }
}
