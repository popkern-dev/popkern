use crate::ExactBoundaryCandidateOutcome;
use crate::types::{BlendResult, CornerMethod, CornerPatchInfo};
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{FaceHandle, VertexHandle};

pub(super) fn gregory_corners_for_vertex(
    result: &BlendResult,
    vertex: VertexHandle,
) -> Vec<&CornerPatchInfo> {
    result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex && p.method == CornerMethod::Gregory)
        .collect()
}

pub(super) fn rectangular_gregory_corner_faces(
    result: &BlendResult,
    vertex: VertexHandle,
) -> Vec<FaceHandle> {
    gregory_corners_for_vertex(result, vertex)
        .into_iter()
        .filter_map(|p| {
            let geometry = &result.brep.faces.get(p.face)?.surface.as_ref()?.geometry;
            matches!(
                geometry,
                SurfaceGeometry::Gregory(
                    GregoryPatch::Rectangular(_) | GregoryPatch::RectangularDeg5(_)
                )
            )
            .then_some(p.face)
        })
        .collect()
}

pub(super) fn assert_logical_gregory_corner(
    result: &BlendResult,
    vertex: VertexHandle,
    n_sides: usize,
    context: &str,
) {
    let patches = gregory_corners_for_vertex(result, vertex);
    if patches.is_empty()
        && result.repair_report.apex_routes.iter().any(|route| {
            route.exact_boundary_report.candidate_status == ExactBoundaryCandidateOutcome::Accepted
                && !route.exact_boundary_report.committed_faces.is_empty()
        })
    {
        return;
    }
    assert!(
        patches.len() == 1 || patches.len() == n_sides,
        "{context}: expected either one logical Gregory corner patch or one per-sector Gregory \
         record when the rectangular subpatch path is not committed, got {}",
        patches.len()
    );
    for patch in patches {
        assert_eq!(patch.n_sides, n_sides, "{context}: n_sides drifted");
        assert!(
            result.brep.faces[patch.face].surface.is_some(),
            "{context}: logical Gregory corner face must carry a surface"
        );
    }
}
