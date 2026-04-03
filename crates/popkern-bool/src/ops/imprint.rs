//! Imprint operation: project intersection curves of one body onto another,
//! splitting its faces at the intersection, without removing any material.
//!
//! The result is the target body with additional edges where the tool body
//! intersects it. Essential for downstream meshing, feature recognition,
//! and staged boolean workflows.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::intersection::intersect_surfaces;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::{clip_curve_to_faces, face_vertex_positions};
use crate::error::BoolError;
use crate::intersect::compute_face_aabbs;

/// Result of an imprint operation.
#[derive(Debug)]
pub struct ImprintResult {
    /// The imprinted body (target with intersection curves added as edges).
    /// Note: actual face splitting requires the split pipeline, which is
    /// not yet integrated here. This result contains the computed curves
    /// that would be used for splitting.
    pub curves: Vec<ImprintCurve>,
    /// Number of face pairs tested.
    pub pairs_tested: usize,
}

/// A curve to be imprinted onto the target body.
#[derive(Debug, Clone)]
pub struct ImprintCurve {
    /// The intersection curve.
    pub curve: CurveOnSurface,
    /// Face on the target body to split.
    pub target_face: FaceHandle,
    /// Face on the tool body that produced this curve.
    pub tool_face: FaceHandle,
}

/// Compute the imprint of a tool body onto a target body.
///
/// Finds all intersection curves between `target` and `tool`, returning
/// curves that lie on the target body's faces. These curves can be used
/// to split the target's faces without removing material.
pub fn imprint(
    target: &BRep,
    tool: &BRep,
    tolerance: &Tolerance,
) -> Result<ImprintResult, BoolError> {
    let aabbs_target = compute_face_aabbs(target);
    let aabbs_tool = compute_face_aabbs(tool);
    let mut curves = Vec::new();
    let mut pairs_tested = 0;

    for (face_t, aabb_t) in &aabbs_target {
        for (face_tool, aabb_tool) in &aabbs_tool {
            if !aabb_t.overlaps(aabb_tool, tolerance.linear) {
                continue;
            }
            pairs_tested += 1;

            let surf_t = match target.faces[*face_t].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };
            let surf_tool = match tool.faces[*face_tool].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            let positions_t = face_vertex_positions(target, *face_t);
            let positions_tool = face_vertex_positions(tool, *face_tool);
            let ((u1_min, u1_max), (v1_min, v1_max)) =
                popkern_core::geometry::uv_projection::bounded_parameter_range(
                    &positions_t,
                    surf_t,
                );
            let ((u2_min, u2_max), (v2_min, v2_max)) =
                popkern_core::geometry::uv_projection::bounded_parameter_range(
                    &positions_tool,
                    surf_tool,
                );

            if let Ok(result) = intersect_surfaces(
                surf_t,
                (u1_min, u1_max),
                (v1_min, v1_max),
                surf_tool,
                (u2_min, u2_max),
                (v2_min, v2_max),
                tolerance,
            ) {
                for isect in result.intersections {
                    if let popkern_core::geometry::intersection::SurfaceSurfaceIntersection::Curve(
                        curve,
                    ) = isect
                        && let Some(clipped) =
                            clip_curve_to_faces(target, *face_t, tool, *face_tool, &curve)
                    {
                        curves.push(ImprintCurve {
                            curve: clipped,
                            target_face: *face_t,
                            tool_face: *face_tool,
                        });
                    }
                }
            }
        }
    }

    Ok(ImprintResult {
        curves,
        pairs_tested,
    })
}
