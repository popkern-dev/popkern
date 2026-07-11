//! Imprint operation: project intersection curves of one body onto another,
//! splitting its faces at the intersection, without removing any material.
//!
//! The result is the target body with additional edges where the tool body
//! intersects it. Essential for downstream meshing, feature recognition,
//! and staged boolean workflows.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::edge_triple::{EdgeTripleOptions, compute_edge_geometry_unanchored};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::{clip_curve_to_faces, edge_triple_to_curve_on_surface, face_vertex_positions};
use crate::error::BoolError;
use crate::intersect::{IntersectionRecord, compute_face_aabbs};
use crate::split::split_along_curves;

/// Result of an imprint operation.
#[derive(Debug)]
pub struct ImprintResult {
    /// The target body with intersection curves split into its faces.
    /// Faces that crossed an intersection curve have been replaced by
    /// fragments along the curve; faces that did not cross any curve
    /// are unchanged. The tool body is consumed and not returned.
    pub imprinted_target: BRep,
    /// The intersection curves that drove the imprint, for inspection
    /// or to drive secondary operations.
    pub curves: Vec<ImprintCurve>,
    /// Number of face pairs tested.
    pub pairs_tested: usize,
}

/// A curve produced by the imprint operation.
#[derive(Debug, Clone)]
pub struct ImprintCurve {
    /// The intersection curve.
    pub curve: CurveOnSurface,
    /// Face on the target body that this curve lies on.
    pub target_face: FaceHandle,
    /// Face on the tool body that produced this curve.
    pub tool_face: FaceHandle,
}

/// Compute the imprint of a tool body onto a target body.
///
/// Finds all intersection curves between `target` and `tool`, then runs
/// the split phase of the boolean pipeline against the target body
/// only, replacing each crossed face by its fragments along the
/// intersection curves. The tool body is consumed and not returned;
/// the original `target` is unchanged (the result is built from a
/// clone).
///
/// Returns an [`ImprintResult`] containing the imprinted target, the
/// intersection curves, and a `pairs_tested` counter for telemetry.
pub fn imprint(
    target: &BRep,
    tool: &BRep,
    tolerance: &Tolerance,
) -> Result<ImprintResult, BoolError> {
    let aabbs_target = compute_face_aabbs(target);
    let aabbs_tool = compute_face_aabbs(tool);
    let mut curves = Vec::new();
    let mut records: Vec<IntersectionRecord> = Vec::new();
    let mut pairs_tested = 0;

    for (face_t, aabb_t) in &aabbs_target {
        for (face_tool, aabb_tool) in &aabbs_tool {
            if !aabb_t.overlaps(aabb_tool, tolerance.linear()) {
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

            // Stage 2b.59: migrated from `intersect_surfaces` to the
            // anchorless EdgeTriple API. Same pattern as section.rs
            // (stage 2b.58): face-pair discovery → unanchored
            // `Vec<EdgeTriple>` → unwrap each to `CurveOnSurface`
            // → clip → record.
            let opts = EdgeTripleOptions::default();
            let triples = match compute_edge_geometry_unanchored(
                surf_t,
                (u1_min, u1_max),
                (v1_min, v1_max),
                surf_tool,
                (u2_min, u2_max),
                (v2_min, v2_max),
                tolerance,
                &opts,
            ) {
                Ok(t) => t,
                Err(_) => continue,
            };
            for triple in triples {
                if let Some(curve) = edge_triple_to_curve_on_surface(&triple) {
                    for clipped in clip_curve_to_faces(target, *face_t, tool, *face_tool, &curve) {
                        curves.push(ImprintCurve {
                            curve: clipped.clone(),
                            target_face: *face_t,
                            tool_face: *face_tool,
                        });
                        records.push(IntersectionRecord {
                            face_a: *face_t,
                            face_b: *face_tool,
                            curve: clipped,
                            is_tangent: false,
                            is_coincident: false,
                            start_boundary: Default::default(),
                            end_boundary: Default::default(),
                        });
                    }
                }
            }
        }
    }

    // If no intersections, the target is unchanged.
    if records.is_empty() {
        return Ok(ImprintResult {
            imprinted_target: target.clone(),
            curves,
            pairs_tested,
        });
    }

    // Run the split phase against (target, tool) clones. The split
    // pipeline splits both inputs but only the target side is returned;
    // the tool clone is discarded.
    let split_result = split_along_curves(target.clone(), tool.clone(), &records, tolerance)
        .map_err(|diag| diag.into_error())?;

    Ok(ImprintResult {
        imprinted_target: split_result.a,
        curves,
        pairs_tested,
    })
}
