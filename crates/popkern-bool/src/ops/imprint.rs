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
use crate::split::{SplitSide, split_body_along_curves};

pub use super::imprint_classify::FaceClass;

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
    /// Where each face of `imprinted_target` lies relative to the tool.
    pub classes: Vec<(FaceHandle, FaceClass)>,
}

impl ImprintResult {
    /// The class of `face` of the imprinted target, when decided.
    pub fn class_of(&self, face: FaceHandle) -> Option<FaceClass> {
        self.classes
            .iter()
            .find(|(f, _)| *f == face)
            .map(|(_, c)| *c)
    }

    /// The imprinted target's faces of `class`.
    pub fn faces_classified(&self, class: FaceClass) -> impl Iterator<Item = FaceHandle> + '_ {
        self.classes
            .iter()
            .filter(move |(_, c)| *c == class)
            .map(|(f, _)| *f)
    }
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
                    for clipped in
                        clip_curve_to_faces(target, *face_t, tool, *face_tool, &curve, None)
                    {
                        curves.push(ImprintCurve {
                            curve: clipped.clone(),
                            target_face: *face_t,
                            tool_face: *face_tool,
                        });
                        records.push(IntersectionRecord {
                            residue_owner_a: None,
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
    let imprinted_target = if records.is_empty() {
        target.clone()
    } else {
        imprinted(target, tool, &records, tolerance)?
    };
    let classes = super::imprint_classify::classify_faces(&imprinted_target, tool, tolerance);
    Ok(ImprintResult {
        imprinted_target,
        curves,
        pairs_tested,
        classes,
    })
}

/// The target's faces split by the tool: as the cell complex splits
/// them when both are solids (the same identity booleans use), else (a
/// sheet, or a partition that declines, recorded) along the records'
/// target-side curves.
fn imprinted(
    target: &BRep,
    tool: &BRep,
    records: &[IntersectionRecord],
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    use popkern_core::topology::BodyKind;
    if target.body_kind() == BodyKind::Solid && tool.body_kind() == BodyKind::Solid {
        match crate::partition::partition(&[target, tool], *tolerance).and_then(|p| p.imprinted(0))
        {
            Ok(body) => return Ok(body),
            Err(e) => crate::complex_shadow::append_log_line(&format!("ImprintFallback {e}")),
        }
    }
    let (imprinted_target, _) =
        split_body_along_curves(target.clone(), SplitSide::A, records, tolerance)
            .map_err(|diag| diag.into_error())?;
    Ok(imprinted_target)
}
