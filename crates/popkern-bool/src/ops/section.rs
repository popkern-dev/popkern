//! Section operation: compute intersection curves between two bodies
//! without creating a new solid.
//!
//! The output is a set of curves (and possibly surface patches for
//! coincident regions) embedded in 3D. Used for interference analysis,
//! cross-section visualization, and as a building block for full booleans.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::edge_triple::{EdgeTripleOptions, compute_edge_geometry_unanchored};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::{clip_curve_to_faces, edge_triple_to_curve_on_surface, face_vertex_positions};
use crate::error::BoolError;
use crate::intersect::compute_face_aabbs;

/// Result of a section operation.
#[derive(Debug)]
pub struct SectionResult {
    /// Intersection curves between the two bodies.
    pub curves: Vec<SectionCurve>,
    /// Number of face pairs tested.
    pub pairs_tested: usize,
}

/// A single curve from a section operation.
#[derive(Debug, Clone)]
pub struct SectionCurve {
    /// The intersection curve with UV data on both surfaces.
    pub curve: CurveOnSurface,
    /// Face from body A that this curve lies on.
    pub face_a: FaceHandle,
    /// Face from body B that this curve lies on.
    pub face_b: FaceHandle,
}

/// Compute the section (intersection curves) between two bodies.
///
/// Returns all intersection curves without creating a new solid.
/// This is the intersection phase of a boolean without the
/// imprint/classify/assemble stages.
pub fn section(a: &BRep, b: &BRep, tolerance: &Tolerance) -> Result<SectionResult, BoolError> {
    let aabbs_a = compute_face_aabbs(a);
    let aabbs_b = compute_face_aabbs(b);
    let mut curves = Vec::new();
    let mut pairs_tested = 0;

    for (face_a, aabb_a) in &aabbs_a {
        for (face_b, aabb_b) in &aabbs_b {
            if !aabb_a.overlaps(aabb_b, tolerance.linear()) {
                continue;
            }
            pairs_tested += 1;

            let surf_a = match a.faces[*face_a].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };
            let surf_b = match b.faces[*face_b].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            let positions_a = face_vertex_positions(a, *face_a);
            let positions_b = face_vertex_positions(b, *face_b);
            let ((u1_min, u1_max), (v1_min, v1_max)) =
                popkern_core::geometry::uv_projection::bounded_parameter_range(
                    &positions_a,
                    surf_a,
                );
            let ((u2_min, u2_max), (v2_min, v2_max)) =
                popkern_core::geometry::uv_projection::bounded_parameter_range(
                    &positions_b,
                    surf_b,
                );

            // Stage 2b.58: migrated from `intersect_surfaces` to the
            // anchorless EdgeTriple API. Section discovery has no
            // natural anchors (the curve is being computed from
            // scratch), so we use `compute_edge_geometry_unanchored`
            // which mirrors `intersect_surfaces`'s UV-range form
            // and returns one `EdgeTriple` per intersection branch.
            let opts = EdgeTripleOptions::default();
            let triples = match compute_edge_geometry_unanchored(
                surf_a,
                (u1_min, u1_max),
                (v1_min, v1_max),
                surf_b,
                (u2_min, u2_max),
                (v2_min, v2_max),
                tolerance,
                &opts,
            ) {
                Ok(t) => t,
                Err(_) => continue,
            };
            for triple in triples {
                if let Some(cos) = edge_triple_to_curve_on_surface(&triple) {
                    for clipped in clip_curve_to_faces(a, *face_a, b, *face_b, &cos) {
                        curves.push(SectionCurve {
                            curve: clipped,
                            face_a: *face_a,
                            face_b: *face_b,
                        });
                    }
                }
            }
        }
    }

    Ok(SectionResult {
        curves,
        pairs_tested,
    })
}
