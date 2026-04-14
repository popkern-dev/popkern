//! Partial-overlap classification for fragments lying in the overlap
//! region of partial coincident faces.
//!
//! Split out of `fragment_classify.rs` to keep the parent under the
//! 500-line structural cap.

use popkern_core::topology::BRep;

use super::fragment_sample::compute_centroid;
use super::graph::IntersectionGraph;
use super::walk::Fragment;
use crate::coincidence::CoincidenceMap;
use crate::coincidence::types::CoincidenceOrientation;
use crate::error::{BoolOperation, Classification};

/// Check if a fragment lies in the overlap region of a partial coincident face.
///
/// For partial overlaps, `classify_on_face` returns `None` (it only handles
/// full overlaps). This function checks if the fragment's centroid is in the
/// overlap region (inside both face polygons in shared UV space). If so,
/// returns the coincidence classification. Otherwise returns `None` to fall
/// through to normal ray-casting.
pub(super) fn classify_partial_overlap(
    graph: &IntersectionGraph,
    fragment: &Fragment,
    operand_is_a: bool,
    operation: BoolOperation,
    coincidence: &CoincidenceMap,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<Classification> {
    use crate::coincidence::overlap_boundary::point_in_polygon_2d;
    use popkern_core::geometry::surface::SurfaceGeometry;

    // Find a partial overlap for this fragment's face.
    let overlap = coincidence.face_overlaps.iter().find(|fo| {
        let matches = if operand_is_a {
            fo.face_a == fragment.face.face
        } else {
            fo.face_b == fragment.face.face
        };
        matches && !fo.full_overlap
    })?;

    // Get both face surfaces — must be planar for UV projection.
    let surf_a = brep_a.faces.get(overlap.face_a)?.surface.as_ref()?;
    let SurfaceGeometry::Plane {
        origin,
        u_axis,
        v_axis,
        ..
    } = &surf_a.geometry
    else {
        return None;
    };

    // Project both face boundaries into shared UV space.
    let project = |brep: &BRep, fh: popkern_core::topology::FaceHandle| -> Vec<(f64, f64)> {
        crate::clip::face_vertex_positions(brep, fh)
            .iter()
            .map(|p| {
                let d = *p - *origin;
                (d.dot(u_axis), d.dot(v_axis))
            })
            .collect()
    };
    let poly_a = project(brep_a, overlap.face_a);
    let poly_b = project(brep_b, overlap.face_b);

    if poly_a.len() < 3 || poly_b.len() < 3 {
        return None;
    }

    // Get fragment centroid and project to UV.
    let centroid = compute_centroid(graph, &fragment.outer_loop);
    let d = centroid - *origin;
    let u = d.dot(u_axis);
    let v = d.dot(v_axis);

    // Check if centroid is inside BOTH polygons (= in the overlap region).
    if point_in_polygon_2d(u, v, &poly_a) && point_in_polygon_2d(u, v, &poly_b) {
        // In the overlap region — apply coincidence rules.
        let orientation = overlap.orientation;
        Some(match (operation, orientation) {
            (BoolOperation::Subtraction, CoincidenceOrientation::Same) => Classification::Inside,
            (BoolOperation::Subtraction, CoincidenceOrientation::Opposite) => {
                if operand_is_a {
                    Classification::OnBoundaryCrossing
                } else {
                    Classification::Inside
                }
            }
            (BoolOperation::Union, CoincidenceOrientation::Same) => {
                if operand_is_a {
                    Classification::OnBoundaryCooriented
                } else {
                    Classification::Inside
                }
            }
            (BoolOperation::Union, CoincidenceOrientation::Opposite) => Classification::Inside,
            (BoolOperation::Intersection, CoincidenceOrientation::Same) => {
                if operand_is_a {
                    Classification::OnBoundaryCooriented
                } else {
                    Classification::Outside
                }
            }
            (BoolOperation::Intersection, CoincidenceOrientation::Opposite) => {
                Classification::Outside
            }
        })
    } else {
        None // Not in overlap region — fall through to ray casting.
    }
}
