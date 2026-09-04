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
    // Find a partial overlap for this fragment's face.
    let overlap = coincidence.face_overlaps.iter().find(|fo| {
        let matches = if operand_is_a {
            fo.face_a == fragment.face.face
        } else {
            fo.face_b == fragment.face.face
        };
        matches && !fo.full_overlap
    })?;

    // The overlap region is where BOTH coincident faces cover the same
    // plane. Test the fragment's 3D centroid for containment in EACH
    // face's own 2D region via the curve-aware `point_inside_face_uv`
    // (PK-435: the earlier polygon-projection test bailed on a disc —
    // fewer than three boundary vertices — so a coincident disc/rect
    // overlap never applied the keep-one-skin rule and the disc's lens
    // fragment fell to ray-casting on the coincident plane, where ray
    // parity is ill-defined and dropped the kept union skin).
    let centroid = compute_centroid(graph, &fragment.outer_loop);
    if !point_in_face_plane(brep_a, overlap.face_a, centroid)
        || !point_in_face_plane(brep_b, overlap.face_b, centroid)
    {
        return None; // Not in the overlap region — fall through to ray casting.
    }

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
        (BoolOperation::Intersection, CoincidenceOrientation::Opposite) => Classification::Outside,
    })
}

/// Curve-aware test: does `point` (a 3D position on or near the shared
/// plane) lie inside `face`'s 2D region? Projects `point` into the
/// face's planar UV and defers to `point_inside_face_uv`, which walks
/// the face's actual boundary curves (so a disc's circular rim is
/// handled, unlike a vertex-polygon projection).
fn point_in_face_plane(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
    point: popkern_core::geometry::Point3,
) -> bool {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::surface::SurfaceGeometry;
    let Some(sof) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let SurfaceGeometry::Plane {
        origin,
        u_axis,
        v_axis,
        ..
    } = &sof.geometry
    else {
        return false;
    };
    let d = point - *origin;
    let uv = Point3::new(d.dot(u_axis), d.dot(v_axis), 0.0);
    crate::split_inner_loop::point_inside_face_uv(brep, face, uv, &brep.tolerance)
}
