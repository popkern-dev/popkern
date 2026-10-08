//! Far-fragment verdicts (PK-731): a fragment the intersection never
//! touched, lying wholly outside the other operand's bounding box, is
//! Outside the other operand — decided without a sample point or a ray.
//!
//! On a body cut many times by small tools (a 64-slot encoder disc)
//! almost every fragment of the target is such a face: before this,
//! each one paid the interior sampler and three ray votes per boolean,
//! the largest share of a late slot cut's merge.
//!
//! The verdict is exact. Every point of the fragment lies on its own
//! face, inside that face's containing box (`compute_face_aabbs` boxes
//! CONTAIN the face — they gate the face-pair SSI); a point outside the
//! other operand's box is outside that closed solid. The box margin
//! (`probe_eps_coarse`, the classifier's largest probe step, plus ten
//! linear tolerances) keeps every probe the ray path would have taken
//! outside it too. Only fragments with no record (`Shared`) half-edge
//! qualify: those are the only ones that carry record-side witness
//! votes, the one other consumer of a fragment's sample.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::topology::FaceHandle;

use crate::merge::graph::{EdgeSource, IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// The other operand's extent, padded, for each side.
pub(super) struct FarGate {
    /// Operand B's padded box (what A-side fragments test against).
    other_of_a: Option<Aabb>,
    /// Operand A's padded box.
    other_of_b: Option<Aabb>,
}

fn union_box(boxes: &[(FaceHandle, Aabb)], pad: f64) -> Option<Aabb> {
    let corners: Vec<Point3> = boxes.iter().flat_map(|(_, b)| [b.min, b.max]).collect();
    let bb = Aabb::from_points(&corners)?;
    Some(Aabb {
        min: Point3::new(bb.min.x - pad, bb.min.y - pad, bb.min.z - pad),
        max: Point3::new(bb.max.x + pad, bb.max.y + pad, bb.max.z + pad),
    })
}

impl FarGate {
    /// The gate for operands with face boxes `aabbs_a` / `aabbs_b`.
    pub(super) fn new(
        aabbs_a: &[(FaceHandle, Aabb)],
        aabbs_b: &[(FaceHandle, Aabb)],
        pad: f64,
    ) -> Self {
        FarGate {
            other_of_a: union_box(aabbs_b, pad),
            other_of_b: union_box(aabbs_a, pad),
        }
    }

    /// Whether `fragment` is Outside the other operand by extent alone
    /// (see the module docs). `own_boxes` are the fragment operand's
    /// face boxes.
    pub(super) fn is_far(
        &self,
        graph: &IntersectionGraph,
        fragment: &Fragment,
        own_boxes: &[(FaceHandle, Aabb)],
    ) -> bool {
        let other = match fragment.face.operand {
            Operand::A => self.other_of_a.as_ref(),
            Operand::B => self.other_of_b.as_ref(),
        };
        let Some(other) = other else {
            return false;
        };
        let cut = fragment
            .outer_loop
            .iter()
            .chain(fragment.inner_loops.iter().flatten())
            .any(|&h| matches!(graph.half_edges[h].source, EdgeSource::Shared { .. }));
        if cut {
            return false;
        }
        own_boxes
            .iter()
            .find(|(f, _)| *f == fragment.face.face)
            .is_some_and(|(_, own)| !own.overlaps(other, 0.0))
    }
}
