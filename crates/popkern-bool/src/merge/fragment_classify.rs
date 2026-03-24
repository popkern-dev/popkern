//! Fragment classification and selection for boolean operations.
//!
//! Classifies each face fragment as inside, outside, or on-boundary of the
//! opposite operand, then selects which fragments to keep (and whether to
//! flip their normals) based on the boolean operation type.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use popkern_core::{geometry::Point3, tolerance::Tolerance, topology::BRep};

use super::graph::{EdgeSource, IntersectionGraph, Operand};
use super::walk::Fragment;
use crate::classify::point_in_solid;
use crate::error::{BoolError, BoolOperation, PointClassification};
use crate::intersect::compute_face_aabbs;

/// Selection decision for a fragment.
#[derive(Debug)]
pub(crate) struct FragmentSelection {
    /// Whether to include this fragment in the result.
    pub keep: bool,
    /// Whether to flip the surface normal of this fragment.
    pub flip: bool,
}

/// Classify and select fragments based on the boolean operation.
///
/// For each fragment, computes a sample point (centroid of outer loop vertices),
/// tests it against the opposite operand via point-in-solid classification,
/// and applies the selection table for the given operation.
///
/// # Errors
///
/// Returns [`BoolError`] if point-in-solid classification fails for any
/// fragment.
pub(crate) fn classify_and_select_fragments(
    fragments: &[Fragment],
    graph: &IntersectionGraph,
    brep_a: &BRep,
    brep_b: &BRep,
    op: BoolOperation,
    tolerance: &Tolerance,
) -> Result<Vec<FragmentSelection>, BoolError> {
    // Precompute face AABBs for both operands (needed by point_in_solid).
    let aabbs_a = compute_face_aabbs(brep_a);
    let aabbs_b = compute_face_aabbs(brep_b);

    let mut selections = Vec::with_capacity(fragments.len());

    for fragment in fragments {
        // Use the midpoint of a non-shared (original/split) edge as the
        // sample point.  The centroid of the entire loop can land inside the
        // opposite solid for *both* fragments of a split face, causing wrong
        // classification.  An original edge's midpoint is guaranteed to be
        // well inside the fragment and away from the intersection curve.
        let sample = sample_point_for_fragment(graph, &fragment.outer_loop);

        // Determine which operand this fragment belongs to.
        let (own_operand, other_brep, other_aabbs) = match fragment.face.operand {
            Operand::A => (Operand::A, brep_b, &aabbs_b),
            Operand::B => (Operand::B, brep_a, &aabbs_a),
        };

        // Classify the sample point against the other operand.
        let classification = match point_in_solid(&sample, other_brep, other_aabbs, tolerance) {
            Ok(c) => c,
            Err(_) => {
                // Fall back to Outside if classification fails (e.g. no
                // conclusive ray). This is conservative — fragments default
                // to being kept for union, discarded for intersection.
                PointClassification::Outside
            }
        };

        // Apply selection table.
        let selection = apply_selection_table(op, own_operand, classification);
        selections.push(selection);
    }

    Ok(selections)
}

/// Pick a representative sample point for a fragment.
///
/// Prefers the midpoint of a non-shared (original or split) edge, which is
/// guaranteed to lie well inside the fragment and away from the intersection
/// curve.  Falls back to the loop centroid if every edge is shared.
fn sample_point_for_fragment(graph: &IntersectionGraph, loop_hes: &[usize]) -> Point3 {
    for &he_idx in loop_hes {
        let he = &graph.half_edges[he_idx];
        match &he.source {
            EdgeSource::OriginalA(_)
            | EdgeSource::OriginalB(_)
            | EdgeSource::SplitA { .. }
            | EdgeSource::SplitB { .. } => {
                let start = graph.nodes[he.origin].position;
                let twin_idx = he.twin;
                let end = graph.nodes[graph.half_edges[twin_idx].origin].position;
                return Point3::new(
                    (start.x + end.x) * 0.5,
                    (start.y + end.y) * 0.5,
                    (start.z + end.z) * 0.5,
                );
            }
            EdgeSource::Shared { .. } => continue,
        }
    }
    // Fallback: all edges are shared (shouldn't normally happen).
    compute_centroid(graph, loop_hes)
}

/// Compute the centroid of a loop's vertex positions.
fn compute_centroid(graph: &IntersectionGraph, loop_hes: &[usize]) -> Point3 {
    if loop_hes.is_empty() {
        return Point3::new(0.0, 0.0, 0.0);
    }

    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;

    for &he_idx in loop_hes {
        let origin = graph.half_edges[he_idx].origin;
        let pos = graph.nodes[origin].position;
        sum_x += pos.x;
        sum_y += pos.y;
        sum_z += pos.z;
    }

    let n = loop_hes.len() as f64;
    Point3::new(sum_x / n, sum_y / n, sum_z / n)
}

/// Apply the boolean selection table.
///
/// Returns a [`FragmentSelection`] indicating whether to keep the fragment
/// and whether to flip its normal.
fn apply_selection_table(
    op: BoolOperation,
    operand: Operand,
    class: PointClassification,
) -> FragmentSelection {
    match op {
        BoolOperation::Union => {
            // Keep outside fragments from both operands.
            // Keep on-boundary from A only.
            match (operand, class) {
                (_, PointClassification::Outside) => FragmentSelection {
                    keep: true,
                    flip: false,
                },
                (Operand::A, PointClassification::OnBoundary) => FragmentSelection {
                    keep: true,
                    flip: false,
                },
                _ => FragmentSelection {
                    keep: false,
                    flip: false,
                },
            }
        }
        BoolOperation::Intersection => {
            // Keep inside fragments from both operands.
            // Keep on-boundary from A only.
            match (operand, class) {
                (_, PointClassification::Inside) => FragmentSelection {
                    keep: true,
                    flip: false,
                },
                (Operand::A, PointClassification::OnBoundary) => FragmentSelection {
                    keep: true,
                    flip: false,
                },
                _ => FragmentSelection {
                    keep: false,
                    flip: false,
                },
            }
        }
        BoolOperation::Subtraction => {
            // Keep outside from A, inside from B (flipped).
            // On-boundary from A kept.
            match (operand, class) {
                (Operand::A, PointClassification::Outside) => FragmentSelection {
                    keep: true,
                    flip: false,
                },
                (Operand::A, PointClassification::OnBoundary) => FragmentSelection {
                    keep: true,
                    flip: false,
                },
                (Operand::B, PointClassification::Inside) => FragmentSelection {
                    keep: true,
                    flip: true,
                },
                _ => FragmentSelection {
                    keep: false,
                    flip: false,
                },
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::PointClassification;
    use crate::merge::graph::Operand;

    #[test]
    fn union_keeps_outside_both_operands() {
        let sel = apply_selection_table(
            BoolOperation::Union,
            Operand::A,
            PointClassification::Outside,
        );
        assert!(sel.keep);
        assert!(!sel.flip);

        let sel = apply_selection_table(
            BoolOperation::Union,
            Operand::B,
            PointClassification::Outside,
        );
        assert!(sel.keep);
        assert!(!sel.flip);
    }

    #[test]
    fn union_discards_inside() {
        let sel = apply_selection_table(
            BoolOperation::Union,
            Operand::A,
            PointClassification::Inside,
        );
        assert!(!sel.keep);
    }

    #[test]
    fn intersection_keeps_inside_both() {
        let sel = apply_selection_table(
            BoolOperation::Intersection,
            Operand::A,
            PointClassification::Inside,
        );
        assert!(sel.keep);
        assert!(!sel.flip);

        let sel = apply_selection_table(
            BoolOperation::Intersection,
            Operand::B,
            PointClassification::Inside,
        );
        assert!(sel.keep);
    }

    #[test]
    fn subtraction_keeps_outside_a_and_flips_inside_b() {
        let sel = apply_selection_table(
            BoolOperation::Subtraction,
            Operand::A,
            PointClassification::Outside,
        );
        assert!(sel.keep);
        assert!(!sel.flip);

        let sel = apply_selection_table(
            BoolOperation::Subtraction,
            Operand::B,
            PointClassification::Inside,
        );
        assert!(sel.keep);
        assert!(sel.flip);
    }

    #[test]
    fn subtraction_discards_outside_b() {
        let sel = apply_selection_table(
            BoolOperation::Subtraction,
            Operand::B,
            PointClassification::Outside,
        );
        assert!(!sel.keep);
    }

    #[test]
    fn on_boundary_kept_for_a_in_all_ops() {
        for op in &[
            BoolOperation::Union,
            BoolOperation::Intersection,
            BoolOperation::Subtraction,
        ] {
            let sel = apply_selection_table(*op, Operand::A, PointClassification::OnBoundary);
            assert!(sel.keep, "OnBoundary from A should be kept for {op:?}");
        }
    }
}
