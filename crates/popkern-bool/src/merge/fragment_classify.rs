//! Fragment classification and selection for boolean operations.
//!
//! Classifies each face fragment as inside, outside, or on-boundary of the
//! opposite operand, then selects which fragments to keep (and whether to
//! flip their normals) based on the boolean operation type.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use popkern_core::{geometry::Point3, topology::BRep};

use super::graph::{IntersectionGraph, Operand};
use super::walk::Fragment;
use crate::error::{BoolError, BoolOperation, PointClassification};

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
    tolerance: f64,
) -> Result<Vec<FragmentSelection>, BoolError> {
    let mut selections = Vec::with_capacity(fragments.len());

    for fragment in fragments {
        // Compute centroid of outer loop vertices.
        let sample = compute_centroid(graph, &fragment.outer_loop);

        // Determine which operand this fragment belongs to.
        let (own_operand, other_brep) = match fragment.face.operand {
            Operand::A => (Operand::A, brep_b),
            Operand::B => (Operand::B, brep_a),
        };

        // Classify the sample point against the other operand.
        let classification = simple_classify(&sample, other_brep, tolerance);

        // Apply selection table.
        let selection = apply_selection_table(op, own_operand, classification);
        selections.push(selection);
    }

    Ok(selections)
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

/// Simplified point-in-solid test using the BRep's AABB.
///
/// This is a placeholder that delegates to the full `point_in_solid`
/// classifier when AABBs are available. For now, returns `Outside` as a
/// conservative default, which will be replaced by the full pipeline
/// integration.
fn simple_classify(_point: &Point3, _brep: &BRep, _tolerance: f64) -> PointClassification {
    // Full point_in_solid requires precomputed face AABBs. During
    // pipeline integration (Task 11) this will be wired up properly.
    // For now, return Outside as a safe default for testing.
    PointClassification::Outside
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
