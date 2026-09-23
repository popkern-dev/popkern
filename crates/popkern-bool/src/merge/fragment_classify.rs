//! Fragment classification and selection for boolean operations.
//!
//! Classifies each face fragment as inside, outside, or on-boundary of the
//! opposite operand, then selects which fragments to keep (and whether to
//! flip their normals) based on the boolean operation type.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use popkern_core::{tolerance::Tolerance, topology::BRep};

use super::fragment_classify_overlap::classify_partial_overlap;
use super::fragment_sample::{sample_point_for_fragment, sample_point_uv};

#[path = "fragment_classify_deferred.rs"]
mod fragment_classify_deferred;
use super::graph::{IntersectionGraph, Operand};
use super::walk::Fragment;
use crate::classify::{point_in_solid_final, point_in_solid_final_with_confidence};
use crate::classify_on_face::classify_on_face;
use crate::coincidence::CoincidenceMap;
use crate::error::{BoolError, BoolOperation, Classification, PointClassification};
use crate::intersect::compute_face_aabbs;
use fragment_classify_deferred::apply_deferred_on_overrides;

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
    coincidence: &CoincidenceMap,
) -> Result<Vec<FragmentSelection>, BoolError> {
    // Precompute face AABBs for both operands (needed by point_in_solid).
    let aabbs_a = compute_face_aabbs(brep_a);
    let aabbs_b = compute_face_aabbs(brep_b);

    let mut selections = Vec::with_capacity(fragments.len());

    let mut witness_pending: Vec<(
        usize,
        Operand,
        PointClassification,
        Vec<crate::record_side_witness::RecordVote>,
    )> = Vec::new();
    for fragment in fragments {
        let own_operand = fragment.face.operand;
        let is_a = own_operand == Operand::A;

        // First check: is this fragment on a coincident face?
        // If so, use the ON-face rules from the spec (§5.1).
        if let Some(on_class) = classify_on_face(fragment.face.face, is_a, op, coincidence) {
            let selection = on_face_to_selection(on_class, own_operand, op);
            selections.push(selection);
            continue;
        }

        // Second check: is this fragment in the overlap region of a
        // PARTIAL coincident face? If the fragment's sample point is
        // inside both face polygons, apply coincidence rules.
        if let Some(on_class) =
            classify_partial_overlap(graph, fragment, is_a, op, coincidence, brep_a, brep_b)
        {
            let selection = on_face_to_selection(on_class, own_operand, op);
            selections.push(selection);
            continue;
        }

        // Normal path: sample point + ray casting.
        // Try the split-edge midpoint first, then fall back to the loop
        // centroid. Split-edge midpoints can land on the opposite operand's
        // boundary (where two faces share an edge), giving ambiguous results.
        let (_, other_brep, other_aabbs) = match own_operand {
            Operand::A => (Operand::A, brep_b, &aabbs_b),
            Operand::B => (Operand::B, brep_a, &aabbs_a),
        };

        let sample =
            sample_point_for_fragment(graph, &fragment.outer_loop, &fragment.face, brep_a, brep_b);
        // FINAL variant: fragments classify against the pristine INPUT
        // operand, so the input-body gates apply — plane hole hits are
        // genuine voids, and the cylinder angular/axial gates reject
        // phantom hits beyond a wall's banded span (measured, PK-341:
        // a ray from the base wall's mid-band fragment phantom-hit the
        // tube's outer wall ABOVE its z-band, parity flipped, the
        // fragment was kept, and the subtraction served an interior
        // membrane). The cellular pipeline already classifies with the
        // final variant; the graph-merge fallback was the outlier.
        let (mut classification, contested) = match point_in_solid_final_with_confidence(
            &sample,
            other_brep,
            other_aabbs,
            tolerance,
        ) {
            Ok(c) => c,
            Err(e) => {
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!(
                        "  [classify] {:?} sample ({:.3},{:.3},{:.3}) point_in_solid ERR {e:?} -> Outside",
                        fragment.face, sample.x, sample.y, sample.z
                    );
                }
                (PointClassification::Outside, false)
            }
        };
        // CONTESTED-VERDICT RESAMPLE: a split ray vote signals a
        // systematic per-ray miss (PK-460: a thread tool's groove
        // strip voted 3–2 Outside from strictly inside the tool).
        // The interior-step sample is an independent witness far
        // from the split edges; its verdict replaces the primary's
        // only when IT is uncontested.
        if contested {
            let resample = super::fragment_sample_interior::interior_step_sample(
                graph,
                &fragment.outer_loop,
                &fragment.face,
                brep_a,
                brep_b,
            );
            let reverdict = resample.as_ref().map(|p| {
                point_in_solid_final_with_confidence(p, other_brep, other_aabbs, tolerance)
            });
            if popkern_core::trace::trace_enabled("boolean") {
                eprintln!(
                    "  [classify] {:?} contested {classification:?} at ({:.3},{:.3},{:.3}) resample {:?} -> {reverdict:?}",
                    fragment.face,
                    sample.x,
                    sample.y,
                    sample.z,
                    resample.map(|p| (p.x, p.y, p.z)),
                );
            }
            if let Some(Ok((reclass, false))) = reverdict {
                classification = reclass;
            }
        }
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!(
                "  [classify] {:?} sample ({:.3},{:.3},{:.3}) -> {classification:?}",
                fragment.face, sample.x, sample.y, sample.z
            );
        }

        // When the primary sample results in rejection, retry with a
        // UV-space interior point for PLANAR faces only. On planar faces,
        // the UV↔3D mapping is linear, so the UV centroid reliably maps
        // to the face interior. Split-edge midpoints can land exactly on
        // the opposite operand's surface (at shared intersection edges),
        // giving ambiguous Outside results.
        // On curved faces (spheres, cylinders), UV centroids can map to
        // misleading 3D positions — skip the retry for those.
        // UV retry: only for B-fragments in subtraction on planar faces.
        // This targets the specific case where B's inner wall faces have
        // split-edge midpoints on A's surface (shared intersection edges).
        let eligible_for_uv_retry = !apply_selection_table(op, own_operand, classification).keep
            && own_operand == Operand::B
            && op == BoolOperation::Subtraction;
        if eligible_for_uv_retry {
            let own_brep = match own_operand {
                Operand::A => brep_a,
                Operand::B => brep_b,
            };
            let is_planar = own_brep
                .faces
                .get(fragment.face.face)
                .and_then(|f| f.surface.as_ref())
                .is_some_and(|s| {
                    matches!(
                        s.geometry,
                        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
                    )
                });
            if is_planar
                && let Some(uv_sample) =
                    sample_point_uv(graph, &fragment.outer_loop, &fragment.face, brep_a, brep_b)
            {
                let uv_class =
                    match point_in_solid_final(&uv_sample, other_brep, other_aabbs, tolerance) {
                        Ok(c) => c,
                        Err(_) => classification,
                    };
                if apply_selection_table(op, own_operand, uv_class).keep {
                    classification = uv_class;
                }
            }
        }

        // Record-side witness votes (PK-494 layer 3): collected per
        // OPERAND-A fragment, applied post-loop after per-face
        // outward-sign calibration — see `record_side_witness.rs`.
        // Collection matches the override scoping (A-only): a B-side
        // fragment's votes reference operand-A faces, whose sign
        // tallies no override ever consults — dead weight whose
        // per-vote obstruction gates priced O(records²) on multi-wrap
        // tool flanks (PK-496; each flank fragment carries a vote per
        // wrap record, each vote a dense surface×segment scan).
        if own_operand == Operand::A {
            let votes = crate::record_side_witness::record_votes(
                graph,
                fragment,
                &sample,
                brep_a,
                brep_b,
                other_aabbs,
                tolerance,
            );
            if !votes.is_empty() {
                witness_pending.push((selections.len(), own_operand, classification, votes));
            }
        }

        let selection = apply_selection_table(op, own_operand, classification);

        selections.push(selection);
    }

    // Post-pass: calibrate per-face outward signs from the ray-verdict
    // consensus and override minority mismatches (PK-494 layer 3).
    if !witness_pending.is_empty() {
        let all: Vec<(
            PointClassification,
            Vec<crate::record_side_witness::RecordVote>,
        )> = witness_pending
            .iter()
            .map(|(_, _, c, v)| (*c, v.clone()))
            .collect();
        let signs = crate::record_side_witness::calibrate_face_signs(&all);
        for (idx, operand, class, votes) in &witness_pending {
            // OPERAND-A-ONLY overrides: tool-side fragments carry walk
            // phantom semantics their ray verdicts encode deliberately
            // (PK-428 keyed-bore class — both measured regressions were
            // B-side, and confidence does not discriminate: the keyed
            // phantom was contested while the coil-end void-backed
            // sliver reads confidently wrong at 50+ turns). Workpiece
            // fragments carry no such semantics — the unrestricted
            // sweep regressed only B-side cases.
            if *operand != Operand::A {
                continue;
            }
            let Some(witness) = crate::record_side_witness::calibrated_verdict(votes, &signs)
            else {
                continue;
            };
            if witness == *class {
                continue;
            }
            if popkern_core::trace::trace_enabled("boolean")
                || popkern_core::trace::trace_enabled("uv-scramble")
            {
                eprintln!(
                    "  [record-side] fragment #{idx} ray {class:?} -> calibrated witness {witness:?}"
                );
            }
            selections[*idx] = apply_selection_table(op, *operand, witness);
        }
    }

    // §6.4: deferred ON resolution auto-override pass.
    //
    // Walk every fragment whose face appears in the coincidence map's
    // partial-overlap (ambiguous-on) list and ask
    // `classify_on_deferred::disambiguate_with_neighbors` whether the
    // immediate ON resolution should be replaced with a different
    // classification based on neighbour majority. If it would, swap
    // the selection.
    //
    // This is a no-op when the coincidence map has zero partial
    // overlaps (the existing test suite). When future inputs DO
    // produce partial overlaps, this pass picks up the more global
    // classification automatically.
    apply_deferred_on_overrides(fragments, &mut selections, brep_a, brep_b, op, coincidence);

    Ok(selections)
}

/// Convert a `Classification` (which may be ON-with-orientation) to
/// a `PointClassification` (which has only Inside/Outside/OnBoundary).
pub(super) fn classification_to_point(c: Classification) -> PointClassification {
    match c {
        Classification::Inside => PointClassification::Inside,
        Classification::Outside => PointClassification::Outside,
        Classification::OnBoundaryCooriented | Classification::OnBoundaryCrossing => {
            PointClassification::OnBoundary
        }
    }
}

fn on_face_to_selection(
    class: Classification,
    _operand: Operand,
    _op: BoolOperation,
) -> FragmentSelection {
    match class {
        Classification::Outside => FragmentSelection {
            keep: true,
            flip: false,
        },
        Classification::Inside => FragmentSelection {
            keep: false,
            flip: false,
        },
        Classification::OnBoundaryCooriented => FragmentSelection {
            keep: true,
            flip: false,
        },
        Classification::OnBoundaryCrossing => FragmentSelection {
            keep: true,
            flip: true,
        },
    }
}

/// Pick a representative sample point for a fragment.
///
/// Prefers the midpoint of a non-shared (original or split) edge, which is
/// guaranteed to lie well inside the fragment and away from the intersection
/// curve.  Falls back to the loop centroid if every edge is shared.
pub(super) fn apply_selection_table(
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
