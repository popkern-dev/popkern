//! Fragment classification for the cellular pipeline (PK-496 split
//! from `boolean.rs` under the 500-line budget): per-fragment
//! sampling via `fragment_classify_one`, the PK-494 record-side
//! witness post-pass, and the A-side verdict-propagation planner
//! (`fragment_propagation`) that replaces O(fragments) rays with
//! O(components + audits) on multi-wrap workpiece faces.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::{BoolDiagnostics, BoolOperation, PointClassification};
use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::{IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// Keep-one-skin selection rules with contact ORIENTATION where
/// known: COORIENTED contact (both outward normals agree) means the
/// two boundaries coincide as the same skin — union/intersect keep
/// exactly one copy (operand A's), and subtract drops A's piece (the
/// removed region sits directly behind it, opening the mouth).
/// ANTI-oriented contact means the bodies face each other across the
/// surface — the union of the two materials makes it interior
/// (drop), while a subtract leaves the wall intact (keep). Unknown
/// orientation keeps the historical cooriented assumption.
fn selection_for(
    op: BoolOperation,
    own_operand: Operand,
    class: PointClassification,
    cooriented: Option<bool>,
) -> FragmentSelection {
    let (keep, flip) = match op {
        BoolOperation::Union => (
            class == PointClassification::Outside
                || (own_operand == Operand::A
                    && class == PointClassification::OnBoundary
                    && cooriented != Some(false)),
            false,
        ),
        BoolOperation::Intersection => (
            class == PointClassification::Inside
                || (own_operand == Operand::A
                    && class == PointClassification::OnBoundary
                    && cooriented != Some(false)),
            false,
        ),
        BoolOperation::Subtraction => match own_operand {
            Operand::A => (
                class == PointClassification::Outside
                    || (class == PointClassification::OnBoundary && cooriented != Some(true)),
                false,
            ),
            Operand::B => (class == PointClassification::Inside, true),
        },
    };
    FragmentSelection { keep, flip }
}

/// DEGENERATE loops (a lone non-closed half-edge, or a he/twin
/// bigon) enclose no area — boundary-graze residue. Keeping one
/// produced a dangling face in the assembly.
fn is_degenerate_fragment(
    graph: &IntersectionGraph,
    fragment: &Fragment,
    tolerance: &Tolerance,
) -> bool {
    match fragment.outer_loop.len() {
        1 => {
            let he = &graph.half_edges[fragment.outer_loop[0]];
            if he.origin != graph.half_edges[he.twin].origin {
                true
            } else {
                // A lone SELF-LOOP is a real fragment boundary when
                // it is a closed intersection curve (a circle cut
                // through a face). A ZERO-EXTENT self-loop is a
                // tangency-point record (a surface grazing another
                // at one point) — keeping it assembled a point-face
                // that broke the solid's Euler count. Cull by
                // sampled curve extent against the tolerance.
                let extent = he
                    .curve
                    .as_ref()
                    .or_else(|| graph.half_edges[he.twin].curve.as_ref())
                    .map(|c| {
                        use popkern_core::geometry::curve::CurveEval;
                        let origin = graph.nodes[he.origin].position;
                        (1..8)
                            .map(|i| {
                                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                                c.geometry.point_at(t).distance_to(&origin)
                            })
                            .fold(0.0, f64::max)
                    })
                    .unwrap_or(0.0);
                extent < tolerance.linear() * 10.0
            }
        }
        2 => {
            let a = fragment.outer_loop[0];
            if graph.half_edges[a].twin == fragment.outer_loop[1] {
                // A he/twin bigon is graze residue only for RECORD
                // pairs. A BOUNDARY pair walked as one loop is a
                // periodic seam spine — the outer loop of a full
                // periodic chart (a full sphere's two meridian
                // copies bound the whole wall face, with the
                // imprint loops as its holes).
                matches!(
                    graph.half_edges[a].source,
                    crate::merge::graph::EdgeSource::Shared { .. }
                )
            } else {
                // An OPEN two-half-edge walk with (sampled) zero
                // area is residue of a seeded walk that swallowed
                // a tangency self-loop and dead-ended — e.g. the
                // unconsumed direction of a boundary-collinear
                // corner record plus the zero-extent point-contact
                // record at its far end. It can never assemble
                // (measured: kept as OnBoundary, it became an
                // unclosed 20 mm loop and the whole flush-arm
                // chain fell to graph-merge). Real two-arc bigons
                // (a split circle) have genuine sampled area and
                // stay.
                let start = graph.half_edges[fragment.outer_loop[0]].origin;
                let last = fragment.outer_loop[1];
                let end = graph.half_edges[graph.half_edges[last].twin].origin;
                let open = start != end;
                let area =
                    crate::merge::walk::walk_geometry::loop_newell(graph, &fragment.outer_loop)
                        .length();
                open && area < (tolerance.linear() * 10.0).powi(2)
            }
        }
        _ => false,
    }
}

/// The opposite plain verdict (propagation flips across transversal
/// records; only Inside/Outside ever propagate).
fn flipped(class: PointClassification) -> PointClassification {
    match class {
        PointClassification::Inside => PointClassification::Outside,
        PointClassification::Outside => PointClassification::Inside,
        other => other,
    }
}

/// Classify fragments using cellular logic: determine which volumetric
/// cell each fragment borders, then select based on cell membership.
///
/// For each fragment:
/// 1. Sample a point on the fragment
/// 2. Test it against the OTHER body via point-in-solid — or, on an
///    audited A-side propagation component, derive the verdict from
///    the component seed's ray with flip parity (PK-496)
/// 3. Apply the operation predicate
#[allow(clippy::too_many_arguments)]
pub(super) fn classify_fragments_cellular(
    fragments: &[Fragment],
    graph: &IntersectionGraph,
    brep_a: &BRep,
    brep_b: &BRep,
    op: BoolOperation,
    tolerance: &Tolerance,
    face_overlaps: &[crate::coincidence::types::FaceOverlap],
) -> Result<Vec<FragmentSelection>, BoolDiagnostics> {
    use crate::intersect::compute_face_aabbs;

    let aabbs_a = compute_face_aabbs(brep_a);
    let aabbs_b = compute_face_aabbs(brep_b);

    // Faces participating in a COINCIDENT overlap pair: fragments of
    // these lie in the other body's boundary plane, so a raw sample
    // has ill-defined ray parity (measured: two mirror-image cap
    // fragments classified Inside on one side and Outside on the
    // other — the kept duplicate skin went non-manifold). They are
    // routed to the offset-pair classifier, whose (in, out) verdict
    // maps to OnBoundary and the operation's keep-one-skin rules.
    let overlap_a: std::collections::HashSet<_> =
        face_overlaps.iter().map(|fo| fo.face_a).collect();
    let overlap_b: std::collections::HashSet<_> =
        face_overlaps.iter().map(|fo| fo.face_b).collect();

    // Probe offset for side/coorientation sampling. It must sit
    // decades ABOVE the linear tolerance band (so boundary parity
    // noise cannot flip a verdict) yet decades BELOW feature scale
    // (so the probe point stays within the local cell). The old
    // `linear × 1e3` scaling honored only the first bound: at a
    // coarse linear tolerance of 1e-4 the ±16× retry probes stepped
    // 1.6 — clear outside a unit-sized model — and every
    // coincident-skin verdict degraded to parity noise (measured:
    // flush cube∪cube at 1e-4 kept both overlap skins, χ = −10,
    // 16 open edges). Cap by the combined model diagonal; the floor
    // stays above the tolerance band even for tiny models.
    let model_diag = {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for brep in [brep_a, brep_b] {
            for (_, v) in &brep.vertices {
                for (i, c) in [v.position.x, v.position.y, v.position.z]
                    .into_iter()
                    .enumerate()
                {
                    lo[i] = lo[i].min(c);
                    hi[i] = hi[i].max(c);
                }
            }
        }
        let d = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    };
    let probe_eps = if model_diag.is_finite() && model_diag > 0.0 {
        (tolerance.linear() * 1.0e3)
            .max(1.0e-3)
            .min(model_diag * 2.0e-3)
            .max(tolerance.linear() * 8.0)
    } else {
        (tolerance.linear() * 1.0e3).max(1.0e-3)
    };
    // The coarser retry scale (corner-graze disambiguation) obeys the
    // same upper bound: an unclamped ×16 retry at linear=1e-2 probed
    // 1.28 — half the model away — read the far side of the solid,
    // and the trust-the-coarser-scale rule overrode a CORRECT
    // fine-scale verdict (measured: cube∪sphere at 1e-2 dropped the
    // sphere wall, χ = 8). When the clamp collapses the ladder the
    // retry degenerates to the fine scale, which is simply a no-op.
    let probe_eps_coarse = if model_diag.is_finite() && model_diag > 0.0 {
        (probe_eps * 16.0).min(model_diag * 5.0e-2).max(probe_eps)
    } else {
        probe_eps * 16.0
    };

    let ctx = super::fragment_classify_one::ClassifyCtx {
        graph,
        brep_a,
        brep_b,
        aabbs_a: &aabbs_a,
        aabbs_b: &aabbs_b,
        overlap_a: &overlap_a,
        overlap_b: &overlap_b,
        probe_eps,
        probe_eps_coarse,
        tolerance,
    };

    let n = fragments.len();
    let degenerate: Vec<bool> = fragments
        .iter()
        .map(|f| is_degenerate_fragment(graph, f, tolerance))
        .collect();

    // A-side verdict propagation plan (PK-496): audited components
    // classify from one seed ray + flip parity instead of a ray per
    // wrap-band fragment.
    let plan = if super::fragment_propagation::propagation_enabled() {
        super::fragment_propagation::build_plan(
            graph,
            fragments,
            &degenerate,
            &overlap_a,
            brep_a,
            brep_b,
        )
    } else {
        None
    };
    let deferred: Vec<Option<usize>> = match &plan {
        Some(p) => p.deferred.clone(),
        None => vec![None; n],
    };

    // Pass 1: sample-and-classify everything that is not deferred to
    // a component verdict (seeds and audits are NOT deferred).
    let mut results: Vec<Option<(Point3, PointClassification, Option<bool>)>> =
        (0..n).map(|_| None).collect();
    for i in 0..n {
        if degenerate[i] || deferred[i].is_some() {
            continue;
        }
        results[i] = Some(super::fragment_classify_one::sample_and_classify(
            &ctx,
            &fragments[i],
        ));
    }

    // Pass 2: resolve components. A seed that is not plain In/Out, or
    // any audit ray disagreeing with its propagated expectation,
    // declines the component back to the per-fragment path.
    let mut propagated = vec![false; n];
    if let Some(plan) = &plan {
        for comp in &plan.components {
            let seed_class = results[comp.seed].as_ref().map(|r| r.1);
            let propagatable = matches!(
                seed_class,
                Some(PointClassification::Inside | PointClassification::Outside)
            );
            let audits_ok = propagatable
                && comp.audits.iter().all(|&(a, par)| {
                    let expected = seed_class.map(|c| if par { flipped(c) } else { c });
                    results[a].as_ref().map(|r| r.1) == expected
                });
            if audits_ok {
                let mut filled = 0usize;
                for &(m, par) in &comp.members {
                    if deferred[m].is_some()
                        && results[m].is_none()
                        && let Some(seed) = seed_class
                    {
                        let class = if par { flipped(seed) } else { seed };
                        let sample = super::fragment_sample_helpers::sample_fragment_point(
                            graph,
                            &fragments[m].outer_loop,
                        );
                        results[m] = Some((sample, class, None));
                        propagated[m] = true;
                        filled += 1;
                    }
                }
                if popkern_core::trace::trace_enabled("cellular") {
                    eprintln!(
                        "  [propagate] component of {}: seed {:?} {seed_class:?}, {} audits agree, {filled} propagated",
                        comp.members.len(),
                        fragments[comp.seed].face,
                        comp.audits.len()
                    );
                }
            } else {
                if popkern_core::trace::trace_enabled("cellular") {
                    eprintln!(
                        "  [propagate] component of {} DECLINED (seed {seed_class:?}, audit mismatch or non-plain) — per-fragment",
                        comp.members.len()
                    );
                }
                for &(m, _) in &comp.members {
                    if results[m].is_none() {
                        results[m] = Some(super::fragment_classify_one::sample_and_classify(
                            &ctx,
                            &fragments[m],
                        ));
                    }
                }
            }
        }
    }

    // Pass 3: selections in fragment order + witness vote collection.
    let mut selections = Vec::with_capacity(n);
    let mut witness_pending: Vec<(
        usize,
        PointClassification,
        Option<bool>,
        Vec<crate::record_side_witness::RecordVote>,
    )> = Vec::new();
    for (i, fragment) in fragments.iter().enumerate() {
        if degenerate[i] {
            if popkern_core::trace::trace_enabled("cellular") {
                eprintln!(
                    "  [sel] {:?} outer={} DEGENERATE — dropped",
                    fragment.face,
                    fragment.outer_loop.len()
                );
            }
            selections.push(FragmentSelection {
                keep: false,
                flip: false,
            });
            continue;
        }
        let (sample, class, cooriented) = match results[i] {
            Some(r) => r,
            // Defensive: a planner/member bookkeeping miss falls back
            // to the per-fragment path rather than a wrong verdict.
            None => super::fragment_classify_one::sample_and_classify(&ctx, fragment),
        };
        let own_operand = fragment.face.operand;
        let other_aabbs: &[_] = match own_operand {
            Operand::A => &aabbs_b,
            Operand::B => &aabbs_a,
        };
        // Record-side witness votes (PK-494 layer 3): collected per
        // RAY-CLASSIFIED fragment, applied post-loop after per-face
        // outward-sign calibration (the stored orientation flags lie
        // on some construction paths — the ray-verdict consensus
        // recovers each face's honest sign; only minority mismatches
        // under a clear majority get overridden). Scoped to plain
        // In/Out. PROPAGATED fragments are excluded: their verdicts
        // are structural (seed ray ⊕ flip parity, component-audited),
        // not rays for the witness to check — and each vote's
        // segment-obstruction gate prices a dense surface×segment
        // scan over the multi-wrap walls, which re-created the
        // quadratic merge cost propagation had just removed
        // (measured: len-80 merge 91 s with votes on all fragments).
        // Feeding their vote copies into the sign-calibration
        // consensus would also multiply the seed's one ray verdict
        // by the component size.
        //
        // Collection is OPERAND-A-ONLY like the overrides: a B-side
        // fragment's votes reference operand-A faces, and signs for
        // those keys are never consulted (overrides are A-only and
        // read the B-face signs that A-side votes calibrate) — dead
        // tallies whose obstruction gates priced O(records²) on the
        // tool flanks (each multi-wrap flank fragment carries a vote
        // per wrap record, each vote a dense surface×segment scan).
        if own_operand == Operand::A
            && !propagated[i]
            && matches!(
                class,
                PointClassification::Inside | PointClassification::Outside
            )
        {
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
                witness_pending.push((selections.len(), class, cooriented, votes));
            }
        }

        if popkern_core::trace::trace_enabled("cellular") {
            eprintln!(
                "  [sel] {:?} outer={} sample=({:.2},{:.2},{:.2}) class={class:?} co={cooriented:?}",
                fragment.face,
                fragment.outer_loop.len(),
                sample.x,
                sample.y,
                sample.z
            );
        }
        selections.push(selection_for(op, own_operand, class, cooriented));
    }

    // Post-pass: calibrate per-face outward signs from the ray-verdict
    // consensus and override minority mismatches (PK-494 layer 3).
    if !witness_pending.is_empty() {
        let all: Vec<(
            PointClassification,
            Vec<crate::record_side_witness::RecordVote>,
        )> = witness_pending
            .iter()
            .map(|(_, c, _, v)| (*c, v.clone()))
            .collect();
        let signs = crate::record_side_witness::calibrate_face_signs(&all);
        for (idx, class, cooriented, votes) in &witness_pending {
            // OPERAND-A-ONLY overrides (see the collection comment):
            // tool-side fragments carry walk phantom semantics that
            // their ray verdicts encode deliberately (PK-428
            // keyed-bore class — both measured regressions were
            // B-side, and confidence does not discriminate: the
            // keyed phantom was contested while the coil-end
            // void-backed sliver reads confidently wrong at 50+
            // turns). Workpiece-side fragments have no such
            // semantics.
            if fragments[*idx].face.operand != Operand::A {
                continue;
            }
            let Some(witness) = crate::record_side_witness::calibrated_verdict(votes, &signs)
            else {
                continue;
            };
            if witness == *class {
                continue;
            }
            let fragment = &fragments[*idx];
            if popkern_core::trace::trace_enabled("cellular")
                || popkern_core::trace::trace_enabled("boolean")
            {
                eprintln!(
                    "  [record-side] {:?} ray {class:?} -> calibrated witness {witness:?}",
                    fragment.face
                );
            }
            selections[*idx] = selection_for(op, fragment.face.operand, witness, *cooriented);
        }
    }

    Ok(selections)
}
