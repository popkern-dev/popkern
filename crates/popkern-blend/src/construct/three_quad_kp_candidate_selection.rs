//! Boundary-fidelity candidate selection for the KP binonic
//! promotion (PK-249 observability).
//!
//! The promotion builds up to three candidates — the guide-derived
//! ring, the exact-injection ring, and the transactionally relocked
//! ring — and serves the strongest boundary-faithful one. The
//! selection outcome (which candidate won, and why each stronger
//! path fell short when it did) is recorded as a structured fact so
//! the extreme-asymmetry decline class stays attributable.
//!
//! Nested inside `three_quad_kp_promote` via `#[path]`.

use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GregoryPatch, RectangularGregoryBiNonic,
};

use super::{
    PROMOTED_BOUNDARY_RESIDUAL_GATE, RING_G0_GATE, RING_G1_GATE, optimize_kp_binonic_candidate,
    worst_ring_g0_gap, worst_ring_normal_sine,
};

/// The candidate the promotion will gate and commit, with its
/// provenance and the recorded selection fact.
pub(super) struct SelectedCandidate {
    pub binonic: Vec<RectangularGregoryBiNonic>,
    pub provenance: &'static str,
    pub selection_fact: Option<crate::types::ApexConstructionDiagnostic>,
    /// PK-275: the exact candidate's armed outer-fillet rung
    /// outcome, for the promotion's rung fact (None when unarmed or
    /// when the exact candidate never ran the bi-nonic optimizer).
    pub exact_rung_outcome: Option<String>,
}

/// Choose between the guide candidate, the exact-injection
/// candidate (n < 6), and the relocked candidate — strongest
/// boundary fidelity first — recording each path's outcome.
pub(super) fn select_boundary_candidate(
    guide_binonic: Vec<RectangularGregoryBiNonic>,
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
    outer_fillets: Option<&[popkern_core::geometry::surface::SurfaceGeometry]>,
    gate_fillets: &[popkern_core::geometry::surface::SurfaceGeometry],
    pre_worst: f64,
    n: usize,
) -> SelectedCandidate {
    let Some(boundaries) = input.hex_boundaries.as_deref() else {
        return SelectedCandidate {
            binonic: guide_binonic,
            provenance: "guide",
            selection_fact: None,
            exact_rung_outcome: None,
        };
    };
    let mut binonic = guide_binonic.clone();
    let mut provenance = "guide";
    let exact_outcome;
    let mut relock_outcome = String::from("not_attempted");
    let mut exact_rung_outcome: Option<String> = None;
    let exact_candidate = if n < 6 {
        match optimize_kp_binonic_candidate(input, options, outer_fillets) {
            None => {
                exact_outcome = String::from("unavailable");
                None
            }
            Some((candidate, rung_outcome)) => {
                exact_rung_outcome = rung_outcome;
                // PK-249: a fold is the severe end of the chart
                // orientation defect — attempt the transactional
                // unfold before rejecting the exact candidate; the
                // repaired ring must then pass the full gate chain
                // like any other exact candidate. PK-544: a G1 failure
                // can be a fold across a spoke, which the repair now
                // clears; a plain crease declines there as
                // `already_clean` at no cost.
                let (candidate, repair_tag) = match promotion_gate_failure(&candidate, pre_worst) {
                    Some(failure)
                        if failure.starts_with("self_intersection")
                            || failure.starts_with("g1(") =>
                    {
                        match super::super::three_quad_chart_orientation::repair_chart_orientation(
                            &candidate,
                            Some(gate_fillets),
                        ) {
                            Ok(outcome) => (Some(outcome.ring), format!("{failure}+repaired;")),
                            Err(reason) => (None, format!("{failure}+repair_declined({reason})")),
                        }
                    }
                    _ => (Some(candidate), String::new()),
                };
                match candidate {
                    None => {
                        // Repair declined: the tag is the whole
                        // outcome; fall through to the relock like
                        // any other exact-candidate failure.
                        exact_outcome = repair_tag;
                        None
                    }
                    Some(candidate) => {
                        if let Some(gate_failure) = promotion_gate_failure(&candidate, pre_worst) {
                            exact_outcome = format!("{repair_tag}{gate_failure}");
                            None
                        } else {
                            match super::super::three_quad_kp_boundary_relock::boundary_residual_for_binonic(
                                &candidate, boundaries,
                            ) {
                                Some(residual)
                                    if residual <= PROMOTED_BOUNDARY_RESIDUAL_GATE =>
                                {
                                    exact_outcome =
                                        format!("{repair_tag}accepted({residual:.3e})");
                                    Some(candidate)
                                }
                                Some(residual) => {
                                    exact_outcome = format!(
                                        "{repair_tag}boundary_residual({residual:.3e})"
                                    );
                                    None
                                }
                                None => {
                                    exact_outcome =
                                        format!("{repair_tag}residual_unmeasurable");
                                    None
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        exact_outcome = format!("ring_too_large(n={n})");
        None
    };
    if let Some(exact) = exact_candidate {
        binonic = exact;
        provenance = "exact";
    } else {
        let relock = super::super::three_quad_kp_boundary_relock::attempted_boundary_relock(
            &guide_binonic,
            boundaries,
            pre_worst,
        );
        relock_outcome = match (&relock.accepted, relock.best_residual) {
            (Some(_), Some(residual)) => format!("accepted({residual:.3e})"),
            (None, Some(residual)) => format!("boundary_residual({residual:.3e})"),
            (_, None) => String::from("no_mode_scored"),
        };
        if let Some(relocked) = relock.accepted {
            binonic = relocked;
            provenance = "relocked";
        }
    }
    let selection_fact = crate::construct::three_quad_kp_diagnostics::boundary_candidate_selection(
        provenance,
        &exact_outcome,
        &relock_outcome,
        super::super::three_quad_kp_boundary_relock::boundary_residual_for_binonic(
            &binonic, boundaries,
        ),
        n,
    );
    SelectedCandidate {
        binonic,
        provenance,
        selection_fact: Some(selection_fact),
        exact_rung_outcome,
    }
}

/// PK-249 band-clearing hook for the gate-accepted candidate:
/// attempt the transactional orientation repair when the committed
/// ring carries chart handedness-reversal bands, re-applying the
/// promotion's own pre-rung quality guard. Mutates `binonic` to the
/// repaired ring on acceptance; always returns the recorded fact
/// when the hook engaged (`None` for band-free candidates).
pub(super) fn clear_committed_candidate_bands(
    binonic: &mut Vec<RectangularGregoryBiNonic>,
    gate_fillets: &[popkern_core::geometry::surface::SurfaceGeometry],
    pre_worst: f64,
    post_worst: f64,
) -> Option<crate::types::ApexConstructionDiagnostic> {
    let chart = &crate::construct::three_quad_kp_diagnostics::chart_orientation_fact;
    let pre_repair_bands =
        super::super::three_quad_chart_orientation::binonic_reversal_band_count(binonic);
    if pre_repair_bands == 0 {
        return None;
    }
    match super::super::three_quad_chart_orientation::repair_chart_orientation(
        binonic,
        Some(gate_fillets),
    ) {
        Ok(outcome) => {
            // Re-apply the promotion's own quality guard against the
            // pre-rung baseline — the repair's internal gate only
            // bounds regression against the accepted candidate.
            if pre_worst <= 1e-12 || outcome.post_g1 <= 1.2 * pre_worst {
                let fact = chart(
                    "repaired",
                    outcome.pre_bands,
                    0,
                    outcome.pre_g1,
                    outcome.post_g1,
                    Some(outcome.post_g2),
                );
                *binonic = outcome.ring;
                Some(fact)
            } else {
                Some(chart(
                    "declined(quality_regression_vs_prerung)",
                    outcome.pre_bands,
                    outcome.pre_bands,
                    outcome.pre_g1,
                    outcome.post_g1,
                    Some(outcome.post_g2),
                ))
            }
        }
        Err(reason) => Some(chart(
            &format!("declined({reason})"),
            pre_repair_bands,
            pre_repair_bands,
            post_worst,
            post_worst,
            None,
        )),
    }
}

/// First promotion gate the candidate fails, described with its
/// measured value — `None` when every gate passes. The description
/// feeds the boundary-candidates fact so exact-injection declines
/// on the extreme-asymmetry class stay attributable.
fn promotion_gate_failure(
    candidate: &[RectangularGregoryBiNonic],
    pre_worst: f64,
) -> Option<String> {
    let binonic_gp: Vec<GregoryPatch> = candidate
        .iter()
        .cloned()
        .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
        .collect();
    if let Some(i) = binonic_gp
        .iter()
        .position(|gp| gp.self_intersects_with_tolerance(1e-6))
    {
        return Some(format!("self_intersection(sub_patch={i})"));
    }
    let g0 = worst_ring_g0_gap(&binonic_gp);
    if g0 > RING_G0_GATE {
        return Some(format!("g0({g0:.3e})"));
    }
    let post_worst = worst_ring_normal_sine(&binonic_gp);
    if post_worst > RING_G1_GATE {
        return Some(format!("g1({post_worst:.3e})"));
    }
    if pre_worst > 1e-12 && post_worst > 1.2 * pre_worst {
        return Some(format!(
            "quality_regression(pre={pre_worst:.3e},post={post_worst:.3e})"
        ));
    }
    None
}
