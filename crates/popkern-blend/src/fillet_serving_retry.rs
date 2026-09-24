//! PK-250 serving-boundary self-intersection retry for the fillet
//! pipeline — the wrapper every public fillet entry point calls.
//!
//! See [`fillet_edges_internal`] for the contract; the pipeline body
//! itself is [`crate::fillet_pipeline::fillet_edges_attempt`].

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BlendError;
use crate::fillet_pipeline::fillet_edges_attempt;
use crate::types::{BlendResult, FilletOptions, FilletSpec};
use crate::{cancellation, progress};

/// Entry for [`fillet_edges`] + cancellable/progress variants: one
/// pipeline attempt, then the PK-250 serving-boundary check — each
/// served direct-apex corner region must be free of face-face
/// self-intersection (measured AFTER every tail pass; a mid-pipeline
/// check is premature by the PK-165/PK-220 lessons and flipped 20
/// clean serving fixtures when tried) — and, since PK-260, no corner
/// may carry a recorded serving-shape decline (the mid-pipeline
/// fallback serves a degraded vertex-blend body built on
/// direct-apex-convention anchors). Either trigger forces ONE re-run
/// with the affected vertices excluded from direct-apex promotion so
/// the vertex-blend dispatch serves those corners from clean vertex-blend
/// anchors; the retry is recorded on the repair report. A failing
/// retry keeps the first result — never serve worse than the
/// measured first pass. Default ON since the PK-245 arming campaign
/// (per-fixture triage of the pre-existing debt, 2026-08-07).
pub(crate) fn fillet_edges_internal(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    mut tracking: Option<&mut TrackingContext>,
    cancel: Option<&cancellation::CancellationToken>,
    mut reporter: Option<&mut dyn progress::ProgressReporter>,
) -> Result<BlendResult, BlendError> {
    let mut gated = Vec::new();
    let first = fillet_edges_attempt(
        brep,
        specs,
        options,
        tolerance,
        &mut tracking,
        cancel,
        &mut reporter,
        &std::collections::BTreeSet::new(),
        &mut gated,
    );
    let mut first = match first {
        Ok(first) => first,
        // PK-533: a first pass the served-body gate declined still
        // gets the exclusion re-run for its direct-apex corners. The
        // gate fires inside the attempt, so without this the decline
        // pre-empted the very re-run that serves those corners clean.
        Err(error) if options.serving_self_intersection_retry && !gated.is_empty() => {
            return retry_gated_decline(
                brep,
                specs,
                options,
                tolerance,
                &mut tracking,
                cancel,
                &mut reporter,
                error,
                gated,
            );
        }
        Err(error) => return Err(error),
    };
    if !options.serving_self_intersection_retry {
        return Ok(first);
    }
    let offending = served_direct_apex_self_intersections(&first);
    // PK-260: serving-shape declines trigger the same exclusion
    // re-run. The mid-pipeline fallback runs the vertex-blend builder on
    // direct-apex-convention anchors and serves a degraded body
    // (measured 8.5e-2 edge/surface deviations on the oblique
    // base-corner class); excluding the corner BEFORE anchor splitting
    // re-places perp-walk anchors and the vertex-blend dispatch serves clean.
    let declined = serving_shape_declined_corners(&first);
    if offending.is_empty() && declined.is_empty() {
        return Ok(first);
    }
    // PK-251: before downgrading the corner to vertex-blend, attempt the
    // in-place footprint repair — validator-driven interior
    // control-point pushes on the served cap. The repair runs on a
    // CLONE and is adopted only when every offending region clears,
    // so the "never serve worse than the measured first pass"
    // contract is preserved; anything else falls through to the
    // PK-250 exclusion retry unchanged (with the attempt recorded).
    let mut repaired_body = first.brep.clone();
    let mut repairs = Vec::new();
    for (vertex, _) in &offending {
        let Some(built_faces) = first
            .repair_report
            .apex_routes
            .iter()
            .find(|route| route.vertex == *vertex)
            .map(|route| route.built_faces.clone())
        else {
            continue;
        };
        repairs.push(
            crate::fillet_serving_footprint_repair::repair_served_direct_apex_footprint(
                &mut repaired_body,
                *vertex,
                &built_faces,
            ),
        );
    }
    let all_cleared = !repairs.is_empty() && repairs.iter().all(|repair| repair.cleared);
    // A serving-shape decline cannot be footprint-repaired (the
    // declined corner has no served direct cap) — it always takes the
    // exclusion re-run.
    let recheck_clean = declined.is_empty()
        && all_cleared
        && first
            .repair_report
            .apex_routes
            .iter()
            .filter(|route| {
                matches!(route.route, crate::types::ApexConstructionRoute::DirectApex)
                    && !route.built_faces.is_empty()
            })
            .all(|route| {
                crate::apex_construction::exact_boundary_cap::cap_quality_postcondition::
                    validate_candidate_region_self_intersection(&repaired_body, &route.built_faces)
                    .is_ok()
            });
    first.repair_report.serving_footprint_repairs = repairs;
    if recheck_clean {
        first.brep = repaired_body;
        // The repair moved cap interior control points, so the served
        // mid-seam numbers recorded earlier are stale — re-measure on
        // the repaired body (quality stays honestly reported).
        crate::apex_construction::exact_boundary_cap::cap_quality_postcondition::
            refresh_served_mid_seam_diagnostics(&first.brep, &mut first.repair_report);
        return Ok(first);
    }
    let excluded: std::collections::BTreeSet<popkern_core::topology::VertexHandle> = offending
        .iter()
        .map(|(vertex, _)| *vertex)
        .chain(declined.iter().map(|(vertex, _)| *vertex))
        .collect();
    match fillet_edges_attempt(
        brep,
        specs,
        options,
        tolerance,
        &mut tracking,
        cancel,
        &mut reporter,
        &excluded,
        &mut Vec::new(),
    ) {
        Ok(mut retried) => {
            // PK-258 keep-better: adopt the retried body only when it is
            // the better SERVE. Crossing-freedom is the hard floor (a
            // self-intersecting body cannot pass the PK-165 serving
            // gate, so no quality number buys it back — measured: the
            // recordless modify fixture shipped a 9v3x7v5 crossing when
            // bars outranked validity); served cap-outer quality decides
            // between validity-EQUAL candidates, where the blind
            // pre-contract adoption traded a g1 0.755 vertex-blend fan for a
            // served direct cap without ever measuring either. Ties
            // adopt the retry. An unmeasurable candidate (no served
            // corner faces) fails the bars and ranks worst on quality —
            // quality claims need evidence.
            // A serving-shape decline counts against the candidate's
            // validity floor exactly like a crossing: the declined
            // corner's construction failed its postconditions and the
            // mid-pipeline fallback body is degraded.
            let first_score = serve_score(
                &first,
                &excluded,
                !offending.is_empty() || !declined.is_empty(),
            );
            let retried_offending = served_direct_apex_self_intersections(&retried);
            let retried_declined = serving_shape_declined_corners(&retried);
            let retried_score = serve_score(
                &retried,
                &excluded,
                !retried_offending.is_empty() || !retried_declined.is_empty(),
            );
            let adopt_retry = retried_score >= first_score;
            let decision_detail = format!(
                "; keep-better: first (bars={}, crossing_free={}, g0={:.3e}, g1={:.3e}) vs \
                 retry (bars={}, crossing_free={}, g0={:.3e}, g1={:.3e}) -> {}",
                first_score.meets_bars,
                first_score.crossing_free,
                first_score.worst_g0,
                first_score.worst_g1,
                retried_score.meets_bars,
                retried_score.crossing_free,
                retried_score.worst_g0,
                retried_score.worst_g1,
                if adopt_retry { "retry" } else { "first" }
            );
            let retries: Vec<crate::types::ServingSelfIntersectionRetry> = offending
                .into_iter()
                .chain(declined)
                .map(
                    |(vertex, detail)| crate::types::ServingSelfIntersectionRetry {
                        vertex,
                        detail: format!("{detail}{decision_detail}"),
                    },
                )
                .collect();
            if !adopt_retry {
                first.repair_report.serving_self_intersection_retries = retries;
                return Ok(first);
            }
            retried.repair_report.serving_self_intersection_retries = retries;
            // Keep the failed repair attempts visible on the served
            // report — the retried pipeline run knows nothing of them.
            retried.repair_report.serving_footprint_repairs =
                std::mem::take(&mut first.repair_report.serving_footprint_repairs);
            Ok(retried)
        }
        Err(_) => Ok(first),
    }
}

/// The exclusion re-run for a first pass the served-body gate
/// declined (PK-533). The first pass has no serve to keep, so the
/// re-run is adopted whenever it passes every gate itself, with each
/// excluded corner recorded on the report under the gate's reason.
/// Otherwise the first pass's typed decline stands.
#[allow(clippy::too_many_arguments)]
fn retry_gated_decline(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: &mut Option<&mut TrackingContext>,
    cancel: Option<&cancellation::CancellationToken>,
    reporter: &mut Option<&mut dyn progress::ProgressReporter>,
    first_decline: BlendError,
    gated: Vec<(popkern_core::topology::VertexHandle, String)>,
) -> Result<BlendResult, BlendError> {
    let excluded = gated.iter().map(|(vertex, _)| *vertex).collect();
    match fillet_edges_attempt(
        brep,
        specs,
        options,
        tolerance,
        tracking,
        cancel,
        reporter,
        &excluded,
        &mut Vec::new(),
    ) {
        Ok(mut retried) => {
            retried.repair_report.serving_self_intersection_retries = gated
                .into_iter()
                .map(
                    |(vertex, detail)| crate::types::ServingSelfIntersectionRetry {
                        vertex,
                        detail: format!("{detail}; exclusion retry served"),
                    },
                )
                .collect();
            Ok(retried)
        }
        Err(_) => Err(first_decline),
    }
}

/// Keep-better ranking for a served candidate over the retried corner
/// set: crossing-freedom first (the hard serving-validity floor),
/// commercial cap-outer bars second, then measured quality. Ordered
/// so `>=` means "at least as good".
#[derive(Clone, Copy)]
struct ServeScore {
    meets_bars: bool,
    crossing_free: bool,
    worst_g0: f64,
    worst_g1: f64,
}

impl ServeScore {
    const MAX_CAP_OUTER_G0: f64 = 5.0e-3;
    const MAX_CAP_OUTER_G1: f64 = 5.0e-2;
}

impl PartialEq for ServeScore {
    fn eq(&self, other: &Self) -> bool {
        self.partial_cmp(other) == Some(std::cmp::Ordering::Equal)
    }
}

impl PartialOrd for ServeScore {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        (
            self.crossing_free,
            self.meets_bars,
            std::cmp::Reverse(ordered_float(self.worst_g1)),
            std::cmp::Reverse(ordered_float(self.worst_g0)),
        )
            .partial_cmp(&(
                other.crossing_free,
                other.meets_bars,
                std::cmp::Reverse(ordered_float(other.worst_g1)),
                std::cmp::Reverse(ordered_float(other.worst_g0)),
            ))
    }
}

fn ordered_float(value: f64) -> u64 {
    // Non-negative continuity measures order correctly through their
    // IEEE bit patterns; NaN (unmeasured) ranks worst of all.
    if value.is_nan() {
        u64::MAX
    } else {
        value.to_bits()
    }
}

/// Measure a candidate's serve over the retried corner set: worst
/// served cap-outer continuity across the corners' recorded cap faces
/// (direct AND vertex-blend routes record them), plus crossing state. A
/// corner with no measurable seam contributes NaN, which fails the
/// bars and ranks worst — quality claims need evidence.
fn serve_score(
    result: &BlendResult,
    corners: &std::collections::BTreeSet<popkern_core::topology::VertexHandle>,
    has_crossings: bool,
) -> ServeScore {
    let mut worst_g0 = 0.0f64;
    let mut worst_g1 = 0.0f64;
    for route in result
        .repair_report
        .apex_routes
        .iter()
        .filter(|route| corners.contains(&route.vertex))
    {
        match crate::construct::cap_outer_acceptance::worst_cap_support_continuity(
            &result.brep,
            &route.built_faces,
        ) {
            Some((g0, g1)) => {
                worst_g0 = worst_g0.max(g0);
                worst_g1 = worst_g1.max(g1);
            }
            None => {
                worst_g0 = f64::NAN;
                worst_g1 = f64::NAN;
                break;
            }
        }
    }
    ServeScore {
        meets_bars: worst_g0 <= ServeScore::MAX_CAP_OUTER_G0
            && worst_g1 <= ServeScore::MAX_CAP_OUTER_G1,
        crossing_free: !has_crossings,
        worst_g0,
        worst_g1,
    }
}

/// Corners whose direct-apex construction declined at the
/// serving-shape gate (PK-260), with the recorded decline reason. The
/// mid-pipeline fallback serves a degraded vertex-blend body (direct-apex
/// anchor convention); the exclusion re-run serves clean.
fn serving_shape_declined_corners(
    result: &BlendResult,
) -> Vec<(popkern_core::topology::VertexHandle, String)> {
    result
        .repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            route.route == crate::types::ApexConstructionRoute::DirectApexFallbackToVertexBlend
        })
        .filter_map(|route| {
            route
                .fallback_reason
                .as_ref()
                .filter(|reason| reason.contains(crate::types::SERVING_SHAPE_DECLINE_PREFIX))
                .map(|reason| (route.vertex, format!("serving-shape decline: {reason}")))
        })
        .collect()
}

/// Served direct-apex corners whose region self-intersects, with the
/// validator's face-pair detail.
fn served_direct_apex_self_intersections(
    result: &BlendResult,
) -> Vec<(popkern_core::topology::VertexHandle, String)> {
    result
        .repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            matches!(route.route, crate::types::ApexConstructionRoute::DirectApex)
                && !route.built_faces.is_empty()
        })
        .filter_map(|route| {
            crate::apex_construction::exact_boundary_cap::cap_quality_postcondition::
                validate_candidate_region_self_intersection(&result.brep, &route.built_faces)
                .err()
                .map(|detail| (route.vertex, detail))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::ServeScore;

    fn score(crossing_free: bool, worst_g0: f64, worst_g1: f64) -> ServeScore {
        ServeScore {
            meets_bars: worst_g0 <= ServeScore::MAX_CAP_OUTER_G0
                && worst_g1 <= ServeScore::MAX_CAP_OUTER_G1,
            crossing_free,
            worst_g0,
            worst_g1,
        }
    }

    /// PK-258 keep-better ranking: crossing-freedom is the hard
    /// serving-validity floor, commercial bars second, measured
    /// quality third, and unmeasured (NaN) ranks worst of all.
    #[test]
    fn keep_better_ranking_orders_validity_then_bars_then_quality() {
        // A crossing-free serve beats a crossing-flagged one no matter
        // the quality gap (the flagged body cannot pass the serving
        // gate) — the 0.755 vertex-blend fan vs 1e-4 flagged direct case.
        assert!(score(true, 1.0e-2, 0.755) > score(false, 1.0e-8, 1.0e-4));
        // Between validity-equal candidates the commercial bars, then
        // the measured quality, decide — the quality-blind adoption
        // this contract replaces.
        assert!(score(true, 1.0e-6, 1.0e-4) > score(true, 1.0e-2, 0.755));
        assert!(score(false, 1.0e-6, 1.0e-4) > score(false, 1.0e-6, 2.0e-4));
        // Unmeasured quality (no served corner faces) ranks worst
        // among validity-equal candidates: quality claims need
        // evidence.
        assert!(score(true, 1.0e-2, 0.755) > score(true, f64::NAN, f64::NAN));
        // Ties adopt the retry via `>=`.
        assert!(score(true, 1.0e-6, 1.0e-4) >= score(true, 1.0e-6, 1.0e-4));
    }
}
