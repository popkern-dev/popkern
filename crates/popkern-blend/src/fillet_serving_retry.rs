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

#[path = "fillet_serving_retry_probes.rs"]
mod probes;
use probes::{
    served_direct_apex_creases, served_direct_apex_self_intersections,
    served_direct_apex_unbounded_trims, serving_shape_declined_corners,
};

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
fn fillet_edges_serving(
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
    let mut declined = serving_shape_declined_corners(&first);
    declined.extend(served_direct_apex_unbounded_trims(&first));
    declined.extend(served_direct_apex_creases(&first));
    if offending.is_empty() && declined.is_empty() {
        return Ok(first);
    }
    // PK-251's in-place footprint repair used to run here. Retired
    // (PK-823): across 82 cube and pyramid corners it cleared its own
    // region on 12, and on every one the recheck refused the repaired
    // body and this exclusion re-run served the corner anyway — its up to
    // 64 rounds of pair intersection were always discarded.
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
            let mut retried_declined = serving_shape_declined_corners(&retried);
            retried_declined.extend(served_direct_apex_unbounded_trims(&retried));
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

/// [`fillet_edges_serving`] plus the result's setback-corner
/// compatibility record (PK-575): for every committed setback cap, how
/// well the two neighbours its boundary turns between agree at each
/// corner, against the commercial gates.
pub(crate) fn fillet_edges_internal(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: Option<&cancellation::CancellationToken>,
    reporter: Option<&mut dyn progress::ProgressReporter>,
) -> Result<BlendResult, BlendError> {
    let filleted: Vec<_> = specs.iter().map(|s| s.edge).collect();
    options
        .corner_setback
        .validate(brep, &filleted, &options.vertex_mode, options.corner_cap)?;
    let mut result =
        fillet_edges_serving(brep, specs, options, tolerance, tracking, cancel, reporter)?;
    if options.corner_cap.construction == crate::CornerConstruction::RollingBall {
        // PK-861: eligible corners trade their setback cap for the ball.
        crate::rolling_ball_corner::settle_rolling_ball_corners(
            &mut result,
            &tolerance,
            options.corner_cap.strict,
        )
        .map_err(|(vertex, reason)| {
            let edge = specs.first().map(|s| s.edge).unwrap_or_default();
            BlendError::partial_failure(
                specs.len(),
                brep.clone(),
                vec![(
                    edge,
                    crate::error::BlendFailure::CornerCapUnavailable {
                        vertex,
                        construction: crate::CornerConstruction::RollingBall,
                        reason,
                        location: brep.vertices.get(vertex).map(|v| v.position),
                    },
                )],
            )
        })?;
    }
    let caps: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.setback_cap)
        .map(|p| p.face)
        .collect();
    if !caps.is_empty() {
        result.repair_report.setback_corners =
            crate::setback_corner_compatibility::setback_corner_compatibility(
                &result.brep,
                &caps,
                &crate::blend_continuity_evidence::ContinuityGates::commercial(),
            );
        result.repair_report.setback_rims =
            crate::setback_rim_landing::setback_rim_landing(&result.brep, &caps);
    }
    Ok(result)
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
