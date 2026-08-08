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
/// clean serving fixtures when tried). A hit forces ONE re-run with
/// the offending vertices excluded from direct-apex promotion so the
/// legacy dispatch serves those corners; the retry is recorded on the
/// repair report. A failing retry keeps the first result — never
/// serve worse than the measured first pass. Opt-in via
/// [`FilletOptions::serving_self_intersection_retry`] until the
/// PK-245 arming campaign triages the pre-existing debt an always-on
/// retry exposes (9 fixtures measured 2026-08-07).
pub(crate) fn fillet_edges_internal(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    mut tracking: Option<&mut TrackingContext>,
    cancel: Option<&cancellation::CancellationToken>,
    mut reporter: Option<&mut dyn progress::ProgressReporter>,
) -> Result<BlendResult, BlendError> {
    let mut first = fillet_edges_attempt(
        brep,
        specs,
        options,
        tolerance,
        &mut tracking,
        cancel,
        &mut reporter,
        &std::collections::BTreeSet::new(),
    )?;
    if !options.serving_self_intersection_retry {
        return Ok(first);
    }
    let offending = served_direct_apex_self_intersections(&first);
    if offending.is_empty() {
        return Ok(first);
    }
    // PK-251: before downgrading the corner to legacy, attempt the
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
    let recheck_clean = all_cleared
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
    let excluded: std::collections::BTreeSet<popkern_core::topology::VertexHandle> =
        offending.iter().map(|(vertex, _)| *vertex).collect();
    match fillet_edges_attempt(
        brep,
        specs,
        options,
        tolerance,
        &mut tracking,
        cancel,
        &mut reporter,
        &excluded,
    ) {
        Ok(mut retried) => {
            retried.repair_report.serving_self_intersection_retries = offending
                .into_iter()
                .map(
                    |(vertex, detail)| crate::types::ServingSelfIntersectionRetry {
                        vertex,
                        detail,
                    },
                )
                .collect();
            // Keep the failed repair attempts visible on the served
            // report — the retried pipeline run knows nothing of them.
            retried.repair_report.serving_footprint_repairs =
                std::mem::take(&mut first.repair_report.serving_footprint_repairs);
            Ok(retried)
        }
        Err(_) => Ok(first),
    }
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
