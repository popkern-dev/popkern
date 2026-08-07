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
    let first = fillet_edges_attempt(
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
