//! PK-533: the served-body gate and the serving retry.
//!
//! The territory gate runs inside each pipeline attempt, so a gate
//! decline once pre-empted the PK-250/PK-260 exclusion retry that
//! serves the same corner clean. A steep triangular pyramid apex
//! declined under default options: its first pass falls back from
//! direct-apex to the vertex-blend builder on direct-apex anchors,
//! and that body has 24 edges off their faces (worst 9.2e-2). A
//! gate-declined first pass now takes the exclusion retry for its
//! direct-apex corners.

use popkern_blend::error::{BlendError, BlendFailure};
use popkern_blend::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_blend::{BlendResult, fillet_edges};
use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Fillet the apex of a steep triangular pyramid (height 3 over a
/// unit base) with mixed radii.
fn steep_apex(options: &FilletOptions) -> Result<BlendResult, BlendError> {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(3, 1.0, 3.0, tol).expect("pyramid fixture builds");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(vertex, _)| vertex)
        .expect("pyramid has an apex");
    let specs: Vec<FilletSpec> = pyramid
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .zip([0.04, 0.06, 0.05])
        .map(|(edge, r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    assert_eq!(specs.len(), 3);
    fillet_edges(&pyramid, &specs, options, tol, None)
}

/// Default options serve the apex check-9 clean through the retry,
/// and the rejected first pass is recorded with the gate's reason.
#[test]
fn gate_declined_first_pass_is_served_by_the_exclusion_retry() {
    let result = steep_apex(&FilletOptions::default()).expect("the steep apex serves");
    let errors = result.brep.validate_all();
    assert!(
        !errors
            .iter()
            .any(|e| matches!(e, ValidationError::CurveSurfaceDeviation { .. })),
        "the served apex must hold check 9"
    );
    let retries = &result.repair_report.serving_self_intersection_retries;
    assert_eq!(retries.len(), 1, "one corner retried, got {retries:?}");
    assert!(
        retries[0]
            .detail
            .contains("DirectApexFallbackToVertexBlend")
            && retries[0].detail.contains("off their faces"),
        "the record must name the rejected route and the gate, got {}",
        retries[0].detail
    );
    assert!(
        result.corner_patches.iter().all(|patch| patch.setback_cap),
        "the retry serves the apex through setback"
    );
}

/// With the retry off, the first pass's inconsistency is the typed
/// decline: the rescue is disabled, never bypassed.
#[test]
fn without_the_retry_the_first_pass_declines_typed() {
    let options = FilletOptions {
        serving_self_intersection_retry: false,
        ..FilletOptions::default()
    };
    let Err(BlendError::PartialFailure { failures, .. }) = steep_apex(&options) else {
        panic!("the first pass must not serve");
    };
    assert!(
        failures.iter().any(|(_, failure)| matches!(
            failure,
            BlendFailure::InconsistentGeometry { reason, location: Some(_) }
                if reason.contains("off their faces")
        )),
        "the decline must be the typed territory inconsistency, got {failures:?}"
    );
}
