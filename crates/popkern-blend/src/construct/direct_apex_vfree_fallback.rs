//! Measured fallback boundary for the promoted V-free repeated-apex
//! default (PK-124).
//!
//! V-free topology is the production default, but its Stage-D contract is
//! only proven for the pyramid repeated-apex class; cube-class apexes
//! interact with dispatch-time pocket completion (PK-074) and fail native
//! acceptance under it. When the default (not an explicit strict
//! `dispatch_overrides::V_FREE_REPEATED_APEX` request) fails acceptance,
//! the dispatch retries the apex-retaining path once and records the
//! fallback as a structured diagnostic, so an explicit strict request
//! still fails loudly while production inputs keep their strongest
//! passing construction.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::apex_construction::construct::ApexBuildResult;
use crate::apex_construction::constructive_builder::construct_apex_constructively;
use crate::apex_construction::constructive_builder_options::ConstructiveApexExecutionOptions;
use crate::fillet_surface::FilletSurfaceResult;

/// Whether the active V-free contract came from the promoted default
/// rather than an explicit override request.
pub(super) fn v_free_promoted_by_default() -> bool {
    crate::dispatch_overrides::V_FREE_REPEATED_APEX
        .get()
        .is_none()
}

/// Measured fallback decision for the promoted V-free default: when the
/// strict result failed structural/G0 acceptance OR the cap-outer G1
/// commercial gate, retry the apex-retaining path. The strict topology
/// contract stays available under the explicit strict override
/// (`dispatch_overrides::V_FREE_REPEATED_APEX`), which never
/// falls back. Returns the accepted retry result, or `None` when the
/// strict result stands (or the retry also declined).
#[allow(clippy::too_many_arguments)]
pub(super) fn resolve_promoted_v_free_decline(
    brep: &mut BRep,
    snapshot: &BRep,
    plan: &crate::apex_construction::plan::ApexPlan,
    surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
    execution_options: ConstructiveApexExecutionOptions,
    constructive_err: Option<&str>,
    result: &ApexBuildResult,
    acceptance_error: impl Fn(&BRep, &ApexBuildResult) -> Option<String>,
) -> Option<ApexBuildResult> {
    if !execution_options.require_v_free_repeated_apex_for(plan.n_fillets())
        || !v_free_promoted_by_default()
    {
        return None;
    }
    if popkern_core::trace::trace_enabled("cap-g1") {
        eprintln!("[capg1] resolve_promoted_v_free_decline: constructive_err={constructive_err:?}");
    }
    let g1_declined = constructive_err.is_none();
    let decline_reason = match constructive_err {
        Some(err) => err.to_string(),
        None => {
            super::direct_apex_acceptance::native_gregory_cap_outer_g1_diagnostic(
                brep, result, 5.0e-3, 5.0e-2,
            )?;
            "cap-outer G1 exceeded the commercial gate under the promoted V-free default"
                .to_string()
        }
    };
    // When only the G1 gate declined, the strict result is still a valid
    // fallback target: keep whichever construction actually measures
    // better instead of unconditionally preferring the retry.
    let strict_state = g1_declined.then(|| brep.clone());
    let retry = retry_apex_retaining_after_v_free_decline(
        brep,
        snapshot,
        plan,
        surfaces,
        tolerance,
        execution_options,
        &decline_reason,
        acceptance_error,
    );
    match (retry.as_ref(), strict_state) {
        (Some(retry_result), Some(strict_state)) => {
            let retry_also_fails_g1 =
                super::direct_apex_acceptance::native_gregory_cap_outer_g1_diagnostic(
                    brep,
                    retry_result,
                    5.0e-3,
                    5.0e-2,
                )
                .is_some();
            if retry_also_fails_g1 {
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace] apex-retaining retry also fails the cap-outer G1 \
                         gate; keeping the strict V-free result"
                    );
                }
                *brep = strict_state;
                return None;
            }
        }
        (None, Some(strict_state)) => {
            // The retry declined or ERRORED with `brep` restored to
            // the pre-construction snapshot; the caller keeps the
            // strict RESULT, so the strict BODY must be restored too
            // — leaking the snapshot here served a carve-era body
            // with a committed-cap report (measured pyramid3_uniform:
            // 10 faces, plane cap, 2/3 seams behind Accepted/Coons).
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace] apex-retaining retry unavailable; keeping the \
                     strict V-free result (G1 diagnostic recorded)"
                );
            }
            *brep = strict_state;
            return None;
        }
        _ => {}
    }
    retry
}

/// Retry the constructive build with the apex-retaining options after a
/// V-free decline. Returns the accepted retry result, or `None` (with
/// `brep` restored to `snapshot`) when the retry also declines.
#[allow(clippy::too_many_arguments)]
pub(super) fn retry_apex_retaining_after_v_free_decline(
    brep: &mut BRep,
    snapshot: &BRep,
    plan: &crate::apex_construction::plan::ApexPlan,
    surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
    execution_options: ConstructiveApexExecutionOptions,
    v_free_error: &str,
    acceptance_error: impl Fn(&BRep, &ApexBuildResult) -> Option<String>,
) -> Option<ApexBuildResult> {
    *brep = snapshot.clone();
    let retry_options = execution_options.with_v_free_repeated_apex_disabled();
    let mut retry_result =
        match construct_apex_constructively(brep, plan, surfaces, tolerance, None, retry_options) {
            Ok(result) => result,
            Err(err) => {
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace] apex-retaining retry after V-free decline \
                         errored: {err}"
                    );
                }
                *brep = snapshot.clone();
                return None;
            }
        };
    if let Some(retry_error) = acceptance_error(brep, &retry_result) {
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace] apex-retaining retry after V-free decline also failed: \
                 {retry_error}"
            );
        }
        *brep = snapshot.clone();
        return None;
    }
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] V-free default declined ({v_free_error}); apex-retaining \
             retry accepted"
        );
    }
    retry_result
        .diagnostics
        .push(crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "direct_apex.v_free_promotion".to_string(),
                code: "direct_apex.v_free_promotion.apex_retaining_fallback".to_string(),
                detail_code: Some(v_free_error.to_string()),
                sector_name: "direct_apex".to_string(),
                edge: None,
                support_half_edge: None,
                support_face: None,
                residual: None,
                outer_side: None,
                sample_parameter: None,
                tolerance: None,
                shared_fact_ids: Vec::new(),
            },
        ));
    Some(retry_result)
}

/// PK-153/PK-167: retire the V-free vestigial old shell (flaps,
/// slivers, mirror winding, unglued membranes) and weld the exposed
/// seams. Self-guarding no-op off the V-free path; transactional
/// inside the consolidation itself. Since the constructive in-place
/// path became primary, the winding reconciliation here is a
/// VALIDATOR BACKSTOP: firing means an upstream constructive
/// invariant broke, so the caller surfaces a structured diagnostic
/// (the aggregate fallback path may legitimately still use it).
pub(super) fn consolidate_served_shell(
    brep: &mut popkern_core::topology::BRep,
) -> (bool, Option<crate::types::ApexConstructionDiagnostic>) {
    let consolidation =
        crate::apex_construction::v_free_shell_consolidation::consolidate_vestigial_support_shell(
            brep,
        );
    if consolidation.applied() && crate::direct_apex_trace::enabled() {
        eprintln!("[direct-apex-trace] vestigial shell consolidation: {consolidation:?}");
    }
    let did_work = consolidation.applied();
    if consolidation.reversed_faces == 0 && consolidation.rebuilt_rings == 0 {
        return (did_work, None);
    }
    (
        did_work,
        Some(crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "v_free_consolidation".to_string(),
                code: "v_free_consolidation.winding_backstop_fired".to_string(),
                detail_code: Some(format!(
                    "reversed_faces={};rebuilt_rings={}",
                    consolidation.reversed_faces, consolidation.rebuilt_rings
                )),
                sector_name: "direct_apex".to_string(),
                edge: None,
                support_half_edge: None,
                support_face: None,
                residual: None,
                outer_side: None,
                sample_parameter: None,
                tolerance: None,
                shared_fact_ids: Vec::new(),
            },
        )),
    )
}
