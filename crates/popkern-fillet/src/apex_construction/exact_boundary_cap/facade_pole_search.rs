//! PK-255 pole-placement search (design:
//! docs/architecture/exact-boundary-pole-placement.md).
//!
//! The exact-boundary pole defaults to the ring shared-corner
//! centroid, which at high radius asymmetry lands ON strip
//! end-section planes and fans sector sheets into a corner-anchored
//! cap×strip crossing no interior-push repair can clear. Plane-margin
//! heuristics were measured non-monotone across the fixture×corner
//! matrix, so placement is a SEARCH accepted by measurement on the
//! solved sheets: candidate poles along the blend-ward ray, each
//! evaluated by the validator's own crossing verdict plus the
//! cap-outer G0/G1 bars, smallest displacement wins, explicit
//! centroid fallback otherwise. Kill-switch: `POPKERN_POLE_SEARCH=0`.

use popkern_core::topology::{BRep, FaceHandle};

use super::cap_quality_postcondition::validate_candidate_region_self_intersection;
use super::facade::{
    ExactBoundaryCapCandidateResult, ExactBoundaryCapCandidateStatus,
    run_exact_boundary_cap_candidate_attempt,
};
use crate::apex_construction::error::ApexBuildError;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;

/// Candidate stations along the blend-ward ray, as fractions of the
/// search scale (the largest fillet-end arc chord). 0.15 was the
/// measured winner on cube [0.05, 0.1, 0.2]; the ladder brackets it.
const CANDIDATE_FRACTIONS: [f64; 4] = [0.075, 0.15, 0.225, 0.3];
/// Commercial cap-outer bars a winning candidate must hold.
const MAX_CAP_OUTER_G0: f64 = 5.0e-3;
const MAX_CAP_OUTER_G1: f64 = 5.0e-2;

pub(super) fn run_with_pole_search(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    sectors: &[GregoryCapSector],
) -> Result<ExactBoundaryCapCandidateResult, ApexBuildError> {
    let entry_snapshot = brep.clone();
    // Every attempt below starts from `entry_snapshot`, so the
    // fast-validation baseline is computed ONCE here, scoped to the
    // cap's one-ring (PK-273): the whole-model sweep is priced by
    // committed Gregory count and was previously recomputed
    // identically per attempt.
    let baseline_validation_errors =
        brep.validate_region_collect_fast(&brep.validation_region(&[cap_face]));
    let baseline = run_exact_boundary_cap_candidate_attempt(
        brep,
        cap_face,
        plan,
        sectors,
        &baseline_validation_errors,
        None,
    )?;
    if !search_enabled()
        || baseline.status != ExactBoundaryCapCandidateStatus::Accepted
        || committed_region_is_clean(brep, &baseline.committed_faces)
    {
        return Ok(baseline);
    }
    let Some((pole, Some((direction, scale)))) = baseline.pole else {
        return Ok(baseline);
    };
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] pole search: baseline commit is crossing-flagged; \
             trying {} candidates along ({:.3},{:.3},{:.3}) scale {scale:.4}",
            CANDIDATE_FRACTIONS.len(),
            direction.x,
            direction.y,
            direction.z,
        );
    }
    for fraction in CANDIDATE_FRACTIONS {
        *brep = entry_snapshot.clone();
        let candidate_pole = pole + direction * (scale * fraction);
        let candidate = run_exact_boundary_cap_candidate_attempt(
            brep,
            cap_face,
            plan,
            sectors,
            &baseline_validation_errors,
            Some(candidate_pole),
        )?;
        if candidate.status == ExactBoundaryCapCandidateStatus::Accepted
            && committed_region_is_clean(brep, &candidate.committed_faces)
            && committed_quality_holds(brep, &candidate.committed_faces)
        {
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace] pole search: candidate at fraction {fraction} \
                     serves crossing-free within the commercial bars"
                );
            }
            return Ok(candidate);
        }
    }
    // Centroid fallback: reproduce the baseline deterministically.
    *brep = entry_snapshot;
    run_exact_boundary_cap_candidate_attempt(
        brep,
        cap_face,
        plan,
        sectors,
        &baseline_validation_errors,
        None,
    )
}

fn search_enabled() -> bool {
    !std::env::var("POPKERN_POLE_SEARCH").is_ok_and(|value| value == "0")
}

fn committed_region_is_clean(brep: &BRep, committed_faces: &[FaceHandle]) -> bool {
    !committed_faces.is_empty()
        && validate_candidate_region_self_intersection(brep, committed_faces).is_ok()
}

fn committed_quality_holds(brep: &BRep, committed_faces: &[FaceHandle]) -> bool {
    crate::construct::cap_outer_acceptance::cap_support_acceptance_error(
        brep,
        committed_faces,
        MAX_CAP_OUTER_G0,
        MAX_CAP_OUTER_G1,
    )
    .is_none()
}
