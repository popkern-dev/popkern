//! The setback vertex-blend rung of `insert_vertex_blend` (PK-094
//! surgery entry, PK-225 N=3 production gate, PK-574 fallback ladder).
//! See the call site for the rung ordering contract and
//! `docs/architecture/variational-corner-fit.md` for the ladder.

use super::InsertVertexBlendOutput;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

use crate::construct::setback_cap_fairness::{CAP_FAIRNESS_BOUND, CapInterior, cap_interior};
use crate::construct::setback_decline::{SetbackCapStrategy, SetbackDecline, SetbackStage};
use crate::construct::setback_vertex_topology::{
    SETBACK_TOPOLOGY_LAST, SetbackTopologyOutcome, setback_fan_commit::fan_seam_g2_failure,
    try_setback_boundary_surgery,
};
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerCapPreference, CornerConstruction,
    CornerMethod, CornerPatchInfo, CornerSubPatchProvenance, SetbackCapFairness,
    SetbackFitDiagnostic, SetbackFitMode, SetbackFitOutcome,
};
use crate::vertex_classifier::VertexBlendClassification;

#[path = "blend_insert_setback_fit_step.rs"]
mod fit_step;

/// The rung's result: `served` is `Some` when the rung decided (served
/// or hard-errored), `None` to fall through to the vertex-blend corner
/// machinery; `diagnostics` carry the variational-fit record whenever
/// the fit ran, served or not.
pub(super) type SetbackRungOutcome = (
    Option<Result<InsertVertexBlendOutput, String>>,
    Vec<ApexConstructionDiagnostic>,
);

/// One ladder step's verdict: served (or hard-errored), declined at the
/// shared carve (no construction on this carve can serve), or declined
/// with a reason the next step may overcome.
enum RungStep {
    Decided(Result<InsertVertexBlendOutput, String>),
    CarveDeclined(String),
    Declined(String),
}

/// The per-corner inputs every ladder step shares.
struct SetbackRung<'a> {
    vertex: VertexHandle,
    n_sides: usize,
    fillet_surfaces_pre: &'a [FilletSurfaceResult],
    tolerance: &'a Tolerance,
    requested_continuity: ContinuityLevel,
    classification: VertexBlendClassification,
    cap_presolve: Option<&'a crate::construct::setback_cap_presolve::SetbackCapPresolveEntry>,
    commit_snapshot: BRep,
    /// The corner's smallest fillet radius, the fairness measure's unit.
    r_min: f64,
}

/// Attempt the setback rung, each construction transactional on the
/// same carve. The order follows the corner-cap preference (PK-603):
///
/// - `Auto` / `Setback`: the KP20 cap, its G1 downgrade on a
///   seam-G2-only failure, then the variational fit (curvature-weighted
///   under a G2 request, then tangent);
/// - `Fitted`: the variational fit first, then the KP20 cap and its
///   downgrade;
/// - `Classic`: the rung is skipped.
///
/// A strict `Setback` / `Fitted` preference tries only its own
/// construction; when it declines the rung returns an empty serve
/// carrying [`ApexConstructionDiagnostic::CornerCapUnavailable`], which
/// the apply phase turns into a failure instead of a fallback.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_setback_rung(
    brep: &mut BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    fillet_surfaces_pre: &[FilletSurfaceResult],
    tolerance: &Tolerance,
    requested_continuity: ContinuityLevel,
    classification: VertexBlendClassification,
    setback_min_n: usize,
    cap_presolve: Option<&crate::construct::setback_cap_presolve::SetbackCapPresolveEntry>,
    corner_cap: CornerCapPreference,
    r_min: f64,
) -> SetbackRungOutcome {
    if corner_cap.construction == CornerConstruction::Classic
        || !(crate::construct::setback_vertex_topology::setback_topology_commit_enabled()
            && edge_indices.len() >= setback_min_n)
    {
        return (None, Vec::new());
    }
    // On any surgery error the snapshot restore inside the surgery has
    // already run; the next attempt starts from the same model.
    let rung = SetbackRung {
        vertex,
        n_sides: edge_indices.len(),
        fillet_surfaces_pre,
        tolerance,
        requested_continuity,
        classification,
        cap_presolve,
        commit_snapshot: brep.clone(),
        r_min,
    };
    let strict = corner_cap.declines_when_unavailable();
    let fitted_first = corner_cap.construction == CornerConstruction::Fitted;
    let try_kp20 = !(strict && fitted_first);
    let try_fit = !(strict && corner_cap.construction == CornerConstruction::Setback);
    // A fitted-first record is written before the KP20 cap runs.
    let mut record = SetbackFitDiagnostic {
        vertex,
        kp20_decline: match (fitted_first, try_kp20) {
            (true, false) => "not attempted: the fitted cap is strict".to_string(),
            (true, true) => "not attempted: the fitted cap is preferred".to_string(),
            _ => String::new(),
        },
        attempts: Vec::new(),
        served: None,
        fairness: None,
    };
    let mut declines: Vec<String> = Vec::new();
    let order: [bool; 2] = if fitted_first {
        [true, false]
    } else {
        [false, true]
    };
    for fit in order {
        let step = match (fit, try_fit, try_kp20) {
            (true, true, _) => {
                // Modes the KP20 step already tried (PK-612) are not rerun.
                let modes: Vec<SetbackFitMode> = rung
                    .fit_modes()
                    .iter()
                    .copied()
                    .filter(|&mode| !record.attempts.iter().any(|a| a.mode == mode))
                    .collect();
                if modes.is_empty() {
                    continue;
                }
                rung.fit(brep, &mut record, &modes)
            }
            (false, _, true) => {
                let step = rung.kp20(brep, &mut record, try_fit);
                if let RungStep::Declined(reason) | RungStep::CarveDeclined(reason) = &step {
                    record.kp20_decline = reason.clone();
                }
                step
            }
            _ => continue,
        };
        match step {
            RungStep::Decided(Ok((infos, boundary, mut diagnostics))) => {
                // A KP20 serve after declined fit attempts keeps their record.
                if record.served.is_none()
                    && (!record.attempts.is_empty() || record.fairness.is_some())
                {
                    diagnostics.push(ApexConstructionDiagnostic::SetbackFit(record));
                }
                return (Some(Ok((infos, boundary, diagnostics))), Vec::new());
            }
            RungStep::Decided(Err(e)) => return (Some(Err(e)), Vec::new()),
            RungStep::CarveDeclined(reason) => {
                declines.push(reason);
                break;
            }
            RungStep::Declined(reason) => declines.push(reason),
        }
    }
    let reason = declines.join("; ");
    record_decline(vertex, &reason);
    let mut diagnostics = vec![setback_decline_diagnostic(&reason)];
    if !record.attempts.is_empty() {
        diagnostics.push(ApexConstructionDiagnostic::SetbackFit(record));
    }
    if strict {
        diagnostics.push(ApexConstructionDiagnostic::CornerCapUnavailable {
            vertex,
            construction: corner_cap.construction,
            reason,
        });
        return (Some(Ok((Vec::new(), None, diagnostics))), Vec::new());
    }
    (None, diagnostics)
}

/// The setback rung's decline, named for the user (PK-594) and carried
/// into the corner's decline verdicts
/// (`corner_resolution::route_decline_verdicts`), first: the setback
/// cap is the corner's primary construction, and the fallbacks' own
/// verdicts alone (PnS cap-support gaps, Gregory creases) described
/// the last route tried, not why the corner has no cap.
fn setback_decline_diagnostic(reason: &str) -> ApexConstructionDiagnostic {
    let class = if reason.contains(crate::construct::setback_decline::ARC_OVERRUN) {
        "setback arcs overrun a narrow strip: a neighbouring fillet's \
         setback carries this strip's setback arc across its own \
         boundary near the apex (reduce the radius disparity at this \
         corner)"
    } else if reason.contains("fold") {
        "every setback cap folds at this corner"
    } else {
        "the setback cap declined"
    };
    ApexConstructionDiagnostic::ExactBoundaryCap(crate::types::ExactBoundaryCapDiagnostic {
        stage: "setback_rung".to_string(),
        code: "setback.declined".to_string(),
        detail_code: Some(format!("{class} ({reason})")),
        sector_name: String::new(),
        edge: None,
        support_half_edge: None,
        support_face: None,
        residual: None,
        outer_side: None,
        sample_parameter: None,
        tolerance: None,
        shared_fact_ids: Vec::new(),
    })
}

impl SetbackRung<'_> {
    fn surgery(
        &self,
        brep: &mut BRep,
        strategy: SetbackCapStrategy,
    ) -> Result<SetbackTopologyOutcome, SetbackDecline> {
        if crate::dispatch_overrides::setback_decline_forced() {
            return Err("forced decline (dispatch override)".into());
        }
        try_setback_boundary_surgery(
            brep,
            &self.commit_snapshot,
            self.vertex,
            self.fillet_surfaces_pre,
            self.tolerance,
            strategy,
            self.cap_presolve,
        )
    }

    fn info(
        &self,
        face: popkern_core::topology::FaceHandle,
        method: CornerMethod,
        achieved: ContinuityLevel,
        provenance: Option<CornerSubPatchProvenance>,
    ) -> CornerPatchInfo {
        CornerPatchInfo {
            vertex: self.vertex,
            n_sides: self.n_sides,
            method,
            achieved_continuity: achieved,
            classification: self.classification,
            face,
            bridge_fallbacks: Vec::new(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: None,
            cap_sphere: None,
            sub_patch_provenance: provenance,
            setback_cap: true,
            exact_boundary_cap: false,
        }
    }

    /// The KP20 cap. It serves both continuity levels: it is the proven
    /// watertight tiling and a G2 cap satisfies G1.
    ///
    /// Under a G2 request a fan whose seams read below G2
    /// ([`fan_seam_g2_failure`]: its dense seam sampling or the
    /// certified blend evidence, PK-596) hands the corner to the
    /// curvature-weighted fit first (`try_fit`, unless it already ran):
    /// the pyramid-4 h 0.25 r 0.05 apex fan certifies g2 2.9e-1 where
    /// the fit certifies G2 at 2.8e-2. When the fit declines, the fan
    /// serves as a measured G1 downgrade — it still beats every other
    /// fallback measured on that class (rim+verticals box: the fan's
    /// seams at g1 2.9e-3 / g2 5.4e-2 against the Gregory fallback's
    /// 27° crease; pyramid-4 h 3 r 0.05: fan g2 8.8e-2 against the
    /// tangent fit's 19) — recorded G1 on every piece, which the corner
    /// downgrade warning names. The fan is kept, not rebuilt: it does
    /// not depend on the request.
    ///
    /// A fan whose interior bends past [`CAP_FAIRNESS_BOUND`] (PK-612,
    /// see `setback_cap_fairness`) tries every fit the request admits,
    /// and the fit serves when its own interior is fairer. Otherwise
    /// the fan serves, with both readings on the record.
    fn kp20(&self, brep: &mut BRep, record: &mut SetbackFitDiagnostic, try_fit: bool) -> RungStep {
        let outcome = match self.surgery(brep, SetbackCapStrategy::Kp20) {
            Ok(outcome) => outcome,
            Err(d) if d.stage == SetbackStage::Carve => return RungStep::CarveDeclined(d.reason),
            Err(d) => return RungStep::Declined(d.reason),
        };
        let serve = |brep: &mut BRep, achieved| {
            RungStep::Decided(commit(brep, &self.commit_snapshot, &outcome, |f, _| {
                self.info(f, CornerMethod::Gregory, achieved, None)
            }))
        };
        // A KP20 (G2-interior) cap on G1 sweeps yields a G1 blend.
        let g2_failure = (self.requested_continuity == ContinuityLevel::G2)
            .then(|| fan_seam_g2_failure(brep, &outcome.cap_faces))
            .flatten();
        let achieved = if self.requested_continuity == ContinuityLevel::G2 && g2_failure.is_none() {
            ContinuityLevel::G2
        } else {
            ContinuityLevel::G1
        };
        let fan = if crate::dispatch_overrides::cap_fairness_ladder_disabled() {
            CapInterior {
                overshoot: 0.0,
                regularity: f64::INFINITY,
            }
        } else {
            cap_interior(brep, &outcome.cap_faces, self.r_min)
        };
        let (overshoot, unfair) = (fan.overshoot, fan.unfair());
        if !unfair && g2_failure.is_none() {
            return serve(brep, achieved);
        }
        let mut reasons = Vec::new();
        if let Some(reason) = &g2_failure {
            reasons.push(format!("G2 not reached, G1 downgrade held: {reason}"));
        }
        if fan.near_fold() {
            reasons.push(format!(
                "fan nearly folds (interior regularity {:.3e} under {:.0e})",
                fan.regularity,
                crate::construct::setback_cap_fairness::CAP_MIN_REGULARITY
            ));
        }
        if unfair {
            reasons.push(format!(
                "fan interior curvature {overshoot:.2} / r_min over the fairness bound \
                 {CAP_FAIRNESS_BOUND}"
            ));
            record.fairness = Some(SetbackCapFairness {
                kp20: overshoot,
                fit: None,
                bound: CAP_FAIRNESS_BOUND,
            });
        }
        record.kp20_decline = reasons.join("; ");
        // An unfair fan tries every fit the request admits; a fair fan
        // below G2 tries the curvature-weighted fit for a G2 cap.
        let ran = |mode| record.attempts.iter().any(|a| a.mode == mode);
        let modes: Vec<SetbackFitMode> = if unfair {
            self.fit_modes().to_vec()
        } else {
            vec![SetbackFitMode::CurvatureWeighted]
        }
        .into_iter()
        .filter(|&mode| !ran(mode))
        .collect();
        if try_fit && !modes.is_empty() {
            let fan_model = std::mem::replace(brep, self.commit_snapshot.clone());
            if let RungStep::Decided(Ok((infos, boundary, mut diagnostics))) =
                self.fit(brep, record, &modes)
            {
                if !unfair {
                    return RungStep::Decided(Ok((infos, boundary, diagnostics)));
                }
                let faces: Vec<_> = infos.iter().map(|info| info.face).collect();
                let fit = cap_interior(brep, &faces, self.r_min).overshoot;
                if let Some(fairness) = record.fairness.as_mut() {
                    fairness.fit = Some(fit);
                }
                if fan.near_fold() || fit < overshoot {
                    for diagnostic in &mut diagnostics {
                        if let ApexConstructionDiagnostic::SetbackFit(served) = diagnostic {
                            served.fairness = record.fairness;
                        }
                    }
                    return RungStep::Decided(Ok((infos, boundary, diagnostics)));
                }
                // The fit is no fairer: the fan serves.
                record.served = None;
                if let Some(last) = record.attempts.last_mut()
                    && let SetbackFitOutcome::Served(certified) = last.outcome
                {
                    last.outcome = SetbackFitOutcome::Rejected {
                        certified,
                        reason: format!(
                            "interior curvature {fit:.2} / r_min, no fairer than the fan's \
                             {overshoot:.2}"
                        ),
                    };
                }
            }
            *brep = fan_model;
        }
        if fan.near_fold() {
            *brep = self.commit_snapshot.clone();
            return RungStep::Declined(record.kp20_decline.clone());
        }
        serve(brep, achieved)
    }
}

/// Validate a served KP20 surgery and build its corner records; a
/// validation failure restores the snapshot and is a hard error.
fn commit(
    brep: &mut BRep,
    snapshot: &BRep,
    outcome: &SetbackTopologyOutcome,
    info: impl Fn(popkern_core::topology::FaceHandle, usize) -> CornerPatchInfo,
) -> Result<InsertVertexBlendOutput, String> {
    if let Err(e) = brep.validate() {
        *brep = snapshot.clone();
        *SETBACK_TOPOLOGY_LAST
            .lock()
            .unwrap_or_else(|p| p.into_inner()) =
            Some(format!("setback commit failed validation: {e}"));
        return Err(format!("setback vertex blend failed validation: {e}"));
    }
    let infos = outcome
        .cap_faces
        .iter()
        .enumerate()
        .map(|(s, &face)| info(face, s))
        .collect();
    Ok((infos, None, Vec::new()))
}

fn record_decline(vertex: VertexHandle, reason: &str) {
    *SETBACK_TOPOLOGY_LAST
        .lock()
        .unwrap_or_else(|p| p.into_inner()) = Some(format!("setback commit declined: {reason}"));
    if popkern_core::trace::trace_enabled("corner-carve") {
        eprintln!("  [setback] declined at {vertex:?}: {reason}");
    }
}
