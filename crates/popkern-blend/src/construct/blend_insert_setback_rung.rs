//! The setback vertex-blend rung of `insert_vertex_blend` (PK-094
//! surgery entry, PK-225 N=3 production gate, PK-574 fallback ladder).
//! See the call site for the rung ordering contract and
//! `docs/architecture/variational-corner-fit.md` for the ladder.

use super::InsertVertexBlendOutput;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

use crate::construct::setback_decline::{SetbackCapStrategy, SetbackDecline, SetbackStage};
use crate::construct::setback_vertex_topology::{
    SETBACK_TOPOLOGY_LAST, SetbackTopologyOutcome, setback_fan_commit::fan_seam_g2_failure,
    try_setback_boundary_surgery,
};
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerCapPreference, CornerConstruction,
    CornerMethod, CornerPatchInfo, CornerSubPatchFamily, CornerSubPatchProvenance,
    CornerSubPatchRole, SetbackFitAttempt, SetbackFitDiagnostic, SetbackFitMode, SetbackFitOutcome,
};
use crate::vertex_classifier::VertexBlendClassification;

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
    };
    let mut declines: Vec<String> = Vec::new();
    let order: [bool; 2] = if fitted_first {
        [true, false]
    } else {
        [false, true]
    };
    for fit in order {
        let step = match (fit, try_fit, try_kp20) {
            (true, true, _) => rung.fit(brep, &mut record, rung.fit_modes()),
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
                if record.served.is_none() && !record.attempts.is_empty() {
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
    let mut diagnostics = Vec::new();
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
        if self.requested_continuity != ContinuityLevel::G2 {
            // A KP20 (G2-interior) cap on G1 sweeps yields a G1 blend.
            return serve(brep, ContinuityLevel::G1);
        }
        let Some(reason) = fan_seam_g2_failure(brep, &outcome.cap_faces) else {
            return serve(brep, ContinuityLevel::G2);
        };
        record.kp20_decline = format!("G2 not reached, G1 downgrade held: {reason}");
        let fit_ran = record
            .attempts
            .iter()
            .any(|a| a.mode == SetbackFitMode::CurvatureWeighted);
        if try_fit && !fit_ran {
            let fan = std::mem::replace(brep, self.commit_snapshot.clone());
            if let step @ RungStep::Decided(Ok(_)) =
                self.fit(brep, record, &[SetbackFitMode::CurvatureWeighted])
            {
                return step;
            }
            *brep = fan;
        }
        serve(brep, ContinuityLevel::G1)
    }

    /// The fit modes the request admits: curvature-weighted under a G2
    /// request, then tangent.
    fn fit_modes(&self) -> &'static [SetbackFitMode] {
        if self.requested_continuity == ContinuityLevel::G2 {
            &[SetbackFitMode::CurvatureWeighted, SetbackFitMode::Tangent]
        } else {
            &[SetbackFitMode::Tangent]
        }
    }

    /// The variational fit on the carve (PK-574) in `modes` order.
    /// Every attempt is recorded.
    fn fit(
        &self,
        brep: &mut BRep,
        record: &mut SetbackFitDiagnostic,
        modes: &[SetbackFitMode],
    ) -> RungStep {
        let mut reasons = Vec::new();
        for &mode in modes {
            match self.surgery(brep, SetbackCapStrategy::VariationalFit(mode)) {
                Ok(outcome) => {
                    let Some(serve) = outcome.fit.clone() else {
                        *brep = self.commit_snapshot.clone();
                        let reason = "fit surgery returned no fit record".to_string();
                        reasons.push(format!("{mode:?} fit: {reason}"));
                        record.attempts.push(SetbackFitAttempt {
                            mode,
                            outcome: SetbackFitOutcome::Declined(reason),
                        });
                        continue;
                    };
                    if let Err(e) = brep.validate() {
                        *brep = self.commit_snapshot.clone();
                        let reason = format!("validation: {e}");
                        reasons.push(format!("{mode:?} fit: {reason}"));
                        record.attempts.push(SetbackFitAttempt {
                            mode,
                            outcome: SetbackFitOutcome::Rejected {
                                certified: serve.certified,
                                reason,
                            },
                        });
                        continue;
                    }
                    let achieved = match mode {
                        SetbackFitMode::CurvatureWeighted => ContinuityLevel::G2,
                        SetbackFitMode::Tangent => ContinuityLevel::G1,
                    };
                    record.attempts.push(SetbackFitAttempt {
                        mode,
                        outcome: SetbackFitOutcome::Served(serve.certified),
                    });
                    record.served = Some(achieved);
                    let n = outcome.cap_faces.len();
                    let infos = outcome
                        .cap_faces
                        .iter()
                        .enumerate()
                        .map(|(s, &f)| {
                            self.info(
                                f,
                                CornerMethod::PlateMin,
                                achieved,
                                Some(CornerSubPatchProvenance {
                                    family: CornerSubPatchFamily::VariationalFit,
                                    sector_index: s,
                                    sector_count: n,
                                    piece_index: 0,
                                    piece_count: 1,
                                    role: CornerSubPatchRole::Sector,
                                }),
                            )
                        })
                        .collect();
                    let diagnostics = vec![ApexConstructionDiagnostic::SetbackFit(record.clone())];
                    return RungStep::Decided(Ok((infos, None, diagnostics)));
                }
                Err(d) if d.stage == SetbackStage::Carve => {
                    return RungStep::CarveDeclined(d.reason);
                }
                Err(d) => {
                    if popkern_core::trace::trace_enabled("setback-fit") {
                        eprintln!(
                            "  [setback-fit] {mode:?} declined at {:?}: {}",
                            self.vertex, d.reason
                        );
                    }
                    reasons.push(format!("{mode:?} fit: {}", d.reason));
                    record.attempts.push(SetbackFitAttempt {
                        mode,
                        outcome: match d.fit_certified {
                            Some(certified) => SetbackFitOutcome::Rejected {
                                certified,
                                reason: d.reason,
                            },
                            None => SetbackFitOutcome::Declined(d.reason),
                        },
                    });
                }
            }
        }
        RungStep::Declined(reasons.join("; "))
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
