//! The setback rung's variational-fit step (PK-574): the fit modes a
//! request admits, run in order (the two-mode ladder concurrently), each
//! attempt recorded. Split from `blend_insert_setback_rung` for the file
//! budget; it shares that module's private rung state.

use popkern_core::topology::{BRep, FaceHandle};

use super::{RungStep, SetbackRung};
use crate::construct::setback_cap_fairness::{CAP_FAIRNESS_BOUND, cap_interior};
use crate::construct::setback_decline::{SetbackCapStrategy, SetbackDecline, SetbackStage};
use crate::construct::setback_vertex_topology::SetbackTopologyOutcome;
use crate::dispatch_overrides;
use crate::types::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerMethod, CornerSubPatchFamily,
    CornerSubPatchProvenance, CornerSubPatchRole, SetbackFitAttempt, SetbackFitDiagnostic,
    SetbackFitMode, SetbackFitOutcome,
};

/// One setback surgery's outcome.
type SurgeryResult = Result<SetbackTopologyOutcome, SetbackDecline>;

impl SetbackRung<'_> {
    /// The G2 ladder's first two attempts run concurrently (PK-598): the
    /// first on `brep`, the second on a clone of the rung's snapshot,
    /// both under the calling thread's test overrides. Both are pure
    /// functions of the snapshot, so each outcome is the one the serial
    /// ladder computes; the caller still consumes them in order, and a
    /// served first attempt discards the second. The second attempt's
    /// wall time hides behind the first (octagon8 mixed: 0.48 s
    /// curvature-weighted then 1.25 s tangent). The corner-cell tangent
    /// rescue runs serially, only after both fail; a G1 ladder runs
    /// serially (its rescue is rarely needed).
    pub(super) fn concurrent_attempts(
        &self,
        brep: &mut BRep,
        modes: &[SetbackFitMode],
    ) -> (Option<SurgeryResult>, Option<(BRep, SurgeryResult)>) {
        #[cfg(feature = "parallel")]
        if let [a @ SetbackFitMode::CurvatureWeighted, b, ..] = *modes {
            let overrides = crate::dispatch_overrides::ThreadOverrides::capture();
            let mut other = self.commit_snapshot.clone();
            let (first, second) = rayon::join(
                // Called from outside the pool, `join` runs BOTH closures
                // on pool workers: each carries the overrides.
                || overrides.apply(|| self.surgery(brep, self.fit_strategy(a))),
                || overrides.apply(|| self.surgery(&mut other, self.fit_strategy(b))),
            );
            return (Some(first), Some((other, second)));
        }
        let _ = (brep, modes);
        (None, None)
    }

    /// The surgery strategy for a fit in `mode`.
    fn fit_strategy(&self, mode: SetbackFitMode) -> SetbackCapStrategy {
        SetbackCapStrategy::VariationalFit {
            mode,
            r_min: self.r_min,
        }
    }

    /// The fit modes the request admits: curvature-weighted under a G2
    /// request, then tangent, then tangent with corner cells.
    pub(super) fn fit_modes(&self) -> &'static [SetbackFitMode] {
        if self.requested_continuity == ContinuityLevel::G2 {
            &[
                SetbackFitMode::CurvatureWeighted,
                SetbackFitMode::Tangent,
                SetbackFitMode::TangentCornerCells,
            ]
        } else {
            &[SetbackFitMode::Tangent, SetbackFitMode::TangentCornerCells]
        }
    }

    /// The variational fit on the carve (PK-574) in `modes` order.
    /// Every attempt is recorded with its interior reading.
    ///
    /// A served fit whose interior bends past [`CAP_FAIRNESS_BOUND`]
    /// does not end the ladder when a later mode is admitted (PK-614):
    /// the later fit runs and the fairer of the two serves, the other
    /// recorded as rejected. Measured on the 6x-disparity cube corner
    /// (r 0.05 / 0.2 / 0.3, G2 request): the curvature-weighted fit
    /// certifies G2 but bends to 30 / r_min beside the thin strip, the
    /// tangent fit to 5.3. The rule is the fan's (PK-612): past the
    /// bound, fairness outranks the claimed level. The corner-cell
    /// rescue (PK-615) serves only under the bound or when fairer than
    /// a held fit.
    pub(super) fn fit(
        &self,
        brep: &mut BRep,
        record: &mut SetbackFitDiagnostic,
        modes: &[SetbackFitMode],
    ) -> RungStep {
        let mut reasons = Vec::new();
        let (mut first, mut ahead) = self.concurrent_attempts(brep, modes);
        // An unfair served fit held while a later mode runs.
        let mut held: Option<HeldFit> = None;
        for &mode in modes {
            let attempt = if let Some(result) = first.take() {
                result
            } else if let Some((attempted, result)) = ahead.take() {
                *brep = attempted;
                result
            } else {
                if held.is_some() {
                    *brep = self.commit_snapshot.clone();
                }
                self.surgery(brep, self.fit_strategy(mode))
            };
            match attempt {
                Ok(outcome) => {
                    let Some(serve) = outcome.fit.clone() else {
                        *brep = self.commit_snapshot.clone();
                        let reason = "fit surgery returned no fit record".to_string();
                        reasons.push(format!("{mode:?} fit: {reason}"));
                        record.attempts.push(SetbackFitAttempt {
                            mode,
                            outcome: SetbackFitOutcome::Declined(reason),
                            interior: None,
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
                            interior: None,
                        });
                        continue;
                    }
                    let interior = if dispatch_overrides::cap_fairness_ladder_disabled() {
                        0.0
                    } else {
                        cap_interior(brep, &outcome.cap_faces, self.r_min).overshoot
                    };
                    record.attempts.push(SetbackFitAttempt {
                        mode,
                        outcome: SetbackFitOutcome::Served(serve.certified),
                        interior: Some(interior),
                    });
                    let candidate = HeldFit {
                        brep: None,
                        attempt: record.attempts.len() - 1,
                        mode,
                        cap_faces: outcome.cap_faces,
                        interior,
                    };
                    let later = modes.last() != Some(&mode);
                    let chosen = match held.take() {
                        // The corner-cell rescue serves only a fair cap:
                        // past the bound it has moved the corner's
                        // obstruction into the interior (pyramid-4 h3
                        // variable: 726 / r_min), not resolved it.
                        None if mode == SetbackFitMode::TangentCornerCells
                            && interior > CAP_FAIRNESS_BOUND =>
                        {
                            *brep = self.commit_snapshot.clone();
                            let reason = format!(
                                "corner-cell rescue interior curvature {interior:.2} / r_min \
                                 over the fairness bound {CAP_FAIRNESS_BOUND}"
                            );
                            reasons.push(format!("{mode:?} fit: {reason}"));
                            if let Some(last) = record.attempts.last_mut() {
                                last.outcome = SetbackFitOutcome::Rejected {
                                    certified: serve.certified,
                                    reason,
                                };
                            }
                            continue;
                        }
                        None if interior > CAP_FAIRNESS_BOUND && later => {
                            held = Some(HeldFit {
                                brep: Some(std::mem::replace(brep, self.commit_snapshot.clone())),
                                ..candidate
                            });
                            continue;
                        }
                        None => candidate,
                        Some(earlier) => {
                            let (winner, loser) = if candidate.interior < earlier.interior {
                                (candidate, earlier)
                            } else {
                                let mut earlier = earlier;
                                if let Some(model) = earlier.brep.take() {
                                    *brep = model;
                                }
                                (earlier, candidate)
                            };
                            reject_fairer(record, &loser, &winner);
                            winner
                        }
                    };
                    return self.serve_fit(record, chosen);
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
                        interior: None,
                    });
                }
            }
        }
        // Every later mode declined: the held (unfair) fit serves.
        if let Some(mut earlier) = held {
            if let Some(model) = earlier.brep.take() {
                *brep = model;
            }
            return self.serve_fit(record, earlier);
        }
        RungStep::Declined(reasons.join("; "))
    }

    /// The served fit's corner records.
    fn serve_fit(&self, record: &mut SetbackFitDiagnostic, fit: HeldFit) -> RungStep {
        let achieved = fit.mode.level();
        record.served = Some(achieved);
        let n = fit.cap_faces.len();
        let infos = fit
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
        RungStep::Decided(Ok((infos, None, diagnostics)))
    }
}

/// A served fit candidate: its model (`Some` only while held aside),
/// its attempt index in the record, mode, cap faces and interior
/// reading.
struct HeldFit {
    brep: Option<BRep>,
    attempt: usize,
    mode: SetbackFitMode,
    cap_faces: Vec<FaceHandle>,
    interior: f64,
}

/// Records `loser` as rejected for the fairer `winner`.
fn reject_fairer(record: &mut SetbackFitDiagnostic, loser: &HeldFit, winner: &HeldFit) {
    if let Some(attempt) = record.attempts.get_mut(loser.attempt)
        && let SetbackFitOutcome::Served(certified) = attempt.outcome
    {
        attempt.outcome = SetbackFitOutcome::Rejected {
            certified,
            reason: format!(
                "interior curvature {:.2} / r_min; the {:?} fit's {:.2} is fairer \
                 (bound {CAP_FAIRNESS_BOUND})",
                loser.interior, winner.mode, winner.interior
            ),
        };
    }
}
