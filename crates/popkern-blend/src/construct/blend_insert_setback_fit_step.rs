//! The setback rung's variational-fit step (PK-574): the fit modes a
//! request admits, run in order (the two-mode ladder concurrently), each
//! attempt recorded. Split from `blend_insert_setback_rung` for the file
//! budget; it shares that module's private rung state.

use popkern_core::topology::BRep;

use super::{RungStep, SetbackRung};
use crate::construct::setback_decline::{SetbackCapStrategy, SetbackDecline, SetbackStage};
use crate::construct::setback_vertex_topology::SetbackTopologyOutcome;
use crate::types::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerMethod, CornerSubPatchFamily,
    CornerSubPatchProvenance, CornerSubPatchRole, SetbackFitAttempt, SetbackFitDiagnostic,
    SetbackFitMode, SetbackFitOutcome,
};

/// One setback surgery's outcome.
type SurgeryResult = Result<SetbackTopologyOutcome, SetbackDecline>;

impl SetbackRung<'_> {
    /// A two-mode ladder's attempts run concurrently (PK-598): the
    /// first on `brep`, the second on a clone of the rung's snapshot,
    /// both under the calling thread's test overrides. Both are pure
    /// functions of the snapshot, so each outcome is the one the serial
    /// ladder computes; the caller still consumes them in order, and a
    /// served first attempt discards the second. The second attempt's
    /// wall time hides behind the first (octagon8 mixed: 0.48 s
    /// curvature-weighted then 1.25 s tangent).
    pub(super) fn concurrent_attempts(
        &self,
        brep: &mut BRep,
        modes: &[SetbackFitMode],
    ) -> (Option<SurgeryResult>, Option<(BRep, SurgeryResult)>) {
        #[cfg(feature = "parallel")]
        if let [a, b] = *modes {
            let overrides = crate::dispatch_overrides::ThreadOverrides::capture();
            let mut other = self.commit_snapshot.clone();
            let (first, second) = rayon::join(
                // Called from outside the pool, `join` runs BOTH closures
                // on pool workers: each carries the overrides.
                || overrides.apply(|| self.surgery(brep, SetbackCapStrategy::VariationalFit(a))),
                || {
                    overrides
                        .apply(|| self.surgery(&mut other, SetbackCapStrategy::VariationalFit(b)))
                },
            );
            return (Some(first), Some((other, second)));
        }
        let _ = (brep, modes);
        (None, None)
    }

    /// The fit modes the request admits: curvature-weighted under a G2
    /// request, then tangent.
    pub(super) fn fit_modes(&self) -> &'static [SetbackFitMode] {
        if self.requested_continuity == ContinuityLevel::G2 {
            &[SetbackFitMode::CurvatureWeighted, SetbackFitMode::Tangent]
        } else {
            &[SetbackFitMode::Tangent]
        }
    }

    /// The variational fit on the carve (PK-574) in `modes` order.
    /// Every attempt is recorded.
    pub(super) fn fit(
        &self,
        brep: &mut BRep,
        record: &mut SetbackFitDiagnostic,
        modes: &[SetbackFitMode],
    ) -> RungStep {
        let mut reasons = Vec::new();
        let (mut first, mut ahead) = self.concurrent_attempts(brep, modes);
        for &mode in modes {
            let attempt = if let Some(result) = first.take() {
                result
            } else if let Some((attempted, result)) = ahead.take() {
                *brep = attempted;
                result
            } else {
                self.surgery(brep, SetbackCapStrategy::VariationalFit(mode))
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
