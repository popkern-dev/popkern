//! The setback rung's variational-fit step (PK-574): the fit modes a
//! request admits, run in order (the two-mode ladder concurrently), each
//! attempt recorded. Split from `blend_insert_setback_rung` for the file
//! budget; it shares that module's private rung state.

use popkern_core::topology::{BRep, FaceHandle};

use super::{RungStep, SetbackRung};
use crate::construct::setback_cap_fairness::{CAP_FAIRNESS_BOUND, cap_interior};
use crate::construct::setback_decline::{SetbackDecline, SetbackStage};
use crate::construct::setback_vertex_topology::{FitProof, SetbackTopologyOutcome};
use crate::dispatch_overrides;
use crate::types::{
    ApexConstructionDiagnostic, BridgeBow, ContinuityLevel, CornerMethod, CornerSubPatchFamily,
    CornerSubPatchProvenance, CornerSubPatchRole, SetbackFitAttempt, SetbackFitBorder,
    SetbackFitDiagnostic, SetbackFitMode, SetbackFitOutcome,
};

/// One setback surgery's outcome.
type SurgeryResult = Result<SetbackTopologyOutcome, SetbackDecline>;

#[path = "blend_insert_setback_fit_groups.rs"]
mod groups;
use groups::{DeferredFit, outranked};

impl SetbackRung<'_> {
    /// The fit modes the request admits: curvature-weighted under a G2
    /// request, then tangent, then tangent with corner cells, then the
    /// corner-patched fit (PK-619; the ladder runs it after every star
    /// attempt).
    pub(super) fn fit_modes(&self) -> &'static [SetbackFitMode] {
        if self.requested_continuity == ContinuityLevel::G2 {
            &[
                SetbackFitMode::CurvatureWeighted,
                SetbackFitMode::Tangent,
                SetbackFitMode::TangentCornerCells,
                SetbackFitMode::TangentCornerPatched,
            ]
        } else {
            &[
                SetbackFitMode::Tangent,
                SetbackFitMode::TangentCornerCells,
                SetbackFitMode::TangentCornerPatched,
            ]
        }
    }

    /// The ladder over `modes`: every star mode on the corners border,
    /// then on the bridged border (PK-618), then the tangent modes on
    /// the flat and then the round bowed-bridge borders (PK-619); then
    /// the corner-patched fit on the bridged and bowed borders (PK-619
    /// stage 3), less the attempts the record already holds. A later
    /// attempt runs only while every earlier fit is past the fairness
    /// bound, so a star cap fair on any border is never traded for the
    /// patched cap's extra seams (measured: pyramid-4 h3 wide G2 served
    /// the patched cap on the bridged border at 3.46 / r_min, before
    /// the flat bow's star at 1.79). The corner-patched fit never served
    /// on the corners border in the 160-row census, so it is not run
    /// there.
    pub(super) fn fit_ladder(
        &self,
        record: &SetbackFitDiagnostic,
        modes: &[SetbackFitMode],
    ) -> Vec<FitCandidate> {
        const BORDERS: [SetbackFitBorder; 4] = [
            SetbackFitBorder::Corners,
            SetbackFitBorder::Bridged,
            SetbackFitBorder::BowedBridges(BridgeBow::Flat),
            SetbackFitBorder::BowedBridges(BridgeBow::Round),
        ];
        let patched = |mode: &SetbackFitMode| *mode == SetbackFitMode::TangentCornerPatched;
        let star = BORDERS.into_iter().flat_map(|border| {
            modes
                .iter()
                .filter(|m| !patched(m))
                .map(move |&mode| FitCandidate { mode, border })
        });
        let corner_patched = BORDERS[1..].iter().flat_map(|&border| {
            modes
                .iter()
                .filter(|m| patched(m))
                .map(move |&mode| FitCandidate { mode, border })
        });
        star.chain(corner_patched)
            // The bowed border runs the tangent fits only: no curvature-
            // weighted fit served on it in the 160-row census.
            .filter(|c| {
                !(matches!(c.border, SetbackFitBorder::BowedBridges(_))
                    && c.mode == SetbackFitMode::CurvatureWeighted)
            })
            .filter(|c| !c.ran(record))
            .collect()
    }

    /// The variational fit on the carve (PK-574), in `ladder` order.
    /// Every attempt is recorded with its border and interior reading.
    ///
    /// A served fit whose interior bends past [`CAP_FAIRNESS_BOUND`]
    /// does not end the ladder while later attempts remain (PK-614): it
    /// is held, the later attempts run, and the fairest serves, each
    /// fit it beat recorded as rejected. Measured on the 6x-disparity
    /// cube corner (r 0.05 / 0.2 / 0.3, G2 request): the
    /// curvature-weighted fit certifies G2 but bends to 30 / r_min
    /// beside the thin strip, the tangent fit to 5.3. The rule is the
    /// fan's (PK-612): past the bound, fairness outranks the claimed
    /// level. The corner-cell rescue (PK-615) serves only under the
    /// bound, when fairer than `rival` (the reading of an unfair fan the
    /// rung holds), or when there is no rival (PK-621).
    ///
    /// The bridged border's attempts follow the corners border's
    /// (PK-618), so a corner whose corners-border fits are fair never
    /// pays for them. Measured 2026-10-01 on the pentagon h1.2 wide
    /// apex (G2): the corners border's fits read 41 and 277 (the fan
    /// 8.53), the bridged curvature-weighted fit 3.89. A corners-border
    /// attempt that overran a narrow strip already carved the bridged
    /// border (PK-611), and its mode does not rerun there.
    pub(super) fn fit(
        &self,
        brep: &mut BRep,
        record: &mut SetbackFitDiagnostic,
        ladder: &[FitCandidate],
        rival: Option<f64>,
        fan_bar: Option<f64>,
    ) -> RungStep {
        let mut reasons = Vec::new();
        // Each attempt's model and outcome, filled a group at a time.
        let mut pending: Vec<Option<(BRep, SurgeryResult)>> = ladder.iter().map(|_| None).collect();
        // The fairest served fit so far, held while later attempts run;
        // its model is held aside.
        let mut held: Option<HeldFit> = None;
        // Fits past the fairness bound defer their proof to the ladder's
        // end (PK-620, see `conclude_deferred`).
        let defer = !dispatch_overrides::cap_fairness_ladder_disabled();
        let mut deferred: Vec<DeferredFit> = Vec::new();
        for (index, &candidate) in ladder.iter().enumerate() {
            let FitCandidate { mode, .. } = candidate;
            if pending[index].is_none() {
                if candidate.ran(record) {
                    continue;
                }
                // A group's later attempts may already be on the record
                // (a corners attempt that overran ran on the bridged
                // border, PK-611).
                let members: Vec<usize> = (index..self.group_end(ladder, index))
                    .filter(|&i| i == index || !ladder[i].ran(record))
                    .collect();
                let candidates: Vec<FitCandidate> = members.iter().map(|&i| ladder[i]).collect();
                // The reading a fit must beat to serve: the held fit's,
                // and an unfair fan's (PK-620, see `fairness_cut`).
                let bar = match (held.as_ref().map(|h| h.interior), fan_bar) {
                    (Some(h), Some(f)) => Some(h.min(f)),
                    (h, f) => h.or(f),
                };
                let results = self.run_group(&candidates, rival, (bar, defer));
                for (i, result) in members.into_iter().zip(results) {
                    pending[i] = Some(result);
                }
            }
            let Some((model, attempt)) = pending[index].take() else {
                continue;
            };
            *brep = model;
            match attempt {
                Ok(outcome) => {
                    let border = outcome.border;
                    let Some(serve) = outcome.fit.clone() else {
                        *brep = self.commit_snapshot.clone();
                        let reason = "fit surgery returned no fit record".to_string();
                        reasons.push(format!("{mode:?} fit ({border:?}): {reason}"));
                        record.attempts.push(SetbackFitAttempt {
                            mode,
                            border,
                            outcome: SetbackFitOutcome::Declined(reason),
                            interior: None,
                        });
                        continue;
                    };
                    if let Err(e) = brep.validate() {
                        *brep = self.commit_snapshot.clone();
                        let reason = format!("validation: {e}");
                        reasons.push(format!("{mode:?} fit ({border:?}): {reason}"));
                        record.attempts.push(SetbackFitAttempt {
                            mode,
                            border,
                            outcome: match serve.proof {
                                FitProof::Certified(certified) => {
                                    SetbackFitOutcome::Rejected { certified, reason }
                                }
                                FitProof::Deferred(_) => SetbackFitOutcome::Declined(reason),
                            },
                            interior: None,
                        });
                        continue;
                    }
                    let interior = if dispatch_overrides::cap_fairness_ladder_disabled() {
                        0.0
                    } else if let Some(interior) = serve.interior {
                        interior
                    } else {
                        cap_interior(brep, &outcome.cap_faces, self.r_min).overshoot
                    };
                    let certified = match serve.proof {
                        FitProof::Certified(certified) => certified,
                        FitProof::Deferred(proof) => {
                            record.attempts.push(SetbackFitAttempt {
                                mode,
                                border,
                                outcome: SetbackFitOutcome::Declined(
                                    "proof deferred to the ladder's end".to_string(),
                                ),
                                interior: Some(interior),
                            });
                            deferred.push(DeferredFit {
                                // The rescue rule's reach (see below).
                                conditional: matches!(
                                    mode,
                                    SetbackFitMode::TangentCornerCells
                                        | SetbackFitMode::TangentCornerPatched
                                ) && rival.is_some_and(|fan| interior >= fan),
                                fit: HeldFit {
                                    brep: Some(std::mem::replace(
                                        brep,
                                        self.commit_snapshot.clone(),
                                    )),
                                    attempt: record.attempts.len() - 1,
                                    mode,
                                    border,
                                    cap_faces: outcome.cap_faces,
                                    pieces: serve.pieces.clone(),
                                    interior,
                                },
                                proof,
                            });
                            continue;
                        }
                    };
                    record.attempts.push(SetbackFitAttempt {
                        mode,
                        border,
                        outcome: SetbackFitOutcome::Served(certified),
                        interior: Some(interior),
                    });
                    let candidate = HeldFit {
                        brep: None,
                        attempt: record.attempts.len() - 1,
                        mode,
                        border,
                        cap_faces: outcome.cap_faces,
                        pieces: serve.pieces.clone(),
                        interior,
                    };
                    // Past the bound the corner-cell rescue yields to a
                    // fairer rival fan: it can move the corner's
                    // obstruction into the interior (pyramid-4 h3
                    // variable: 726 / r_min against the fan's 27). With
                    // no rival it is held like any other fit and the
                    // fairest serves: declining sends the corner to a
                    // vertex-blend cap that does not certify G1
                    // (pyramid-3 h1.2 wide: a G0 Gregory cap at 418
                    // against the rescue's 10.4; pyramid-5 and -6 h0.4
                    // wide: no cap at all) (PK-621).
                    if held.is_none()
                        && matches!(
                            mode,
                            SetbackFitMode::TangentCornerCells
                                | SetbackFitMode::TangentCornerPatched
                        )
                        && interior > CAP_FAIRNESS_BOUND
                        && rival.is_some_and(|fan| interior >= fan)
                    {
                        *brep = self.commit_snapshot.clone();
                        let reason = format!(
                            "corner-cell rescue interior curvature {interior:.2} / r_min \
                             over the fairness bound {CAP_FAIRNESS_BOUND}, no fairer than \
                             the fan's {:.2}",
                            rival.unwrap_or(f64::INFINITY)
                        );
                        reasons.push(format!("{mode:?} fit ({border:?}): {reason}"));
                        if let Some(last) = record.attempts.last_mut() {
                            last.outcome = SetbackFitOutcome::Rejected { certified, reason };
                        }
                        continue;
                    }
                    // The fairer of the held fit and this one; its model
                    // is left live in `brep`.
                    let best = match held.take() {
                        None => candidate,
                        Some(mut earlier) => {
                            let (winner, loser) = if candidate.interior < earlier.interior {
                                (candidate, earlier)
                            } else {
                                if let Some(model) = earlier.brep.take() {
                                    *brep = model;
                                }
                                (earlier, candidate)
                            };
                            reject_fairer(record, &loser, &winner);
                            winner
                        }
                    };
                    let later = index + 1 < ladder.len();
                    if best.interior > CAP_FAIRNESS_BOUND && later {
                        held = Some(HeldFit {
                            brep: Some(std::mem::replace(brep, self.commit_snapshot.clone())),
                            ..best
                        });
                        continue;
                    }
                    outranked(record, &deferred, &best);
                    return self.serve_fit(record, best);
                }
                Err(d)
                    if d.stage == SetbackStage::Carve && d.border == SetbackFitBorder::Corners =>
                {
                    // No construction on the corners carve can serve; the
                    // bridged border is a different carve.
                    if let Some(earlier) = held.take() {
                        return self.serve_held(brep, record, earlier);
                    }
                    if !deferred.is_empty() {
                        return self.conclude_deferred(brep, record, deferred, &mut reasons);
                    }
                    if !ladder[index..]
                        .iter()
                        .any(|c| c.border != SetbackFitBorder::Corners)
                    {
                        return RungStep::CarveDeclined(d.reason);
                    }
                    reasons.push(format!("{mode:?} fit (Corners): {}", d.reason));
                    record.attempts.push(SetbackFitAttempt {
                        mode,
                        border: d.border,
                        outcome: SetbackFitOutcome::Declined(d.reason),
                        interior: None,
                    });
                    // Every corners-border attempt shares the carve.
                    for (slot, c) in pending.iter_mut().zip(ladder).skip(index + 1) {
                        if c.border == SetbackFitBorder::Corners {
                            *slot = None;
                        }
                    }
                    skip_corners(record, ladder, index);
                }
                Err(d) => {
                    if popkern_core::trace::trace_enabled("setback-fit") {
                        eprintln!(
                            "  [setback-fit] {mode:?} ({:?}) declined at {:?}: {}",
                            d.border, self.vertex, d.reason
                        );
                    }
                    reasons.push(format!("{mode:?} fit ({:?}): {}", d.border, d.reason));
                    record.attempts.push(SetbackFitAttempt {
                        mode,
                        border: d.border,
                        outcome: match d.fit_certified {
                            Some(certified) => SetbackFitOutcome::Rejected {
                                certified: *certified,
                                reason: d.reason,
                            },
                            None => SetbackFitOutcome::Declined(d.reason),
                        },
                        interior: d.interior,
                    });
                }
            }
        }
        // Every later attempt declined: the held (unfair) fit serves.
        if let Some(earlier) = held {
            return self.serve_held(brep, record, earlier);
        }
        if !deferred.is_empty() {
            return self.conclude_deferred(brep, record, deferred, &mut reasons);
        }
        RungStep::Declined(reasons.join("; "))
    }

    /// Serve a held fit, restoring its model.
    fn serve_held(
        &self,
        brep: &mut BRep,
        record: &mut SetbackFitDiagnostic,
        mut fit: HeldFit,
    ) -> RungStep {
        if let Some(model) = fit.brep.take() {
            *brep = model;
        }
        self.serve_fit(record, fit)
    }

    /// The served fit's corner records.
    fn serve_fit(&self, record: &mut SetbackFitDiagnostic, fit: HeldFit) -> RungStep {
        let achieved = fit.mode.level();
        record.served = Some(achieved);
        // A corner-patched cap's faces name their piece (0: the main cap)
        // and their sector within it (PK-619).
        let piece_count = fit.pieces.iter().map(|p| p.piece + 1).max().unwrap_or(1);
        let infos = fit
            .cap_faces
            .iter()
            .zip(&fit.pieces)
            .map(|(&f, place)| {
                self.info(
                    f,
                    CornerMethod::PlateMin,
                    achieved,
                    Some(CornerSubPatchProvenance {
                        family: CornerSubPatchFamily::VariationalFit,
                        sector_index: place.sector,
                        sector_count: place.sectors,
                        piece_index: place.piece,
                        piece_count,
                        role: if place.piece == 0 {
                            CornerSubPatchRole::Sector
                        } else {
                            CornerSubPatchRole::CornerPatchSector
                        },
                    }),
                )
            })
            .collect();
        let diagnostics = vec![ApexConstructionDiagnostic::SetbackFit(record.clone())];
        RungStep::Decided(Ok((infos, None, diagnostics)))
    }
}

/// A served fit: its model (`Some` only while held aside), its attempt
/// index in the record, mode, border, cap faces (with each face's piece
/// and sector) and interior reading.
struct HeldFit {
    brep: Option<BRep>,
    attempt: usize,
    mode: SetbackFitMode,
    border: SetbackFitBorder,
    cap_faces: Vec<FaceHandle>,
    pieces: Vec<crate::construct::setback_vertex_topology::FacePiece>,
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
                "interior curvature {:.2} / r_min; the {:?} fit's {:.2} on the {:?} \
                 border is fairer (bound {CAP_FAIRNESS_BOUND})",
                loser.interior, winner.mode, winner.interior, winner.border
            ),
        };
    }
}

/// One attempt of the fit ladder: a mode on a border.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct FitCandidate {
    pub mode: SetbackFitMode,
    pub border: SetbackFitBorder,
}

impl FitCandidate {
    /// Whether the record already holds this attempt (a corners-border
    /// attempt that overran is recorded on the bridged border).
    fn ran(self, record: &SetbackFitDiagnostic) -> bool {
        record
            .attempts
            .iter()
            .any(|a| a.mode == self.mode && a.border == self.border)
    }
}

/// Record the corners-border attempts after `index` as declined with
/// the shared carve: they cannot carve where it could not.
fn skip_corners(record: &mut SetbackFitDiagnostic, ladder: &[FitCandidate], index: usize) {
    for c in &ladder[index + 1..] {
        if c.border == SetbackFitBorder::Corners && !c.ran(record) {
            record.attempts.push(SetbackFitAttempt {
                mode: c.mode,
                border: c.border,
                outcome: SetbackFitOutcome::Declined(
                    "not attempted: the corners border's carve declined".to_string(),
                ),
                interior: None,
            });
        }
    }
}
