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
    ApexConstructionDiagnostic, BridgeBow, ContinuityLevel, CornerMethod, CornerSubPatchFamily,
    CornerSubPatchProvenance, CornerSubPatchRole, SetbackFitAttempt, SetbackFitBorder,
    SetbackFitDiagnostic, SetbackFitMode, SetbackFitOutcome,
};

/// TEMPORARY CONTAINMENT (PK-619, replaced by PK-620's solver): the
/// bowed-bridge border runs only at corners of at most this many strips.
/// On an octagon its 16-sided fits push the corner to 400-700 s on the
/// dense solve, and every bowed serve on the 160-row census was at six
/// strips or fewer (the octagon h3 rows it reached served at 173-651 /
/// r_min, a near-pinch). PK-620's sector elimination lifts the limit.
const BOWED_BORDER_MAX_STRIPS: usize = 6;

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
    /// curvature-weighted then 1.25 s tangent). The bridged border's
    /// first pair runs the same way when the ladder reaches it; the
    /// corner-cell rescues run serially, and a G1 ladder runs serially
    /// (its rescue is rarely needed).
    pub(super) fn concurrent_attempts(
        &self,
        brep: &mut BRep,
        ladder: &[FitCandidate],
        rival: Option<f64>,
    ) -> (Option<SurgeryResult>, Option<(BRep, SurgeryResult)>) {
        #[cfg(feature = "parallel")]
        if let [
            a @ FitCandidate {
                mode: SetbackFitMode::CurvatureWeighted,
                ..
            },
            b,
            ..,
        ] = *ladder
        {
            let overrides = crate::dispatch_overrides::ThreadOverrides::capture();
            let mut other = self.commit_snapshot.clone();
            let (first, second) = rayon::join(
                // Called from outside the pool, `join` runs BOTH closures
                // on pool workers: each carries the overrides.
                || overrides.apply(|| self.surgery(brep, self.fit_strategy(a, rival))),
                || overrides.apply(|| self.surgery(&mut other, self.fit_strategy(b, rival))),
            );
            return (Some(first), Some((other, second)));
        }
        let _ = (brep, ladder, rival);
        (None, None)
    }

    /// The surgery strategy for a fit candidate.
    fn fit_strategy(&self, candidate: FitCandidate, rival: Option<f64>) -> SetbackCapStrategy {
        SetbackCapStrategy::VariationalFit {
            mode: candidate.mode,
            r_min: self.r_min,
            border: candidate.border,
            rival,
        }
    }

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
            // weighted fit served on it in the 160-row census, and on an
            // octagon that attempt alone costs up to 62 s (PK-620).
            .filter(|c| {
                !(matches!(c.border, SetbackFitBorder::BowedBridges(_))
                    && (c.mode == SetbackFitMode::CurvatureWeighted
                        || self.n_sides > BOWED_BORDER_MAX_STRIPS))
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
    ) -> RungStep {
        let mut reasons = Vec::new();
        let (mut first, mut ahead) = self.concurrent_attempts(brep, ladder, rival);
        // The fairest served fit so far, held while later attempts run;
        // its model is held aside.
        let mut held: Option<HeldFit> = None;
        for (index, &candidate) in ladder.iter().enumerate() {
            let FitCandidate { mode, .. } = candidate;
            // The bridged pass's first pair runs concurrently too, when
            // neither attempt has run (PK-620: the octagon h1.2 variable
            // apex's bridged tangent fit, 10 s, hides behind its
            // curvature-weighted fit, 62 s).
            if first.is_none()
                && ahead.is_none()
                && candidate.border == SetbackFitBorder::Bridged
                && ladder[..index]
                    .iter()
                    .all(|c| c.border == SetbackFitBorder::Corners)
                && ladder[index..]
                    .iter()
                    .take(2)
                    .all(|c| c.border == SetbackFitBorder::Bridged && !c.ran(record))
            {
                (first, ahead) = self.concurrent_attempts(brep, &ladder[index..], rival);
            }
            let attempt = if let Some(result) = first.take() {
                result
            } else if let Some((attempted, result)) = ahead.take() {
                *brep = attempted;
                result
            } else {
                if candidate.ran(record) {
                    continue;
                }
                *brep = self.commit_snapshot.clone();
                self.surgery(brep, self.fit_strategy(candidate, rival))
            };
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
                        border,
                        outcome: SetbackFitOutcome::Served(serve.certified),
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
                            last.outcome = SetbackFitOutcome::Rejected {
                                certified: serve.certified,
                                reason,
                            };
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
                    first = None;
                    ahead = None;
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
        // Every later attempt declined: the held (unfair) fit serves.
        if let Some(earlier) = held {
            return self.serve_held(brep, record, earlier);
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
