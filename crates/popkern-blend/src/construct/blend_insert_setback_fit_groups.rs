//! The fit ladder's attempt groups and its deferred proofs (PK-620),
//! split from `blend_insert_setback_fit_step` for the file budget.

use popkern_core::topology::BRep;

use super::super::{RungStep, SetbackRung};
use super::{FitCandidate, HeldFit, SurgeryResult};
use crate::construct::setback_cap_fairness::CAP_FAIRNESS_BOUND;
use crate::construct::setback_decline::{CapCommitError, SetbackCapStrategy};
use crate::construct::setback_vertex_topology::{DeferredProof, prove_fit};
use crate::types::{
    SetbackFitBorder, SetbackFitCertified, SetbackFitDiagnostic, SetbackFitMode, SetbackFitOutcome,
};

/// A fit past the fairness bound whose proof waits for the ladder's
/// end: its held model and record entry, and what certifies it.
pub(super) struct DeferredFit {
    pub fit: HeldFit,
    pub proof: DeferredProof,
    /// A corner-cell or patched rescue no fairer than the rival fan: the
    /// rescue rule (PK-621) rejects it unless a fit was already held
    /// when it ran, that is, unless an earlier fit certifies.
    pub conditional: bool,
}

impl SetbackRung<'_> {
    /// Where the group of attempts starting at `index` ends. A group's
    /// attempts run concurrently, each on a clone of the rung's snapshot
    /// under the calling thread's test overrides; each is a pure function
    /// of the snapshot, so its outcome is the serial ladder's, and the
    /// ladder still consumes them in order (PK-598). The groups:
    ///
    /// - on the corners border, a G2 ladder's leading curvature-weighted
    ///   and tangent pair, which serves most corners (the corner-cell
    ///   rescue runs alone after it: it is wasted work whenever the pair
    ///   is fair);
    /// - every star attempt on the bridged border;
    /// - every tangent attempt on the bowed borders, flat and round;
    /// - every corner-patched attempt.
    ///
    /// The later groups run only at corners whose earlier fits are all
    /// past the fairness bound, where the serial ladder runs each of
    /// their attempts unless one of them is fair (PK-620: the octagon
    /// h1.2 variable apex reached the patched cap at 73 s serially, its
    /// eight bowed and patched attempts at 7-12 s each).
    pub(super) fn group_end(&self, ladder: &[FitCandidate], index: usize) -> usize {
        let first = ladder[index];
        let patched = |c: &FitCandidate| c.mode == SetbackFitMode::TangentCornerPatched;
        let bowed = |c: &FitCandidate| matches!(c.border, SetbackFitBorder::BowedBridges(_));
        let run = |same: &dyn Fn(&FitCandidate) -> bool| {
            index + ladder[index..].iter().take_while(|c| same(c)).count()
        };
        if patched(&first) {
            run(&patched)
        } else if bowed(&first) {
            run(&|c| bowed(c) && !patched(c))
        } else if first.border == SetbackFitBorder::Bridged {
            run(&|c| c.border == SetbackFitBorder::Bridged && !patched(c))
        } else if first.mode == SetbackFitMode::CurvatureWeighted
            && ladder.get(index + 1).is_some_and(|next| {
                next.border == first.border && next.mode == SetbackFitMode::Tangent
            })
        {
            index + 2
        } else {
            index + 1
        }
    }

    /// Run `candidates`, each on a clone of the snapshot: concurrently
    /// when there are several (see [`Self::group_end`]).
    pub(super) fn run_group(
        &self,
        candidates: &[FitCandidate],
        rival: Option<f64>,
        (held, defer): (Option<f64>, bool),
    ) -> Vec<(BRep, SurgeryResult)> {
        let attempt = |&candidate: &FitCandidate| {
            let mut model = self.commit_snapshot.clone();
            let result = self.surgery(
                &mut model,
                self.fit_strategy(candidate, rival, (held, defer)),
            );
            (model, result)
        };
        #[cfg(feature = "parallel")]
        if candidates.len() > 1 {
            use rayon::prelude::*;
            let overrides = crate::dispatch_overrides::ThreadOverrides::capture();
            // Called from outside the pool, the attempts run on pool
            // workers: each carries the overrides.
            return candidates
                .par_iter()
                .map(|c| overrides.apply(|| attempt(c)))
                .collect();
        }
        candidates.iter().map(attempt).collect()
    }

    /// The surgery strategy for a fit candidate.
    fn fit_strategy(
        &self,
        candidate: FitCandidate,
        rival: Option<f64>,
        (held, defer): (Option<f64>, bool),
    ) -> SetbackCapStrategy {
        SetbackCapStrategy::VariationalFit {
            mode: candidate.mode,
            r_min: self.r_min,
            border: candidate.border,
            rival,
            held,
            defer,
        }
    }

    /// The ladder's end with fits deferred (PK-620): certify them
    /// fairest first (the earlier attempt on a tie, as the ladder's own
    /// comparison keeps it) until one serves. That is the fit the eager
    /// ladder serves, the fairest certified one the rescue rule admits,
    /// and only the fits it beats go uncertified: on the octagon h3
    /// variable apex the eager ladder certified 41 s of fits that a
    /// fairer fit then replaced.
    ///
    /// The rescue rule admits a conditional fit only when a fit was held
    /// when it ran. The first fit the eager ladder holds is never
    /// conditional (nothing is held before it), so a conditional fit is
    /// admitted exactly when an earlier unconditional fit certifies;
    /// those are certified in ladder order only as far as needed.
    pub(super) fn conclude_deferred(
        &self,
        brep: &mut BRep,
        record: &mut SetbackFitDiagnostic,
        mut deferred: Vec<DeferredFit>,
        reasons: &mut Vec<String>,
    ) -> RungStep {
        deferred.sort_by_key(|d| d.fit.attempt);
        let mut fairest: Vec<usize> = (0..deferred.len()).collect();
        fairest.sort_by(|&a, &b| {
            let (a, b) = (&deferred[a].fit, &deferred[b].fit);
            a.interior
                .total_cmp(&b.interior)
                .then(a.attempt.cmp(&b.attempt))
        });
        let mut proofs: Vec<Option<Result<SetbackFitCertified, CapCommitError>>> =
            deferred.iter().map(|_| None).collect();
        let mut winner = None;
        for &k in &fairest {
            if self.proof(&deferred, &mut proofs, k).is_err() {
                continue;
            }
            if deferred[k].conditional
                && !(0..k).any(|e| {
                    !deferred[e].conditional && self.proof(&deferred, &mut proofs, e).is_ok()
                })
            {
                continue;
            }
            winner = Some(k);
            break;
        }
        // The record: each fit's proof outcome, or why it went unproved.
        for (k, d) in deferred.iter().enumerate() {
            let Some(attempt) = record.attempts.get_mut(d.fit.attempt) else {
                continue;
            };
            let (mode, border) = (d.fit.mode, d.fit.border);
            attempt.outcome = match (&proofs[k], winner) {
                (Some(Ok(certified)), Some(w)) if w == k => SetbackFitOutcome::Served(*certified),
                (Some(Ok(certified)), Some(w)) => SetbackFitOutcome::Rejected {
                    certified: *certified,
                    reason: fairer_reason(&d.fit, &deferred[w].fit),
                },
                (Some(Ok(certified)), None) => {
                    let reason = format!(
                        "corner-cell rescue interior curvature {:.2} / r_min over the fairness \
                         bound {CAP_FAIRNESS_BOUND}, no fairer than the rival fan, with no fit \
                         held before it",
                        d.fit.interior
                    );
                    reasons.push(format!("{mode:?} fit ({border:?}): {reason}"));
                    SetbackFitOutcome::Rejected {
                        certified: *certified,
                        reason,
                    }
                }
                (Some(Err(e)), _) => {
                    reasons.push(format!("{mode:?} fit ({border:?}): {}", e.reason));
                    match e.fit_certified {
                        Some(certified) => SetbackFitOutcome::Rejected {
                            certified,
                            reason: e.reason.clone(),
                        },
                        None => SetbackFitOutcome::Declined(e.reason.clone()),
                    }
                }
                (None, Some(w)) => SetbackFitOutcome::Declined(format!(
                    "not certified: {}",
                    fairer_reason(&d.fit, &deferred[w].fit)
                )),
                (None, None) => SetbackFitOutcome::Declined(
                    "not certified: the rescue rule rejects it (no fit held before it)".to_string(),
                ),
            };
        }
        match winner {
            Some(w) => {
                let mut fit = deferred.swap_remove(w).fit;
                if let Some(model) = fit.brep.take() {
                    *brep = model;
                }
                self.serve_fit(record, fit)
            }
            None => {
                *brep = self.commit_snapshot.clone();
                RungStep::Declined(reasons.join("; "))
            }
        }
    }

    /// Deferred fit `k`'s proof, run once on its held model.
    fn proof<'p>(
        &self,
        deferred: &[DeferredFit],
        proofs: &'p mut [Option<Result<SetbackFitCertified, CapCommitError>>],
        k: usize,
    ) -> &'p Result<SetbackFitCertified, CapCommitError> {
        proofs[k].get_or_insert_with(|| {
            let d = &deferred[k];
            match d.fit.brep.as_ref() {
                Some(model) => prove_fit(
                    model,
                    &d.fit.cap_faces,
                    self.fillet_surfaces_pre,
                    d.fit.mode,
                    &d.proof,
                ),
                None => Err(CapCommitError::from("the deferred fit's model is missing")),
            }
        })
    }
}

/// Why `loser` was not served: `winner` is fairer.
fn fairer_reason(loser: &HeldFit, winner: &HeldFit) -> String {
    format!(
        "interior curvature {:.2} / r_min; the {:?} fit's {:.2} on the {:?} border is fairer \
         (bound {CAP_FAIRNESS_BOUND})",
        loser.interior, winner.mode, winner.interior, winner.border
    )
}

/// Record the deferred fits a served fit outranks: never certified.
pub(super) fn outranked(record: &mut SetbackFitDiagnostic, rest: &[DeferredFit], served: &HeldFit) {
    for fit in rest {
        if let Some(attempt) = record.attempts.get_mut(fit.fit.attempt) {
            attempt.outcome = SetbackFitOutcome::Declined(format!(
                "not certified: interior curvature {:.2} / r_min; the {:?} fit's {:.2} on the \
                 {:?} border is fairer and serves",
                fit.fit.interior, served.mode, served.interior, served.border
            ));
        }
    }
}
