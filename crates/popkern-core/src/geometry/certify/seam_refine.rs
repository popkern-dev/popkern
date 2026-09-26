//! Branch-and-bound refinement of one seam under one correspondence
//! mode (PK-499; see [`super::seam::certify_seam`] for the policy).

use super::seam::{CertificateStatus, InconclusiveReason, SeamCertificate, SeamCertifyOptions};
use super::seam_interval::{IntervalBounds, SideB, interval_bounds};
use super::seam_state::SeamSideSource;
use std::collections::HashMap;

/// Enclosures already computed for one seam under fraction
/// correspondence, keyed by the s-interval's exact endpoints. Under
/// fraction correspondence an interval's bounds depend only on its
/// endpoints (the sources cache span nets, nothing else), so a later
/// pass over the same seam — the full-budget run after a short probe
/// ran out — reuses them bit for bit instead of re-evaluating
/// (measured: 27-34% of all evidence certification work on the G2
/// pyramid and boss-plate corpus was such re-evaluation, PK-562).
pub(super) type FractionMemo = HashMap<(u64, u64), Result<IntervalBounds, InconclusiveReason>>;

/// Consecutive non-improving bisection generations after which an
/// interval settles at its bound (see the stall rule in
/// [`super::seam::certify_seam`]).
const STALL_GENERATIONS: u8 = 3;

/// Which failure a certificate reports when several occurred: the
/// most structural one (a side that cannot be enclosed, a degenerate
/// chart, a correspondence failure) outranks running out of
/// resolution or budget, which only says the bounds are loose.
fn severity(reason: InconclusiveReason) -> u8 {
    match reason {
        InconclusiveReason::NotASeam => 5,
        InconclusiveReason::Side(..) => 4,
        InconclusiveReason::Jacobian => 3,
        InconclusiveReason::Correspondence => 2,
        InconclusiveReason::Resolution => 1,
        InconclusiveReason::Budget => 0,
    }
}

/// [`super::seam::certify_seam`] under one correspondence mode.
/// `memo` (fraction correspondence only) reuses and records interval
/// enclosures across passes over the same seam.
pub(super) fn certify_with(
    a: &mut SeamSideSource,
    b: &mut SideB<'_>,
    options: &SeamCertifyOptions,
    mut memo: Option<&mut FractionMemo>,
) -> SeamCertificate {
    debug_assert!(memo.is_none() || matches!(b, SideB::Fraction(_)));
    // Projection starts finer: a Krawczyk step over a quarter of a
    // turning rail never contracts (measured on the boss-plate rings).
    let initial = if matches!(b, SideB::Projected(_)) {
        16
    } else {
        4
    };
    let mut cuts: Vec<f64> = (0..=initial).map(|k| k as f64 / initial as f64).collect();
    cuts.extend(a.fraction_breaks());
    if let SideB::Fraction(b) = b {
        cuts.extend(b.fraction_breaks());
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    // Each queued interval carries its parent's overshoot and the
    // count of consecutive generations that failed to improve on it.
    let mut to_evaluate: Vec<((f64, f64), f64, u8)> = cuts
        .windows(2)
        .map(|w| ((w[0], w[1]), f64::INFINITY, 0))
        .collect();
    let targets = [
        options.targets.gap,
        options.targets.normal_sine,
        options.targets.curvature,
    ];
    let mut lower = [0.0_f64; 3];
    let mut pending: Vec<((f64, f64), IntervalBounds, f64, u8)> = Vec::new();
    let mut settled: Vec<((f64, f64), IntervalBounds)> = Vec::new();
    let mut failure: Option<(InconclusiveReason, (f64, f64))> = None;
    let mut record = |reason: InconclusiveReason, at: (f64, f64)| {
        if failure.is_none_or(|(held, _)| severity(reason) > severity(held)) {
            failure = Some((reason, at));
        }
    };
    let (mut evaluated, mut uncertified) = (0_usize, 0.0_f64);
    let mut structural_stop = false;
    let split = |(lo, hi): (f64, f64)| {
        let m = 0.5 * lo + 0.5 * hi;
        [(lo, m), (m, hi)]
    };
    loop {
        let fraction_mode = matches!(b, SideB::Fraction(_));
        while let Some(((lo, hi), parent_over, parent_stall)) = to_evaluate.pop() {
            // A certified fraction mismatch settles the mode: stop and
            // let the caller switch to projected correspondence.
            if fraction_mode && lower[0] > options.correspondence_tolerance {
                break;
            }
            evaluated += 1;
            let key = (lo.to_bits(), hi.to_bits());
            let result = match memo.as_deref_mut().and_then(|m| m.get(&key)) {
                Some(known) => *known,
                None => {
                    let fresh = interval_bounds(a, b, lo, hi);
                    if let Some(m) = memo.as_deref_mut() {
                        m.insert(key, fresh);
                    }
                    fresh
                }
            };
            match result {
                Ok(bounds) => {
                    for (held, witness) in lower.iter_mut().zip(bounds.lower) {
                        *held = held.max(witness);
                    }
                    pending.push(((lo, hi), bounds, parent_over, parent_stall));
                }
                Err(reason) => {
                    let room = evaluated + to_evaluate.len() + 2 <= options.budget;
                    let floor = if reason == InconclusiveReason::Jacobian {
                        options.singular_min_width
                    } else {
                        options.min_width
                    };
                    if hi - lo > floor && room {
                        to_evaluate.extend(split((lo, hi)).map(|c| (c, f64::INFINITY, 0)));
                        continue;
                    }
                    let reason = match reason {
                        InconclusiveReason::Resolution if hi - lo > options.min_width => {
                            InconclusiveReason::Budget
                        }
                        other => other,
                    };
                    record(reason, (lo, hi));
                    uncertified += hi - lo;
                    if options.stop_at_structural_failure
                        && !matches!(
                            reason,
                            InconclusiveReason::Budget | InconclusiveReason::Resolution
                        )
                    {
                        uncertified += to_evaluate
                            .iter()
                            .map(|(sp, _, _)| sp.1 - sp.0)
                            .sum::<f64>();
                        to_evaluate.clear();
                        structural_stop = true;
                        break;
                    }
                }
            }
        }
        if structural_stop || (fraction_mode && lower[0] > options.correspondence_tolerance) {
            break;
        }
        let overshoot = |bounds: &IntervalBounds| -> f64 {
            (0..3)
                .map(|k| {
                    if let Some(gates) = options.decision_gates
                        && (lower[k] > gates[k] || bounds.upper[k] <= gates[k])
                    {
                        return -1.0; // decided: certified fail or pass
                    }
                    let tol = targets[k] + options.targets.relative * lower[k];
                    (bounds.upper[k] - lower[k] - tol) / tol
                })
                .fold(f64::NEG_INFINITY, f64::max)
        };
        let worst = pending
            .iter()
            .enumerate()
            .map(|(k, (_, bounds, _, _))| (k, overshoot(bounds)))
            .max_by(|x, y| x.1.total_cmp(&y.1));
        let Some((index, over)) = worst.filter(|w| w.1 > 0.0) else {
            break;
        };
        let (span, bounds, parent_over, parent_stall) = pending.swap_remove(index);
        // Stall: a bound that three generations of bisection have not
        // cut by 10% is not converging (a pole or blended corner
        // neighbourhood whose enclosures blow up toward the
        // singularity); it settles at its rigorous bound.
        let stall = if over > 0.9 * parent_over {
            parent_stall + 1
        } else {
            0
        };
        let room = evaluated + 2 <= options.budget;
        if span.1 - span.0 > options.min_width && room && stall < STALL_GENERATIONS {
            to_evaluate.extend(split(span).map(|c| (c, over, stall)));
        } else {
            let reason = if room {
                InconclusiveReason::Resolution
            } else {
                InconclusiveReason::Budget
            };
            record(reason, span);
            settled.push((span, bounds));
        }
    }
    settled.extend(
        pending
            .into_iter()
            .map(|(span, bounds, _, _)| (span, bounds)),
    );
    let mut out = SeamCertificate::inconclusive(InconclusiveReason::Budget);
    out.g0_upper = 0.0;
    out.g1_upper = 0.0;
    out.g2_upper = 0.0;
    out.jacobian_lower = [f64::INFINITY; 2];
    out.lower = lower;
    out.intervals = evaluated;
    out.uncertified_length = uncertified;
    for (span, bounds) in &settled {
        let worst = [&mut out.g0_upper, &mut out.g1_upper, &mut out.g2_upper];
        for (k, slot) in worst.into_iter().enumerate() {
            if bounds.upper[k] > *slot {
                *slot = bounds.upper[k];
                out.worst_at[k] = *span;
            }
        }
        for k in 0..2 {
            out.jacobian_lower[k] = out.jacobian_lower[k].min(bounds.jacobian[k]);
        }
    }
    if settled.is_empty() {
        out.jacobian_lower = [0.0; 2];
    }
    // Correspondence fails only when a gap above tolerance is
    // certified; a loose upper bound is a budget matter.
    if out.lower[0] > options.correspondence_tolerance {
        record(InconclusiveReason::Correspondence, out.worst_at[0]);
    } else if out.g0_upper > options.correspondence_tolerance {
        record(InconclusiveReason::Budget, out.worst_at[0]);
    }
    out.status = match failure {
        None => CertificateStatus::Certified,
        Some((reason, at)) => CertificateStatus::Inconclusive { reason, at },
    };
    out
}
