//! One s-interval's certified bounds (PK-499): the per-branch
//! Taylor-form enclosures of every seam quantity, united, with the
//! true midpoint value as the certified witness. See [`super::seam`].

use super::CertifyDecline;
use super::interval::Iv;
use super::seam::InconclusiveReason;
use super::seam_projection::{ProjectedSide, ProjectionFailure};
use super::seam_state::{SeamQuantities, SeamSideSource, Series, seam_quantities};

/// Side b under one of the two correspondence modes.
pub(super) enum SideB<'x> {
    /// Compared at equal seam fractions.
    Fraction(&'x mut SeamSideSource),
    /// Compared at certified foot points.
    Projected(&'x mut ProjectedSide),
    /// Compared at the fixed phase-shifted fraction
    /// `σ = phase + dir·s (mod 1)` of a closed side b.
    Shifted {
        /// Side b.
        side: &'x mut ProjectedSide,
        /// `σ(0)`.
        phase: f64,
        /// `dσ/ds`, ±1.
        dir: f64,
    },
}

/// Most smooth branch pairs one interval is enclosed through before
/// falling back to the hulled natural enclosure.
const MAX_BRANCH_PAIRS: usize = 16;

/// Bounds of one interval for (gap, sine, curvature): upper bounds
/// over the interval and certified lower witnesses at its midpoint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct IntervalBounds {
    pub(super) upper: [f64; 3],
    pub(super) lower: [f64; 3],
    pub(super) jacobian: [f64; 2],
}

pub(super) fn interval_bounds(
    a: &mut SeamSideSource,
    b: &mut SideB<'_>,
    lo: f64,
    hi: f64,
) -> Result<IntervalBounds, InconclusiveReason> {
    let m = 0.5 * lo + 0.5 * hi;
    let (whole, point) = (Iv::new(lo, hi), Iv::point(m));
    let side = |k: usize| move |e: CertifyDecline| InconclusiveReason::Side(k, e);
    // The true midpoint state (hulled at a knot) supplies the witness
    // and the hemisphere of side b's normal.
    let (a_m, _) = a.state(point).map_err(side(0))?;
    // Each piece: (midpoint quantities, interval quantities, whether
    // only the natural enclosure is sound).
    let mut pieces: Vec<(SeamQuantities, SeamQuantities, bool)> = Vec::new();
    // Certified true values at the midpoint: the witness is the least
    // of them (none certified: no witness, lower bound 0).
    let witness_quantities: Vec<SeamQuantities> = match b {
        SideB::Fraction(src) => {
            let (b_m, _) = src.state(point).map_err(side(1))?;
            let flip = seam_quantities(&a_m, &b_m, false).normal_dot.mid() < 0.0;
            let branches_a = a.state_branches(whole).map_err(side(0))?;
            let branches_b = src.state_branches(whole).map_err(side(1))?;
            if branches_a.len() * branches_b.len() > MAX_BRANCH_PAIRS {
                let (a_t, _) = a.state(whole).map_err(side(0))?;
                let (b_t, _) = src.state(whole).map_err(side(1))?;
                let q_m = seam_quantities(&a_m, &b_m, flip);
                pieces.push((q_m, seam_quantities(&a_t, &b_t, flip), true));
            } else {
                for (key_a, a_t) in &branches_a {
                    let a_mk = a.state_in(point, *key_a).map_err(side(0))?;
                    for (key_b, b_t) in &branches_b {
                        let b_mk = src.state_in(point, *key_b).map_err(side(1))?;
                        pieces.push((
                            seam_quantities(&a_mk, &b_mk, flip),
                            seam_quantities(a_t, b_t, flip),
                            false,
                        ));
                    }
                }
            }
            vec![seam_quantities(&a_m, &b_m, flip)]
        }
        SideB::Projected(_) | SideB::Shifted { .. } => {
            let failure = |f: ProjectionFailure| match f {
                ProjectionFailure::Side(e) => InconclusiveReason::Side(1, e),
                ProjectionFailure::NotContracting => InconclusiveReason::Correspondence,
            };
            let branches_a = a.state_branches(whole).map_err(side(0))?;
            let mut witnesses = Vec::new();
            let mut flip = None;
            for (key_a, a_t) in &branches_a {
                let a_mk = a.state_in(point, *key_a).map_err(side(0))?;
                let states = match &mut *b {
                    SideB::Shifted { side, phase, dir } => side
                        .shifted_states(whole, m, *phase, *dir)
                        .map_err(failure)?,
                    SideB::Projected(projected) => projected
                        .states_along(a_t, whole, &a_mk, m)
                        .map_err(failure)?,
                    SideB::Fraction(_) => return Err(InconclusiveReason::Correspondence),
                };
                for (b_t, b_m) in &states.pieces {
                    let flip = *flip.get_or_insert_with(|| {
                        seam_quantities(&a_mk, b_m, false).normal_dot.mid() < 0.0
                    });
                    pieces.push((
                        seam_quantities(&a_mk, b_m, flip),
                        seam_quantities(a_t, b_t, flip),
                        false,
                    ));
                }
                for truth in &states.witnesses {
                    let flip = flip.unwrap_or(false);
                    witnesses.push(seam_quantities(&a_m, truth, flip));
                }
            }
            witnesses
        }
    };
    let mut jacobian = [f64::INFINITY; 2];
    for (_, q_t, _) in &pieces {
        for (held, side) in jacobian.iter_mut().zip(q_t.jacobian) {
            *held = held.min(side.lo);
        }
    }
    // NaN-safe: a Jacobian bound that is not positive (or undefined)
    // is a degenerate chart.
    if jacobian
        .iter()
        .any(|j| j.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater))
    {
        return Err(InconclusiveReason::Jacobian);
    }
    let offset = whole - point;
    let enclose = |at_m: Series, over: Series, natural_only: bool| -> Iv {
        if natural_only {
            over.c[0]
        } else {
            at_m.enclose_over(&over, offset)
        }
    };
    let norm_upper = |values: [f64; 3]| -> f64 {
        let s: f64 = values.iter().map(|v| v * v).sum();
        s.sqrt().next_up().next_up()
    };
    let max = |v: [f64; 3]| v.iter().copied().fold(0.0_f64, f64::max);
    // Upper bound: the worst branch piece; witness: the true midpoint.
    let bound = |pick: fn(&SeamQuantities) -> [Series; 3], euclidean: bool| {
        let mut upper = 0.0_f64;
        for (q_m, q_t, natural_only) in &pieces {
            let (at_m, over) = (pick(q_m), pick(q_t));
            let parts: [f64; 3] =
                std::array::from_fn(|k| enclose(at_m[k], over[k], *natural_only).mag());
            upper = upper.max(if euclidean {
                norm_upper(parts)
            } else {
                max(parts)
            });
        }
        let lower = witness_quantities
            .iter()
            .map(|q| {
                let at_m = pick(q);
                let witness: [f64; 3] = std::array::from_fn(|k| at_m[k].c[0].mig());
                if euclidean {
                    witness
                        .iter()
                        .map(|v| v * v)
                        .sum::<f64>()
                        .sqrt()
                        .next_down()
                        .max(0.0)
                } else {
                    max(witness)
                }
            })
            .fold(f64::INFINITY, f64::min);
        (upper, if lower.is_finite() { lower } else { 0.0 })
    };
    let (gap, sine, curvature) = (
        bound(|q| q.gap, true),
        bound(|q| q.normal_cross, true),
        bound(|q| q.curvature, false),
    );
    let upper = [gap.0, sine.0, curvature.0];
    if upper.iter().any(|u| !u.is_finite()) {
        return Err(InconclusiveReason::Resolution);
    }
    Ok(IntervalBounds {
        upper,
        lower: [gap.1, sine.1, curvature.1],
        jacobian,
    })
}
