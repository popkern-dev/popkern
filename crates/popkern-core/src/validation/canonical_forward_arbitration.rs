//! Full-period rim-pair arbitration for the canonical-forward audit
//! (PK-478 stage 3).
//!
//! The chain-continuity tiebreak is void at SEAM-ADJACENT rim
//! junctions: on a periodic wall chart every endpoint candidate
//! coincides modulo the period, so the forward and the swapped
//! connection both match and no side convicts (the constructor
//! cylinder's bottom rim, measured in the stage-1 census). The
//! arbitration recovers evidence from a LOOP invariant instead of a
//! junction one: a loop that traverses two full-period
//! constant-parameter rims bounds the material band between them,
//! and a band boundary must run its two rims in OPPOSITE chart
//! directions — a same-direction pair has cancelling shoelace
//! contributions and cannot wind once around the band (such a face
//! would be winding-dishonest). A rim whose edge verdict is already
//! determinate carries a trusted traversal direction; its sibling's
//! pcurve is honest exactly when the two running directions oppose.
//!
//! Tolerances are UV-relative (period-scaled for the running
//! parameter, magnitude-scaled for the held one) — chart space has
//! no global length unit.

use std::collections::HashMap;

use super::canonical_forward::CHAIN_AMBIGUOUS;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};
use crate::validation::{CanonicalForwardState, CanonicalForwardVerdict};

/// A pcurve that runs the full chart period in one parameter while
/// holding the other constant — the UV signature of a wall rim.
struct RimLine {
    /// True when the RUNNING parameter is u (a constant-v rim).
    runs_u: bool,
    /// The held-constant parameter value.
    const_val: f64,
    /// True when the running parameter increases along the pcurve.
    increasing: bool,
}

/// The one-liar decision matrix shared by the chain tiebreak and the
/// rim-pair arbitration. `fc` is the AGREED (twin-matching) claim of
/// both pcurves; `honest_can`/`honest_twin` say whether each side's
/// pcurve is loop-consistent (`None` = no evidence). Exactly one
/// side lies; the trusted side's claim decides the canonical bit.
pub(super) fn resolve_one_liar(
    fc: bool,
    honest_can: Option<bool>,
    honest_twin: Option<bool>,
    canonical: HalfEdgeHandle,
    twin: HalfEdgeHandle,
) -> Option<CanonicalForwardState> {
    use CanonicalForwardState as S;
    match (honest_can, honest_twin) {
        // One-sided evidence is accepted: a chain-backed (or
        // rim-backed) side is trusted when the other is unmeasurable,
        // and a convicted side condemns itself the same way.
        (Some(true), None) | (None, Some(false)) | (Some(true), Some(false)) => Some(if fc {
            S::OkByChain {
                lying_half_edge: twin,
            }
        } else {
            S::ViolatedByChain {
                lying_half_edge: twin,
            }
        }),
        // Canonical side lies: the trusted twin's claim equals fc
        // (agreeing pair), so the twin truly runs fc and the
        // canonical truly runs !fc.
        (None, Some(true)) | (Some(false), None) | (Some(false), Some(true)) => Some(if fc {
            S::ViolatedByChain {
                lying_half_edge: canonical,
            }
        } else {
            S::OkByChain {
                lying_half_edge: canonical,
            }
        }),
        _ => None,
    }
}

impl BRep {
    /// PK-478 stage-3 pass 2: re-resolve `Indeterminate(chain
    /// tiebreak ambiguous)` verdicts by rim-pair arbitration, to a
    /// fixpoint (a newly determinate rim can arbitrate its own
    /// sibling on another loop).
    pub(super) fn arbitrate_rim_pairs(&self, verdicts: &mut [CanonicalForwardVerdict]) {
        loop {
            let states: HashMap<EdgeHandle, CanonicalForwardState> =
                verdicts.iter().map(|v| (v.edge, v.state.clone())).collect();
            let mut changed = false;
            for v in verdicts.iter_mut() {
                let ambiguous = matches!(
                    &v.state,
                    CanonicalForwardState::Indeterminate(r) if *r == CHAIN_AMBIGUOUS
                );
                if !ambiguous {
                    continue;
                }
                let Some(edge) = self.edges.get(v.edge) else {
                    continue;
                };
                let canonical = edge.half_edge;
                let Some(he) = self.half_edges.get(canonical) else {
                    continue;
                };
                let twin = he.twin;
                let Some(fc) = self.pcurve_traversal_claim(canonical) else {
                    continue;
                };
                let honest_can = self.rim_pair_evidence(canonical, &states);
                let honest_twin = self.rim_pair_evidence(twin, &states);
                if let Some(state) = resolve_one_liar(fc, honest_can, honest_twin, canonical, twin)
                {
                    v.state = state;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }

    /// `Some(true)` when `heh`'s pcurve runs its full-period rim in
    /// the loop-consistent direction (opposite its determinate
    /// sibling rim on the same loop), `Some(false)` when it runs
    /// WITH the sibling (the mis-direction signature), `None` when
    /// the side is not a full-period rim or no determinate sibling
    /// exists on its loop.
    fn rim_pair_evidence(
        &self,
        heh: HalfEdgeHandle,
        states: &HashMap<EdgeHandle, CanonicalForwardState>,
    ) -> Option<bool> {
        let line = self.full_period_rim_line(heh)?;
        let he = self.half_edges.get(heh)?;
        let mut cur = he.next;
        let mut steps = 0usize;
        while cur != heh && steps <= self.half_edges.len() {
            steps += 1;
            let sib_he = self.half_edges.get(cur)?;
            if sib_he.edge != he.edge
                && let Some(sib) = self.full_period_rim_line(cur)
                && sib.runs_u == line.runs_u
                && (sib.const_val - line.const_val).abs()
                    > 1e-6 * (1.0 + sib.const_val.abs() + line.const_val.abs())
                && let Some(sib_true_increasing) =
                    trusted_rim_direction(states.get(&sib_he.edge), cur, &sib)
            {
                return Some(line.increasing != sib_true_increasing);
            }
            cur = sib_he.next;
        }
        None
    }

    /// The rim-line signature of `heh`'s pcurve, or `None` when the
    /// pcurve does not run a full chart period with the other
    /// parameter held constant.
    fn full_period_rim_line(&self, heh: HalfEdgeHandle) -> Option<RimLine> {
        let he = self.half_edges.get(heh)?;
        let pc = he.pcurve.as_ref()?;
        let face = self.loops.get(he.loop_)?.face;
        let sof = self.faces.get(face)?.surface.as_ref()?;
        let (u_period, v_period) = sof.geometry.periodicity();
        let (t0, t1) = pc.parameter_range();
        // Decreasing ranges are legal reversed representations; see
        // `pcurve_traversal_claim`.
        if !(t0.is_finite() && t1.is_finite()) || t1 == t0 {
            return None;
        }
        let a = pc.point_at(t0);
        let m = pc.point_at(0.5 * (t0 + t1));
        let b = pc.point_at(t1);
        let line = |run: [f64; 3], held: [f64; 3], period: Option<f64>| -> Option<(f64, bool)> {
            let p = period?;
            if p <= 0.0 {
                return None;
            }
            let span = run[2] - run[0];
            if (span.abs() - p).abs() > 1e-6 * p {
                return None;
            }
            let held_dev = (held[0] - held[1]).abs().max((held[2] - held[1]).abs());
            if held_dev > 1e-6 * (1.0 + p) {
                return None;
            }
            // The midpoint must sit mid-run: a doubled-back pcurve
            // can still land its endpoints one period apart.
            if (run[1] - 0.5 * (run[0] + run[2])).abs() > 1e-3 * p {
                return None;
            }
            Some((held[1], span > 0.0))
        };
        if let Some((const_val, increasing)) = line([a.x, m.x, b.x], [a.y, m.y, b.y], u_period) {
            return Some(RimLine {
                runs_u: true,
                const_val,
                increasing,
            });
        }
        if let Some((const_val, increasing)) = line([a.y, m.y, b.y], [a.x, m.x, b.x], v_period) {
            return Some(RimLine {
                runs_u: false,
                const_val,
                increasing,
            });
        }
        None
    }
}

/// The TRUE running direction of a determinate sibling rim: its
/// stored pcurve direction, negated when that very half-edge was the
/// convicted liar. `None` for indeterminate siblings.
fn trusted_rim_direction(
    state: Option<&CanonicalForwardState>,
    sib_heh: HalfEdgeHandle,
    sib: &RimLine,
) -> Option<bool> {
    use CanonicalForwardState as S;
    match state? {
        S::Ok | S::Violated => Some(sib.increasing),
        S::OkByChain { lying_half_edge } | S::ViolatedByChain { lying_half_edge } => {
            if *lying_half_edge == sib_heh {
                Some(!sib.increasing)
            } else {
                Some(sib.increasing)
            }
        }
        S::Indeterminate(_) => None,
    }
}
