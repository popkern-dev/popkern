//! Canonical-forward audit (PK-478 stage 1).
//!
//! PK-303 convention: an edge's CANONICAL half-edge (`edge.half_edge`)
//! traverses the stored curve in increasing t. On OPEN edges the
//! convention is enforced by endpoints; on CLOSED edges (rim circles,
//! periodic analytic rings) nothing checked it, and producers have
//! been measured violating it (the PK-465 bore wall; the fresh
//! cylinder's top-rim F/F pcurve pair). The boolean's closed-edge
//! side witness exists to paper over those lies; the PK-478 campaign
//! retires it by making this invariant validated and producer-
//! enforced, then registering blindly (see
//! docs/architecture/pk478-sense-aware-merge.md).
//!
//! The audit recovers each half-edge's loop-traversal direction from
//! its OWN pcurve trace mapped through its face's surface, compares
//! twin claims (they must oppose on a manifold), and resolves an
//! agreeing pair (one liar) by UV chain continuity with the loop's
//! neighbouring pcurves. It measures; it does not repair.

use crate::error::ValidationError;
use crate::geometry::SurfaceEval;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

/// Indeterminate reason for an agreeing-claims pair neither the
/// chain tiebreak nor (afterwards) the rim-pair arbitration could
/// resolve; the arbitration pass keys on it.
pub(super) const CHAIN_AMBIGUOUS: &str = "agreeing claims, chain tiebreak ambiguous";

/// Per-closed-edge verdict from [`BRep::audit_canonical_forward`].
#[derive(Debug, Clone)]
pub struct CanonicalForwardVerdict {
    /// The audited closed edge.
    pub edge: EdgeHandle,
    /// The verdict.
    pub state: CanonicalForwardState,
}

/// Audit outcome for one closed edge.
#[derive(Debug, Clone, PartialEq)]
pub enum CanonicalForwardState {
    /// Twin pcurve claims oppose and the canonical half runs the
    /// stored curve forward.
    Ok,
    /// Twin pcurve claims oppose and the canonical half runs the
    /// stored curve BACKWARD — the PK-303 convention is violated.
    Violated,
    /// Twin claims agree (one pcurve lies about its direction); the
    /// chain-continuity tiebreak identified the liar and the trusted
    /// side says the canonical convention holds.
    OkByChain {
        /// The half-edge whose pcurve direction is the liar.
        lying_half_edge: HalfEdgeHandle,
    },
    /// Twin claims agree and the trusted side says the canonical
    /// convention is violated.
    ViolatedByChain {
        /// The half-edge whose pcurve direction is the liar.
        lying_half_edge: HalfEdgeHandle,
    },
    /// Not enough evidence: missing pcurves/curves/surfaces, both
    /// sides chain-ambiguous, or degenerate sampling.
    Indeterminate(&'static str),
}

impl BRep {
    /// Audit the canonical-forward convention on every CLOSED,
    /// non-seam edge. Seam edges (both half-edges on one face) keep
    /// their chart-branch handling and are exempt, matching the
    /// merge's closed-edge witness exemption.
    pub fn audit_canonical_forward(&self) -> Vec<CanonicalForwardVerdict> {
        let winding: std::collections::HashMap<crate::topology::FaceHandle, Option<bool>> = self
            .audit_winding_honesty()
            .into_iter()
            .map(|v| (v.face, v.honest))
            .collect();
        let mut out = Vec::new();
        for (eh, edge) in &self.edges {
            let Some(he) = self.half_edges.get(edge.half_edge) else {
                continue;
            };
            let Some(twin) = self.half_edges.get(he.twin) else {
                continue;
            };
            if he.origin != twin.origin {
                continue; // open edge: endpoint-enforced
            }
            let face_of = |lh| self.loops.get(lh).map(|l| l.face);
            let (f_can, f_twin) = (face_of(he.loop_), face_of(twin.loop_));
            if f_can.is_none() || f_twin.is_none() || f_can == f_twin {
                continue; // dangling or seam
            }
            let state = self.closed_edge_state(edge.half_edge, he.twin, &winding);
            out.push(CanonicalForwardVerdict { edge: eh, state });
        }
        // Pass 2 (stage 3): seam-adjacent rim junctions defeat the
        // chain tiebreak (all endpoints coincide mod period); the
        // rim-pair loop invariant arbitrates them from a determinate
        // sibling rim.
        self.arbitrate_rim_pairs(&mut out);
        out
    }

    /// Repair every measurable canonical-forward violation in place
    /// (PK-478 stage 2 — the explicit-repair arm of the invariant):
    ///
    /// - `Violated` (claims oppose, canonical backward): swap the
    ///   edge's canonical half-edge to its twin. Restores the PK-303
    ///   assumption every canonical consumer makes; trim geometry is
    ///   untouched.
    /// - `ViolatedByChain` (one pcurve lies): reverse the convicted
    ///   pcurve so twin claims oppose, then make the forward-running
    ///   side canonical. The UV image (the trim) is unchanged.
    /// - `Ok`/`OkByChain`/`Indeterminate`: untouched — the repair
    ///   never acts on ambiguity.
    ///
    /// Returns the number of edges repaired. Intended for PRODUCER
    /// exits (the blend's closed-rim strips were the measured lying
    /// population); running it on consumer inputs would mask
    /// producer defects the audit exists to surface.
    pub fn repair_canonical_forward(&mut self) -> usize {
        let verdicts = self.audit_canonical_forward();
        let mut repaired = 0;
        let mut pcurve_reversed = false;
        for v in verdicts {
            match v.state {
                CanonicalForwardState::Violated => {
                    let canonical = self.edges[v.edge].half_edge;
                    let twin = self.half_edges[canonical].twin;
                    self.edges[v.edge].half_edge = twin;
                    repaired += 1;
                }
                CanonicalForwardState::ViolatedByChain { lying_half_edge } => {
                    if let Some(pc) = self
                        .half_edges
                        .get(lying_half_edge)
                        .and_then(|h| h.pcurve.clone())
                    {
                        self.half_edges[lying_half_edge].pcurve = Some(pc.reversed());
                        pcurve_reversed = true;
                    }
                    // After the pcurve repair the claims oppose; the
                    // trusted (non-lying) side kept its claim. Make
                    // whichever side now claims FORWARD canonical.
                    let canonical = self.edges[v.edge].half_edge;
                    let twin = self.half_edges[canonical].twin;
                    let can_claim = self.pcurve_traversal_claim(canonical);
                    match can_claim {
                        Some(true) => {}
                        Some(false) => self.edges[v.edge].half_edge = twin,
                        None => continue,
                    }
                    repaired += 1;
                }
                _ => {}
            }
        }
        // A reversed pcurve flips its loop's measured UV winding;
        // the owning face's orientation flag was calibrated against
        // the LIE and must follow the truth or the face serves
        // winding-dishonest (measured: the bossed housing's shrunk
        // cap after the fillet-rim repair). `reconcile` only touches
        // measured-dishonest faces, so on a producer-fresh body this
        // is scoped to exactly the repaired caps.
        if pcurve_reversed {
            self.reconcile_orientation_flags();
        }
        repaired
    }

    /// Validation-style wrapper: one error per violated edge.
    pub fn validate_canonical_forward_all(&self) -> Vec<ValidationError> {
        self.audit_canonical_forward()
            .into_iter()
            .filter(|v| {
                matches!(
                    v.state,
                    CanonicalForwardState::Violated | CanonicalForwardState::ViolatedByChain { .. }
                )
            })
            .map(|v| ValidationError::CanonicalBackward { edge: v.edge })
            .collect()
    }

    fn closed_edge_state(
        &self,
        canonical: HalfEdgeHandle,
        twin: HalfEdgeHandle,
        winding: &std::collections::HashMap<crate::topology::FaceHandle, Option<bool>>,
    ) -> CanonicalForwardState {
        use CanonicalForwardState as S;
        let claim_can = self.pcurve_traversal_claim(canonical);
        let claim_twin = self.pcurve_traversal_claim(twin);
        match (claim_can, claim_twin) {
            (Some(fc), Some(ft)) if fc != ft => {
                if fc {
                    S::Ok
                } else {
                    S::Violated
                }
            }
            (Some(fc), Some(_ft)) => {
                // Agreeing claims: one pcurve lies. Trust the side
                // whose pcurve is UV-chain-continuous with its loop
                // neighbours; a mis-directed pcurve has both its
                // endpoints against the chain. The shared one-liar
                // decision matrix lives with the rim-pair
                // arbitration, which reuses it in pass 2.
                let cont_can = self.pcurve_chain_continuous(canonical);
                let cont_twin = self.pcurve_chain_continuous(twin);
                // Winding witness: a loop reversed as a whole stays
                // chain-continuous, so the chain test passes it and
                // convicts the honest partner (PK-741: a housing wall
                // whose rims and seams were re-authored around a lying
                // rim convicted both honest caps). Such a face's loops
                // wind against its own orientation flag. When exactly
                // one side's face is winding-dishonest, its pcurve is
                // the liar; otherwise the chain decides.
                let honest_face = |heh: HalfEdgeHandle| {
                    let face = self.loops.get(self.half_edges.get(heh)?.loop_)?.face;
                    winding.get(&face).copied().flatten()
                };
                let (cont_can, cont_twin) = match (honest_face(canonical), honest_face(twin)) {
                    (Some(false), Some(true)) => (Some(false), Some(true)),
                    (Some(true), Some(false)) => (Some(true), Some(false)),
                    _ => (cont_can, cont_twin),
                };
                super::canonical_forward_arbitration::resolve_one_liar(
                    fc, cont_can, cont_twin, canonical, twin,
                )
                .unwrap_or(S::Indeterminate(CHAIN_AMBIGUOUS))
            }
            _ => S::Indeterminate("pcurve traversal unmeasurable"),
        }
    }

    /// `Some(true)` when `heh`'s pcurve traces the stored edge curve
    /// in increasing t (accounting for periodic wrap), `Some(false)`
    /// for decreasing, `None` when unmeasurable.
    pub(crate) fn pcurve_traversal_claim(&self, heh: HalfEdgeHandle) -> Option<bool> {
        let he = self.half_edges.get(heh)?;
        let pc = he.pcurve.as_ref()?;
        let face = self.loops.get(he.loop_)?.face;
        let sof = self.faces.get(face)?.surface.as_ref()?;
        let curve = self.edges.get(he.edge)?.curve.as_ref()?;
        let (t0, t1) = pc.parameter_range();
        // A DECREASING range is a legal reversed representation
        // (evaluation normalizes t into the span, traversal is
        // t_start->t_end either way); only a degenerate span is
        // unmeasurable. Measured: the blend caps author their rim
        // pcurves with decreasing angle ranges.
        if !(t0.is_finite() && t1.is_finite()) || t1 == t0 {
            return None;
        }
        let tm = t0 + (t1 - t0) * 0.5;
        let te = t0 + (t1 - t0) * 0.505;
        let a = pc.point_at(tm);
        let b = pc.point_at(te);
        let p1 = sof.geometry.point_at(a.x, a.y);
        let p2 = sof.geometry.point_at(b.x, b.y);
        let step_len = (p2 - p1).length();
        if step_len <= 0.0 {
            return None;
        }
        // Exact inversion first; an SSI-approximated pcurve trace
        // sits slightly off the stored curve (measured 1.6e-5 on the
        // tilted bore rim, against the inversion's on-curve gate), so
        // fall back to a projection accepted only when the deviation
        // is small against the sampling step — direction ordering
        // survives any deviation well under half a step.
        let param_of = |p: crate::geometry::Point3| -> Option<f64> {
            match crate::geometry::curve_inversion::invert_curve_at_point(curve, p, &self.tolerance)
            {
                Ok(t) => Some(t),
                Err(_) => {
                    let (t, dist) =
                        crate::geometry::curve_nearest_param::nearest_param_on_curve(curve, &p);
                    (dist <= 0.25 * step_len).then_some(t)
                }
            }
        };
        let c1 = param_of(p1)?;
        let c2 = param_of(p2)?;
        let span = curve.t_end - curve.t_start;
        if span.abs() <= 0.0 {
            return None;
        }
        // Minimal-image signed difference (the curve is closed).
        let mut d = c2 - c1;
        d -= span * (d / span).round();
        if d.abs() < 1e-12 * span.abs() {
            return None;
        }
        Some(d > 0.0)
    }

    /// `Some(true)` when `heh`'s pcurve endpoints connect head-to-tail
    /// with its loop neighbours' pcurve endpoints in UV; `Some(false)`
    /// when both connections fail with the endpoints SWAPPED (the
    /// mis-direction signature); `None` for self-loops, missing
    /// neighbour pcurves, or mixed evidence.
    fn pcurve_chain_continuous(&self, heh: HalfEdgeHandle) -> Option<bool> {
        let he = self.half_edges.get(heh)?;
        if he.next == heh {
            return None; // self-loop: no neighbours
        }
        let pc = he.pcurve.as_ref()?;
        let next_pc = self.half_edges.get(he.next)?.pcurve.as_ref()?;
        let prev_pc = self.half_edges.get(he.prev)?.pcurve.as_ref()?;
        let (t0, t1) = pc.parameter_range();
        let start = pc.point_at(t0);
        let end = pc.point_at(t1);
        let (n0, _) = next_pc.parameter_range();
        let next_start = next_pc.point_at(n0);
        let (_, p1) = prev_pc.parameter_range();
        let prev_end = prev_pc.point_at(p1);
        // Periodic-aware UV distance: seam-adjacent pcurves on a
        // periodic chart legitimately connect a full period apart
        // (the aliased-seam conventions).
        let (u_period, v_period) = self
            .loops
            .get(he.loop_)
            .and_then(|l| self.faces.get(l.face))
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.periodicity())
            .unwrap_or((None, None));
        let fold = |mut delta: f64, period: Option<f64>| -> f64 {
            if let Some(p) = period
                && p > 0.0
            {
                delta -= p * (delta / p).round();
            }
            delta
        };
        let d_raw = |a: crate::geometry::Point3, b: crate::geometry::Point3| {
            ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt()
        };
        let d_folded = |a: crate::geometry::Point3, b: crate::geometry::Point3| {
            let du = fold(a.x - b.x, u_period);
            let dv = fold(a.y - b.y, v_period);
            (du * du + dv * dv).sqrt()
        };
        // Chart-scale epsilon from this trio of pcurves.
        let scale = d_raw(start, end).max(d_raw(prev_end, next_start)).max(1e-9);
        let eps = 1e-6 * (1.0 + scale);
        let verdict =
            |d: &dyn Fn(crate::geometry::Point3, crate::geometry::Point3) -> f64| -> (bool, bool) {
                (
                    d(end, next_start) < eps && d(start, prev_end) < eps,
                    d(start, next_start) < eps && d(end, prev_end) < eps,
                )
            };
        // Exact adjacency first — a mis-directed pcurve's endpoint
        // swap often differs by exactly one period, and folding
        // erases that evidence (measured: the chamfer tip rim's
        // convicting swap vanished under an unconditional fold).
        // Fold only when the exact test sees nothing, for the
        // legitimate seam-adjacent period offsets.
        match verdict(&d_raw) {
            (true, false) => Some(true),
            (false, true) => Some(false),
            _ => match verdict(&d_folded) {
                (true, false) => Some(true),
                (false, true) => Some(false),
                _ => None,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::CanonicalForwardState;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// The audit is reversal-STABLE on a constructor cylinder: the
    /// same verdicts before and after `reverse_body_orientation`
    /// (the convention is loop-covariant; reversal swaps canonical
    /// bits and pcurve directions together). Since PK-478 stage 3,
    /// the audit is also fully DETERMINATE here: the bottom rim's
    /// seam-adjacent junction ambiguity (all endpoint candidates
    /// coincide modulo the period) is resolved by full-period
    /// rim-pair arbitration — the determinate top rim decides its
    /// sibling through the wall loop's opposite-direction invariant.
    #[test]
    fn cylinder_audit_is_reversal_stable() {
        let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder");
        let verdicts = cyl.audit_canonical_forward();
        assert!(!verdicts.is_empty(), "cylinder has closed rim edges");
        assert!(
            verdicts
                .iter()
                .any(|v| matches!(v.state, CanonicalForwardState::Ok)),
            "at least one rim resolves Ok: {verdicts:?}"
        );
        assert!(
            !verdicts
                .iter()
                .any(|v| matches!(v.state, CanonicalForwardState::Indeterminate(_))),
            "stage-3 arbitration leaves no indeterminate rims: {verdicts:?}"
        );
        assert!(
            !verdicts.iter().any(|v| matches!(
                v.state,
                CanonicalForwardState::Violated | CanonicalForwardState::ViolatedByChain { .. }
            )),
            "constructor cylinder must not violate canonical-forward: {verdicts:?}"
        );
        let mut reversed = cyl.clone();
        reversed.reverse_body_orientation();
        let rv = reversed.audit_canonical_forward();
        let states: Vec<_> = verdicts.iter().map(|v| format!("{:?}", v.state)).collect();
        let rstates: Vec<_> = rv.iter().map(|v| format!("{:?}", v.state)).collect();
        assert_eq!(states, rstates, "verdicts must be reversal-stable");
    }
}
