//! Direct exact-boundary K-P binonic spoke-row optimizer.

use super::kp_lm_solver::{LmConfig, LmEval, run_lm};
use super::kp_spoke_normal_sine::{
    RING_SPOKE_G1_GATE, lean_penalty, ring_chord_residual, ring_spoke_senses, spoke_fold_count,
    spoke_normal_chord, whole_spoke_parameters, worst_oriented_ring_sine,
};
use popkern_core::geometry::gregory::{
    BoundedRingCrossing, GregoryPatch, RectangularGregoryBiNonic, binonic_ring_crossing_bounded,
    binonic_ring_crossing_witness, binonic_ring_feasible_start_crossing,
};

/// Linear tolerance of the PK-283 ring-level crossing gates — the same
/// value every per-patch fold gate on the promotion ladder uses, so the
/// two gate families measure the same geometry class at the same scale.
pub(in crate::construct) const RING_CROSSING_TOLERANCE: f64 = 1.0e-6;

pub(super) struct ExactBoundaryBinonicResult {
    pub binonic: Vec<RectangularGregoryBiNonic>,
    /// The candidate's residual on the metric it competes by: the
    /// honest [`ring_chord_residual`] for a chord candidate; the triple
    /// product sum for a heuristic one, which is compared like for
    /// like with the phase-2 figure and is honestly gated only after
    /// the repair.
    pub final_residual: f64,
    /// True for the triple-product heuristic's candidate.
    pub heuristic: bool,
}

/// The spoke-row L-M's objective (PK-544).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SpokeObjective {
    /// The oriented spoke normal chord: honest in the normal angle over
    /// the whole half-turn, invariant to chart skew. A step never adds
    /// a fold across a spoke.
    OrientedChord,
    /// The normalized triple product `sine · sin∠(w,u) · sin∠(w,v)`, a
    /// search heuristic only. It rewards skew and reads folds as G1, so
    /// its output reaches the promotion through the chart-orientation
    /// repair and the honest gates. It runs only when the chord
    /// candidate cannot pass the ring G1 gate: on pyramid4_uniform the
    /// chord solve stalls at a creased ring, while the heuristic's
    /// folded ring repairs to G1.
    TripleProductHeuristic,
}

pub(super) fn optimize_exact_boundary_binonic(
    initial: &[RectangularGregoryBiNonic],
) -> Option<ExactBoundaryBinonicResult> {
    if initial.is_empty() || !feasible_start(initial) {
        return None;
    }
    let senses = ring_spoke_senses(initial);
    [
        SpokeObjective::OrientedChord,
        SpokeObjective::TripleProductHeuristic,
    ]
    .into_iter()
    .find_map(|objective| optimize_spoke_rows(initial, &senses, objective))
}

fn optimize_spoke_rows(
    initial: &[RectangularGregoryBiNonic],
    senses: &[bool],
    objective: SpokeObjective,
) -> Option<ExactBoundaryBinonicResult> {
    let eval = DirectSpokeRowEval {
        base: initial.to_vec(),
        senses: senses.to_vec(),
        max_spoke_folds: spoke_fold_count(initial, senses),
        max_bands: reversal_bands(initial),
        objective,
    };
    let baseline = eval.objective_sum(initial);
    let x0 = pack_spoke_rows(initial);
    let result = run_lm(
        eval,
        x0,
        &LmConfig {
            fd_eps: 1e-4,
            lambda_init: 1e-3,
            lambda_max: 1e6,
            max_iter: 15,
            target_residual: 1e-6,
            require_feasible_start: true,
        },
        "exact-boundary-binonic-spoke-row",
    );
    if result.final_residual >= baseline {
        return None;
    }
    let mut binonic = initial.to_vec();
    unpack_spoke_rows(&result.x, &mut binonic);
    // A chord candidate over the G1 gate, or still carrying in-patch
    // reversal bands, goes through the chart repair (PK-548): unfold and
    // polish, or the in-place polish of a clean ring. The polish moves
    // the spoke interiors this solve holds fixed. If the repair fails
    // too, the heuristic objective runs next.
    if objective == SpokeObjective::OrientedChord
        && (worst_oriented_ring_sine(&binonic, senses) > RING_SPOKE_G1_GATE
            || reversal_bands(&binonic) > 0)
    {
        binonic = super::three_quad_apply::three_quad_chart_orientation::repair_chart_orientation(
            &binonic, None,
        )
        .ok()?
        .ring;
    }
    // A proven crossing declines at witness cost; the bounded verdict
    // then settles the rest without a spoke-neighbour surface solve.
    if binonic_ring_crossing_witness(&binonic, RING_CROSSING_TOLERANCE).is_some()
        || !acceptable_ring(&binonic)
    {
        return None;
    }
    let heuristic = objective == SpokeObjective::TripleProductHeuristic;
    let final_residual = if heuristic {
        result.final_residual
    } else {
        ring_chord_residual(&binonic, senses)
    };
    Some(ExactBoundaryBinonicResult {
        binonic,
        final_residual,
        heuristic,
    })
}

struct DirectSpokeRowEval {
    base: Vec<RectangularGregoryBiNonic>,
    /// Per spoke: the relative chart orientation it must keep.
    senses: Vec<bool>,
    /// Spoke folds at the start: an accepted chord step never has more.
    max_spoke_folds: usize,
    /// In-patch chart reversal bands at the start: an accepted chord step
    /// never has more (PK-548; the promotion declines banded rings).
    max_bands: usize,
    objective: SpokeObjective,
}

impl DirectSpokeRowEval {
    fn build_ring(&self, x: &[f64]) -> Vec<RectangularGregoryBiNonic> {
        let mut ring = self.base.clone();
        unpack_spoke_rows(x, &mut ring);
        ring
    }

    fn residual_of(&self, ring: &[RectangularGregoryBiNonic]) -> Vec<f64> {
        match self.objective {
            SpokeObjective::OrientedChord => chord_residual_vector(ring, &self.senses),
            SpokeObjective::TripleProductHeuristic => triple_product_residual_vector(ring),
        }
    }

    fn objective_sum(&self, ring: &[RectangularGregoryBiNonic]) -> f64 {
        self.residual_of(ring).iter().map(|v| v * v).sum()
    }
}

impl LmEval for DirectSpokeRowEval {
    fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>> {
        Some(self.residual_of(&self.build_ring(x)))
    }

    fn fold_safe(&mut self, x: &[f64]) -> bool {
        let ring = self.build_ring(x);
        if self.objective == SpokeObjective::OrientedChord
            && (spoke_fold_count(&ring, &self.senses) > self.max_spoke_folds
                || reversal_bands(&ring) > self.max_bands
                || binonic_ring_crossing_witness(&ring, RING_CROSSING_TOLERANCE).is_some())
        {
            return false;
        }
        !any_self_intersects(&ring)
    }

    // PK-275 measured decision: this eval deliberately KEEPS the FD
    // Jacobian. An analytic selector-chain Jacobian was implemented
    // and FD-agreement-verified here, but (a) it bought no wall time
    // — this optimizer's cost is the fold-safety Gregory
    // self-intersection check, not derivatives (gdb profile,
    // PK-282) — and (b) its more-accurate trajectory converged
    // pyramid4_uniform to a lower-residual ring whose ADJACENT
    // patches intersect each other, a class the per-patch fold gates
    // could not see. PK-283 closed (b): the acceptance is now
    // crossing-honest at ring level (`binonic_ring_crossing` gates
    // entry and accept above), so any trajectory landing in that
    // class is declined rather than served. (a) still holds, so FD
    // remains the keep-better choice on measurement; a future
    // analytic Jacobian is safe to reconsider purely on wall-time
    // evidence. The reverted implementation was never committed
    // (verified against 77894f2a and dangling objects, 2026-08-16).
}

fn pack_spoke_rows(patches: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    let mut x = Vec::with_capacity(patches.len() * 16 * 3);
    for patch in patches {
        for j in 1..=8 {
            push_point(&mut x, patch.grid[8][j]);
        }
        for i in 1..=8 {
            push_point(&mut x, patch.grid[i][1]);
        }
    }
    x
}

fn unpack_spoke_rows(x: &[f64], patches: &mut [RectangularGregoryBiNonic]) {
    let mut idx = 0;
    for patch in patches {
        for j in 1..=8 {
            patch.grid[8][j].x = x[idx];
            patch.grid[8][j].y = x[idx + 1];
            patch.grid[8][j].z = x[idx + 2];
            idx += 3;
        }
        for i in 1..=8 {
            patch.grid[i][1].x = x[idx];
            patch.grid[i][1].y = x[idx + 1];
            patch.grid[i][1].z = x[idx + 2];
            idx += 3;
        }
    }
}

fn push_point(x: &mut Vec<f64>, point: popkern_core::geometry::Point3) {
    x.push(point.x);
    x.push(point.y);
    x.push(point.z);
}

/// Entry feasibility: every patch individually fold-free AND no
/// macro-scale pairwise crossing (PK-283). The coarse feasible-start
/// screen — not the full accept semantics — because raw guide rings
/// legitimately carry spoke-hover contact an order deeper than
/// converged rings (measured: pyramid3_uniform's initial ring at
/// 1.886e-3, polished to 2.38e-5 by this optimizer and served clean).
fn feasible_start(binonic: &[RectangularGregoryBiNonic]) -> bool {
    if any_self_intersects(binonic) {
        return false;
    }
    if let Some(crossing) = binonic_ring_feasible_start_crossing(binonic, RING_CROSSING_TOLERANCE) {
        if popkern_core::trace::trace_enabled("kp-cap") {
            eprintln!(
                "[exact-boundary-binonic] feasible-start ring decline: {}",
                crossing.diagnostic_label()
            );
        }
        return false;
    }
    true
}

/// Accept feasibility: per-patch fold gate plus the ring-level crossing
/// verdict (PK-283 — a class the per-patch fold gates cannot see;
/// bounded under PK-547). Run once per acceptance.
fn acceptable_ring(binonic: &[RectangularGregoryBiNonic]) -> bool {
    !any_self_intersects(binonic) && ring_crossing_decline_label(binonic).is_none()
}

/// PK-283 accept-gate shared arm for the optimizers: the bounded
/// ring-crossing verdict (PK-547), with the `kp-cap` decline trace.
/// A spoke-neighbour pair the separation certificate cannot clear is
/// declined as undecided rather than settled by a surface-surface solve
/// that took up to 223 s. The optimizers hold fallback candidates, and
/// the promotion's full gate stays the acceptance authority. Returns the
/// decline label so optimizer callers can record it as a structured
/// rung outcome.
pub(super) fn ring_crossing_decline_label(binonic: &[RectangularGregoryBiNonic]) -> Option<String> {
    let label = match binonic_ring_crossing_bounded(binonic, RING_CROSSING_TOLERANCE) {
        BoundedRingCrossing::Clean => return None,
        BoundedRingCrossing::Crossing(crossing) => crossing.diagnostic_label(),
        BoundedRingCrossing::Undecided { pair } => {
            format!("ring_crossing_undecided pair=({},{})", pair.0, pair.1)
        }
    };
    if popkern_core::trace::trace_enabled("kp-cap") {
        eprintln!("[kp-ring-gate] accept decline: {label}");
    }
    Some(label)
}

fn any_self_intersects(binonic: &[RectangularGregoryBiNonic]) -> bool {
    binonic.iter().any(|patch| {
        GregoryPatch::RectangularBiNonic(Box::new(patch.clone()))
            .self_intersects_with_tolerance(1.0e-6)
    })
}

/// Per-sample oriented spoke normal chord (PK-544), followed by the
/// chart-lean barrier of both sides' cross tangents, over the whole
/// spoke (PK-547).
fn chord_residual_vector(patches: &[RectangularGregoryBiNonic], senses: &[bool]) -> Vec<f64> {
    let n = patches.len();
    let parameters = whole_spoke_parameters();
    let mut out = Vec::with_capacity(n * (3 + 2 * HOVER_OFFSETS.len()) * parameters.len());
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        let sense = senses.get(k).copied().unwrap_or(true);
        for &t in &parameters {
            let (u, v, w) = (
                a.partial_u(1.0, t),
                b.partial_v(1.0 - t, 0.0),
                a.partial_v(1.0, t),
            );
            out.push(spoke_normal_chord(u, v, w, sense));
            out.push(lean_penalty(u, w));
            out.push(lean_penalty(v, w));
            out.extend(hover_side_penalties(a, b, t, u, v, w));
        }
    }
    out
}

/// Cross-spoke parameter offsets the hover side condition samples.
const HOVER_OFFSETS: [f64; 3] = [0.03, 0.08, 0.16];

/// Hover side condition (PK-548). Near a shared spoke, each sheet's
/// points stay on its own side of the spoke's in-plane normal:
/// `(P_b(δ) − S(t)) · into_a ≤ 0`, and symmetrically for a. For a G1
/// ring this holds to first order (`into_a = −into_b`), so a violation
/// is a sheet bending back across the spoke, the local form of a
/// spoke-neighbour crossing. The objective otherwise sees only the
/// spoke, and on pyramid4_uniform the chord solve ran into such a
/// crossing 1.25e-3 off the spoke, just past the exemption band. Each
/// penalty is the normalized depth into the other side, 0 when on the
/// own side.
fn hover_side_penalties(
    a: &RectangularGregoryBiNonic,
    b: &RectangularGregoryBiNonic,
    t: f64,
    u: popkern_core::geometry::Vec3,
    v: popkern_core::geometry::Vec3,
    w: popkern_core::geometry::Vec3,
) -> impl Iterator<Item = f64> {
    let spoke_point = a.point_at(1.0, t);
    let in_plane = |x: popkern_core::geometry::Vec3| {
        let w_unit = w.normalized()?;
        (x - w_unit * x.dot(&w_unit)).normalized()
    };
    // `u` runs toward the spoke from inside a; `v` runs into b.
    let (into_a, into_b) = (in_plane(u * -1.0), in_plane(v));
    let depth = move |p: popkern_core::geometry::Point3,
                      into: Option<popkern_core::geometry::Vec3>| {
        let offset = p - spoke_point;
        match into {
            Some(into) if offset.length() > 0.0 => (offset.dot(&into) / offset.length()).max(0.0),
            _ => 0.0,
        }
    };
    let (a, b) = (a.clone(), b.clone());
    HOVER_OFFSETS.into_iter().flat_map(move |delta| {
        [
            depth(b.point_at(1.0 - t, delta), into_a),
            depth(a.point_at(1.0 - delta, t), into_b),
        ]
    })
}

/// Per-sample normalized triple product (the search heuristic).
fn triple_product_residual_vector(patches: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    let n = patches.len();
    let mut out = Vec::with_capacity(n * 50);
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            let (u, v, w) = (
                a.partial_u(1.0, t),
                b.partial_v(1.0 - t, 0.0),
                a.partial_v(1.0, t),
            );
            let denom = u.length() * v.length() * w.length();
            let numer = u.cross(&v).dot(&w).abs();
            out.push(if denom > 1.0e-300 { numer / denom } else { 0.0 });
        }
    }
    out
}

#[cfg(test)]
#[path = "kp_exact_boundary_ring_gate_tests.rs"]
mod kp_exact_boundary_ring_gate_tests;

/// In-patch chart reversal bands of a ring, with the repair's chart
/// margin (`CHART_RAIL_MARGIN`).
fn reversal_bands(ring: &[RectangularGregoryBiNonic]) -> usize {
    super::three_quad_apply::three_quad_chart_orientation::binonic_margin_band_count(ring)
}
