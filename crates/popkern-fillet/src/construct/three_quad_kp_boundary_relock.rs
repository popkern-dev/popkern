//! Transactional K-P outer-boundary relock candidates.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{GregoryPatch, HexBoundaryData, RectangularGregoryBiNonic};

use super::three_quad_kp_promote::{worst_ring_coplanarity, worst_ring_g0_gap};

#[derive(Clone, Copy)]
enum RelockMode {
    BoundaryLinearFrame,
    PreserveFirstDerivative,
    PreserveSecondDerivative,
}

pub(super) fn accepted_boundary_relock(
    binonic: &[RectangularGregoryBiNonic],
    boundaries: &[HexBoundaryData],
    pre_worst: f64,
) -> Option<Vec<RectangularGregoryBiNonic>> {
    let mut best_solution: Option<ScoredRelock> = None;
    for mode in [
        RelockMode::BoundaryLinearFrame,
        RelockMode::PreserveFirstDerivative,
        RelockMode::PreserveSecondDerivative,
    ] {
        if let Some(solution) = optimize_boundary_relock(binonic, boundaries, pre_worst, mode)
            && best_solution
                .as_ref()
                .is_none_or(|best| solution.score < best.score)
        {
            best_solution = Some(solution);
        }
    }
    best_solution.map(|solution| solution.patches)
}

#[derive(Clone)]
struct RelockState {
    row0_alpha: Vec<f64>,
    col9_alpha: Vec<f64>,
    ring_recovery: f64,
}

struct ScoredRelock {
    patches: Vec<RectangularGregoryBiNonic>,
    score: f64,
    boundary_residual: f64,
}

fn optimize_boundary_relock(
    baseline: &[RectangularGregoryBiNonic],
    boundaries: &[HexBoundaryData],
    pre_worst: f64,
    mode: RelockMode,
) -> Option<ScoredRelock> {
    if baseline.len() != boundaries.len() {
        return None;
    }
    let mut state = RelockState {
        row0_alpha: vec![0.0; baseline.len()],
        col9_alpha: vec![0.0; baseline.len()],
        ring_recovery: 1.0,
    };
    let mut best = score_candidate(baseline.to_vec(), boundaries, pre_worst)?;
    if best.boundary_residual <= BOUNDARY_RESIDUAL_GATE {
        return Some(best);
    }

    let mut improved = true;
    while improved {
        improved = false;
        for step in [0.25, 0.125, 0.0625, 0.03125, 0.015625] {
            for idx in 0..baseline.len() {
                let mut trial = state.clone();
                trial.row0_alpha[idx] = (trial.row0_alpha[idx] + step).min(1.0);
                if trial.row0_alpha[idx] > state.row0_alpha[idx]
                    && let Some((candidate_state, candidate)) = best_candidate_for_state(
                        baseline, boundaries, pre_worst, trial, mode, &best,
                    )
                    && candidate_improves(&candidate, &best)
                {
                    state = candidate_state;
                    best = candidate;
                    improved = true;
                }

                let mut trial = state.clone();
                trial.col9_alpha[idx] = (trial.col9_alpha[idx] + step).min(1.0);
                if trial.col9_alpha[idx] > state.col9_alpha[idx]
                    && let Some((candidate_state, candidate)) = best_candidate_for_state(
                        baseline, boundaries, pre_worst, trial, mode, &best,
                    )
                    && candidate_improves(&candidate, &best)
                {
                    state = candidate_state;
                    best = candidate;
                    improved = true;
                }
            }
        }
    }

    (best.boundary_residual <= BOUNDARY_RESIDUAL_GATE).then_some(best)
}

fn best_candidate_for_state(
    baseline: &[RectangularGregoryBiNonic],
    boundaries: &[HexBoundaryData],
    pre_worst: f64,
    mut state: RelockState,
    mode: RelockMode,
    current: &ScoredRelock,
) -> Option<(RelockState, ScoredRelock)> {
    let boundary_probe = build_candidate_from_state(baseline, boundaries, &state, mode);
    let probe_residual = boundary_residual_for_binonic(&boundary_probe, boundaries)?;
    if probe_residual > BOUNDARY_RESIDUAL_GATE
        && probe_residual >= current.boundary_residual - 1.0e-9
    {
        return None;
    }

    let mut best: Option<(RelockState, ScoredRelock)> = None;
    for recovery in [1.0, 0.875, 0.5, 0.0] {
        state.ring_recovery = recovery;
        let candidate = build_candidate_from_state(baseline, boundaries, &state, mode);
        if let Some(scored) = score_candidate(candidate, boundaries, pre_worst)
            && best
                .as_ref()
                .is_none_or(|(_, current)| scored.score < current.score)
        {
            best = Some((state.clone(), scored));
        }
    }
    best
}

fn candidate_improves(candidate: &ScoredRelock, current: &ScoredRelock) -> bool {
    candidate.boundary_residual < current.boundary_residual - 1.0e-9
        || candidate.score < current.score * (1.0 - 1.0e-9)
}

fn build_candidate_from_state(
    baseline: &[RectangularGregoryBiNonic],
    boundaries: &[HexBoundaryData],
    state: &RelockState,
    mode: RelockMode,
) -> Vec<RectangularGregoryBiNonic> {
    let mut candidate = baseline.to_vec();
    relock_outer_boundaries(&mut candidate, baseline, boundaries, state, mode);
    recover_shared_spoke_rows(&mut candidate, baseline, state.ring_recovery);
    fill_interiors_from_relocked_frame(&mut candidate);
    candidate
}

fn relock_outer_boundaries(
    binonic: &mut [RectangularGregoryBiNonic],
    baseline: &[RectangularGregoryBiNonic],
    boundaries: &[HexBoundaryData],
    state: &RelockState,
    mode: RelockMode,
) {
    for (idx, (patch, boundary)) in binonic.iter_mut().zip(boundaries.iter()).enumerate() {
        let original = &baseline[idx];
        relock_row0(patch, original, boundary, state.row0_alpha[idx], mode);
        relock_col9(patch, original, boundary, state.col9_alpha[idx], mode);
    }
}

fn relock_row0(
    patch: &mut RectangularGregoryBiNonic,
    original: &RectangularGregoryBiNonic,
    boundary: &HexBoundaryData,
    alpha: f64,
    mode: RelockMode,
) {
    for j in 1..9 {
        let old_row0 = original.grid[0][j];
        let old_row1 = original.grid[1][j];
        let old_row2 = original.grid[2][j];
        let new_row0 = lerp_point(old_row0, boundary.row0[j], alpha);
        patch.grid[0][j] = new_row0;
        if matches!(
            mode,
            RelockMode::PreserveFirstDerivative | RelockMode::PreserveSecondDerivative
        ) {
            let first = point_delta(old_row1, old_row0);
            patch.grid[1][j] = add_delta(new_row0, first);
        }
        if matches!(mode, RelockMode::PreserveSecondDerivative) {
            let second = second_delta(old_row2, old_row1, old_row0);
            patch.grid[2][j] = add_delta(
                add_delta(
                    new_row0,
                    scale_delta(first_delta(patch.grid[1][j], new_row0), 2.0),
                ),
                second,
            );
        } else if matches!(mode, RelockMode::BoundaryLinearFrame) {
            patch.grid[1][j] = lerp_point(patch.grid[0][j], original.grid[9][j], 1.0 / 9.0);
            patch.grid[2][j] = lerp_point(patch.grid[0][j], original.grid[9][j], 2.0 / 9.0);
        }
    }
}

fn relock_col9(
    patch: &mut RectangularGregoryBiNonic,
    original: &RectangularGregoryBiNonic,
    boundary: &HexBoundaryData,
    alpha: f64,
    mode: RelockMode,
) {
    for i in 1..9 {
        let old_col9 = original.grid[i][9];
        let old_col8 = original.grid[i][8];
        let old_col7 = original.grid[i][7];
        let new_col9 = lerp_point(old_col9, boundary.col9[i], alpha);
        patch.grid[i][9] = new_col9;
        if i < 8
            && matches!(
                mode,
                RelockMode::PreserveFirstDerivative | RelockMode::PreserveSecondDerivative
            )
        {
            let first = point_delta(old_col9, old_col8);
            patch.grid[i][8] = subtract_delta(new_col9, first);
        }
        if i < 7 && matches!(mode, RelockMode::PreserveSecondDerivative) {
            let second = second_delta(old_col9, old_col8, old_col7);
            patch.grid[i][7] = add_delta(
                subtract_delta(patch.grid[i][8], first_delta(new_col9, patch.grid[i][8])),
                second,
            );
        } else if i < 8 && matches!(mode, RelockMode::BoundaryLinearFrame) {
            patch.grid[i][8] = lerp_point(original.grid[i][0], patch.grid[i][9], 8.0 / 9.0);
            if i < 7 {
                patch.grid[i][7] = lerp_point(original.grid[i][0], patch.grid[i][9], 7.0 / 9.0);
            }
        }
    }
}

const BOUNDARY_RESIDUAL_GATE: f64 = 5.0e-3;
const RING_G0_GATE: f64 = 1.0e-10;
const RING_G1_GATE: f64 = 5.0e-2;

fn score_candidate(
    candidate: Vec<RectangularGregoryBiNonic>,
    boundaries: &[HexBoundaryData],
    pre_worst: f64,
) -> Option<ScoredRelock> {
    let candidate_gp = candidate
        .iter()
        .cloned()
        .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
        .collect::<Vec<_>>();
    if candidate_gp
        .iter()
        .any(|gp| gp.self_intersects_with_tolerance(1.0e-6))
    {
        return None;
    }
    let g0 = worst_ring_g0_gap(&candidate_gp);
    if g0 > RING_G0_GATE {
        return None;
    }
    let post_worst = worst_ring_coplanarity(&candidate_gp);
    if post_worst > RING_G1_GATE {
        return None;
    }
    if pre_worst > 1.0e-12 && post_worst > 1.2 * pre_worst {
        return None;
    }
    let boundary_residual = sampled_boundary_residual(&candidate_gp, boundaries)?;
    let normalized_boundary = boundary_residual / BOUNDARY_RESIDUAL_GATE;
    let normalized_g1 = post_worst / RING_G1_GATE;
    Some(ScoredRelock {
        patches: candidate,
        score: normalized_boundary * normalized_boundary + 0.05 * normalized_g1 * normalized_g1,
        boundary_residual,
    })
}

fn sampled_boundary_residual(
    patches: &[GregoryPatch],
    boundaries: &[HexBoundaryData],
) -> Option<f64> {
    if patches.len() != boundaries.len() {
        return None;
    }
    let mut worst = 0.0_f64;
    for (patch, boundary) in patches.iter().zip(boundaries.iter()) {
        for sample in 0..=18 {
            let t = (sample as f64) / 18.0;
            let row_target = eval_bezier9(&boundary.row0, t);
            let row_gap = patch.point_at(0.0, t).distance_to(&row_target);
            worst = worst.max(row_gap);
            let col_target = eval_bezier9(&boundary.col9, t);
            let col_gap = patch.point_at(t, 1.0).distance_to(&col_target);
            worst = worst.max(col_gap);
        }
    }
    Some(worst)
}

pub(super) fn boundary_residual_for_binonic(
    patches: &[RectangularGregoryBiNonic],
    boundaries: &[HexBoundaryData],
) -> Option<f64> {
    let patches = patches
        .iter()
        .cloned()
        .map(|patch| GregoryPatch::RectangularBiNonic(Box::new(patch)))
        .collect::<Vec<_>>();
    sampled_boundary_residual(&patches, boundaries)
}

fn eval_bezier9(control_points: &[Point3; 10], t: f64) -> Point3 {
    let mut work = *control_points;
    for level in 1..10 {
        for idx in 0..(10 - level) {
            work[idx] = lerp_point(work[idx], work[idx + 1], t);
        }
    }
    work[0]
}

fn recover_shared_spoke_rows(
    candidate: &mut [RectangularGregoryBiNonic],
    baseline: &[RectangularGregoryBiNonic],
    recovery: f64,
) {
    for (candidate_patch, baseline_patch) in candidate.iter_mut().zip(baseline.iter()) {
        for j in 1..9 {
            candidate_patch.grid[8][j] = lerp_point(
                candidate_patch.grid[8][j],
                baseline_patch.grid[8][j],
                recovery,
            );
            candidate_patch.grid[7][j] = lerp_point(
                candidate_patch.grid[7][j],
                baseline_patch.grid[7][j],
                recovery,
            );
        }
        for i in 1..9 {
            candidate_patch.grid[i][1] = lerp_point(
                candidate_patch.grid[i][1],
                baseline_patch.grid[i][1],
                recovery,
            );
            candidate_patch.grid[i][2] = lerp_point(
                candidate_patch.grid[i][2],
                baseline_patch.grid[i][2],
                recovery,
            );
        }
    }
}

fn fill_interiors_from_relocked_frame(patches: &mut [RectangularGregoryBiNonic]) {
    for patch in patches {
        for i in 3..=6 {
            let alpha = ((i - 2) as f64) / 5.0;
            for j in 3..=6 {
                let beta = ((j - 2) as f64) / 5.0;
                let p2j = patch.grid[2][j];
                let p7j = patch.grid[7][j];
                let pi2 = patch.grid[i][2];
                let pi7 = patch.grid[i][7];
                let p22 = patch.grid[2][2];
                let p72 = patch.grid[7][2];
                let p27 = patch.grid[2][7];
                let p77 = patch.grid[7][7];
                patch.grid[i][j] = Point3::new(
                    coons_coord(
                        alpha, beta, p2j.x, p7j.x, pi2.x, pi7.x, p22.x, p72.x, p27.x, p77.x,
                    ),
                    coons_coord(
                        alpha, beta, p2j.y, p7j.y, pi2.y, pi7.y, p22.y, p72.y, p27.y, p77.y,
                    ),
                    coons_coord(
                        alpha, beta, p2j.z, p7j.z, pi2.z, pi7.z, p22.z, p72.z, p27.z, p77.z,
                    ),
                );
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn coons_coord(
    alpha: f64,
    beta: f64,
    c2j: f64,
    c7j: f64,
    ci2: f64,
    ci7: f64,
    c22: f64,
    c72: f64,
    c27: f64,
    c77: f64,
) -> f64 {
    (1.0 - alpha) * c2j + alpha * c7j + (1.0 - beta) * ci2 + beta * ci7
        - (1.0 - alpha) * (1.0 - beta) * c22
        - alpha * (1.0 - beta) * c72
        - (1.0 - alpha) * beta * c27
        - alpha * beta * c77
}

fn lerp_point(a: Point3, b: Point3, alpha: f64) -> Point3 {
    Point3::new(
        a.x + alpha * (b.x - a.x),
        a.y + alpha * (b.y - a.y),
        a.z + alpha * (b.z - a.z),
    )
}

fn point_delta(a: Point3, b: Point3) -> [f64; 3] {
    [a.x - b.x, a.y - b.y, a.z - b.z]
}

fn first_delta(a: Point3, b: Point3) -> [f64; 3] {
    point_delta(a, b)
}

fn second_delta(a: Point3, b: Point3, c: Point3) -> [f64; 3] {
    [
        a.x - 2.0 * b.x + c.x,
        a.y - 2.0 * b.y + c.y,
        a.z - 2.0 * b.z + c.z,
    ]
}

fn scale_delta(delta: [f64; 3], scale: f64) -> [f64; 3] {
    [delta[0] * scale, delta[1] * scale, delta[2] * scale]
}

fn add_delta(point: Point3, delta: [f64; 3]) -> Point3 {
    Point3::new(point.x + delta[0], point.y + delta[1], point.z + delta[2])
}

fn subtract_delta(point: Point3, delta: [f64; 3]) -> Point3 {
    Point3::new(point.x - delta[0], point.y - delta[1], point.z - delta[2])
}
