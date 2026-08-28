//! Cap-vertex-aware support-cell split parameter computation.
//!
//! Stage A3 of `docs/decisions/2026-06-21-kp20-cap-vertex-aware-
//! support-cell-split.md`. Given a macro piece and the cap-outer
//! vertices that lie on one of its four boundary sides, compute the
//! split parameters that place each sub-cell boundary at the
//! projected vertex parameters — replacing the uniform `i / cells`
//! split that misaligns from the cap-outer topology on variable-
//! radius fixtures.
//!
//! Wired into [`super::split_uniform_u`] / `split_uniform_v` and the
//! macro-piece subdivision in Stage B2; A3 ships the algorithm,
//! diagnostic data, and tests against the algorithm's return value.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::Kp20Bi6Patch;

use super::support_cell_boundary_projection::{
    ProjectionOutcome, SupportCellPatchSide, project_to_boundary,
};
use crate::types::{
    SupportCellSplitAxisDiagnostic, SupportCellSplitAxisStatus, SupportCellSplitDiagnostic,
    SupportCellSplitFallbackReason,
};

/// Cap-outer vertices known to lie on a single boundary side of a
/// macro piece. The caller resolves which `SupportCellPatchSide` they
/// belong to (a Stage B1 helper does this from BRep adjacency) and
/// passes the 3D positions in cap-outer half-edge chain order.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)] // first caller lands in Stage B2
pub(crate) struct CapVertexAxisProjection<'a> {
    pub side: SupportCellPatchSide,
    pub vertices: &'a [Point3],
}

/// Outcome of the cap-vertex-aware split along one axis: the
/// parameter sequence (`cells_axis - 1` entries) plus the diagnostic
/// facts that go into [`SupportCellSplitDiagnostic`].
#[derive(Debug, Clone)]
#[allow(dead_code)] // first caller lands in Stage B2
pub(crate) struct CapVertexSplitOutcome {
    pub params: Vec<f64>,
    pub diagnostic: SupportCellSplitAxisDiagnostic,
}

/// Parametric dedup / endpoint-deadband tolerance (matches the
/// projection helper's `DEDUP_TOL`).
const DEDUP_TOL: f64 = 1.0e-9;

/// Maximum `|projected_t - i / cells|` below which the projected
/// sequence routes through the uniform single-step formula.
const UNIFORM_EQUIVALENCE_TOL: f64 = 1.0e-12;

/// Minimum sub-cell arc-length factor (multiplied by the macro piece's
/// bounding-box diagonal) below which the sequence is rejected as
/// inadmissible — protects the boundary fit from numerically
/// pathological thin sub-cells.
const MIN_SUBCELL_ARC_FACTOR: f64 = 0.05;

/// Boundary samples per sub-cell for arc-length estimation
/// (6 points = 5 chord segments per sub-cell).
const ARC_SAMPLES_PER_SUBCELL: usize = 6;

/// Compute the cap-vertex-aware split parameters along one axis of a
/// macro piece. Returns the parameter sequence plus the diagnostic
/// facts for that axis. When `cap_projection` is `None`, returns the
/// uniform `i / cells_axis` sequence and marks the axis status as
/// `UniformNoProjection`.
#[allow(dead_code)] // first caller lands in Stage B2
pub(crate) fn cap_vertex_splits_along_axis(
    base: &Kp20Bi6Patch,
    cap_projection: Option<&CapVertexAxisProjection>,
    cells_axis: usize,
) -> CapVertexSplitOutcome {
    let Some(cap) = cap_projection else {
        return uniform_outcome(
            cells_axis,
            SupportCellSplitAxisStatus::UniformNoProjection,
            None,
            0,
        );
    };
    if cap.vertices.is_empty() {
        return uniform_outcome(
            cells_axis,
            SupportCellSplitAxisStatus::UniformFallback,
            Some(SupportCellSplitFallbackReason::NoProjection),
            0,
        );
    }
    // Project + filter.
    let (projected, divergence) = project_vertices(base, cap);
    if divergence {
        return uniform_outcome(
            cells_axis,
            SupportCellSplitAxisStatus::UniformFallback,
            Some(SupportCellSplitFallbackReason::DidNotConverge),
            cap.vertices.len(),
        );
    }
    let projected_count = projected.len();
    // Match sequence length to cells_axis - 1.
    let required = cells_axis.saturating_sub(1);
    if projected.len() > required {
        return CapVertexSplitOutcome {
            params: uniform_params(cells_axis),
            diagnostic: SupportCellSplitAxisDiagnostic {
                status: SupportCellSplitAxisStatus::UniformFallback,
                fallback_reason: Some(SupportCellSplitFallbackReason::Overcount),
                cap_vertex_count: cap.vertices.len(),
                projected_count,
                uniform_fill_count: 0,
                min_subcell_arc: 0.0,
                uniform_equivalence_max_delta: 0.0,
            },
        };
    }
    let mut matched = projected.clone();
    let uniform_fill_count = required - projected.len();
    if uniform_fill_count > 0 {
        fill_largest_gaps(&mut matched, uniform_fill_count);
    }
    // Uniform-equivalence shortcut.
    let equivalence_delta = uniform_equivalence_max_delta(&matched, cells_axis);
    if equivalence_delta < UNIFORM_EQUIVALENCE_TOL {
        return CapVertexSplitOutcome {
            params: uniform_params(cells_axis),
            diagnostic: SupportCellSplitAxisDiagnostic {
                status: SupportCellSplitAxisStatus::UniformEquivalent,
                fallback_reason: None,
                cap_vertex_count: cap.vertices.len(),
                projected_count,
                uniform_fill_count,
                min_subcell_arc: 0.0,
                uniform_equivalence_max_delta: equivalence_delta,
            },
        };
    }
    // Admissibility guard.
    let (min_subcell_arc, min_required) = check_admissibility(base, cap.side, &matched);
    if min_subcell_arc < min_required {
        return CapVertexSplitOutcome {
            params: uniform_params(cells_axis),
            diagnostic: SupportCellSplitAxisDiagnostic {
                status: SupportCellSplitAxisStatus::UniformFallback,
                fallback_reason: Some(SupportCellSplitFallbackReason::MinArc),
                cap_vertex_count: cap.vertices.len(),
                projected_count,
                uniform_fill_count,
                min_subcell_arc,
                uniform_equivalence_max_delta: equivalence_delta,
            },
        };
    }
    CapVertexSplitOutcome {
        params: matched,
        diagnostic: SupportCellSplitAxisDiagnostic {
            status: SupportCellSplitAxisStatus::VertexAware,
            fallback_reason: None,
            cap_vertex_count: cap.vertices.len(),
            projected_count,
            uniform_fill_count,
            min_subcell_arc,
            uniform_equivalence_max_delta: equivalence_delta,
        },
    }
}

/// Build the per-macro-piece diagnostic from the per-axis outcomes
/// and the macro piece's identity.
#[allow(dead_code)] // first caller lands in Stage B2
pub(crate) fn support_cell_split_diagnostic(
    sector_index: usize,
    piece_index: usize,
    cells_u: usize,
    cells_v: usize,
    cap_vertex_interior_count: usize,
    u_outcome: &CapVertexSplitOutcome,
    v_outcome: &CapVertexSplitOutcome,
) -> SupportCellSplitDiagnostic {
    SupportCellSplitDiagnostic {
        sector_index,
        piece_index,
        cells_u,
        cells_v,
        cap_vertex_interior_count,
        u_axis: u_outcome.diagnostic,
        v_axis: v_outcome.diagnostic,
    }
}

fn uniform_outcome(
    cells_axis: usize,
    status: SupportCellSplitAxisStatus,
    fallback_reason: Option<SupportCellSplitFallbackReason>,
    cap_vertex_count: usize,
) -> CapVertexSplitOutcome {
    CapVertexSplitOutcome {
        params: uniform_params(cells_axis),
        diagnostic: SupportCellSplitAxisDiagnostic {
            status,
            fallback_reason,
            cap_vertex_count,
            projected_count: 0,
            uniform_fill_count: 0,
            min_subcell_arc: 0.0,
            uniform_equivalence_max_delta: 0.0,
        },
    }
}

fn uniform_params(cells_axis: usize) -> Vec<f64> {
    (1..cells_axis)
        .map(|i| i as f64 / cells_axis as f64)
        .collect()
}

fn project_vertices(base: &Kp20Bi6Patch, cap: &CapVertexAxisProjection) -> (Vec<f64>, bool) {
    let mut projected: Vec<f64> = Vec::with_capacity(cap.vertices.len());
    for vertex in cap.vertices {
        match project_to_boundary(base, cap.side, vertex) {
            ProjectionOutcome::Converged(t) => projected.push(t),
            ProjectionOutcome::OutsideDomain { .. } => {}
            ProjectionOutcome::DidNotConverge { .. } => return (Vec::new(), true),
        }
    }
    projected.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    projected.dedup_by(|a, b| (*a - *b).abs() < DEDUP_TOL);
    (projected, false)
}

fn fill_largest_gaps(params: &mut Vec<f64>, fill_count: usize) {
    for _ in 0..fill_count {
        let mut boundaries: Vec<f64> = Vec::with_capacity(params.len() + 2);
        boundaries.push(0.0);
        boundaries.extend_from_slice(params);
        boundaries.push(1.0);
        let mut largest_gap = 0.0_f64;
        let mut largest_index = 0;
        for index in 0..boundaries.len() - 1 {
            let gap = boundaries[index + 1] - boundaries[index];
            if gap > largest_gap {
                largest_gap = gap;
                largest_index = index;
            }
        }
        let inserted = 0.5 * (boundaries[largest_index] + boundaries[largest_index + 1]);
        params.push(inserted);
        params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    }
}

fn uniform_equivalence_max_delta(params: &[f64], cells_axis: usize) -> f64 {
    params
        .iter()
        .enumerate()
        .map(|(i, t)| (t - (i + 1) as f64 / cells_axis as f64).abs())
        .fold(0.0_f64, f64::max)
}

fn check_admissibility(
    base: &Kp20Bi6Patch,
    side: SupportCellPatchSide,
    params: &[f64],
) -> (f64, f64) {
    let bbox_diag = macro_piece_bbox_diagonal(base);
    let min_required = MIN_SUBCELL_ARC_FACTOR * bbox_diag;
    let mut boundaries: Vec<f64> = Vec::with_capacity(params.len() + 2);
    boundaries.push(0.0);
    boundaries.extend_from_slice(params);
    boundaries.push(1.0);
    let mut min_arc = f64::INFINITY;
    for window in boundaries.windows(2) {
        let start = window[0];
        let end = window[1];
        let arc = sub_cell_arc_length(base, side, start, end);
        if arc < min_arc {
            min_arc = arc;
        }
    }
    (min_arc, min_required)
}

fn sub_cell_arc_length(
    base: &Kp20Bi6Patch,
    side: SupportCellPatchSide,
    start: f64,
    end: f64,
) -> f64 {
    let mut last = boundary_point(base, side, start);
    let mut sum = 0.0_f64;
    for sample in 1..ARC_SAMPLES_PER_SUBCELL {
        let local = sample as f64 / (ARC_SAMPLES_PER_SUBCELL - 1) as f64;
        let t = start + (end - start) * local;
        let point = boundary_point(base, side, t);
        sum += point.distance_to(&last);
        last = point;
    }
    sum
}

fn boundary_point(base: &Kp20Bi6Patch, side: SupportCellPatchSide, t: f64) -> Point3 {
    let cps = match side {
        SupportCellPatchSide::UMin => [
            base.cps[0][0],
            base.cps[0][1],
            base.cps[0][2],
            base.cps[0][3],
            base.cps[0][4],
            base.cps[0][5],
            base.cps[0][6],
        ],
        SupportCellPatchSide::UMax => [
            base.cps[6][0],
            base.cps[6][1],
            base.cps[6][2],
            base.cps[6][3],
            base.cps[6][4],
            base.cps[6][5],
            base.cps[6][6],
        ],
        SupportCellPatchSide::VMin => [
            base.cps[0][0],
            base.cps[1][0],
            base.cps[2][0],
            base.cps[3][0],
            base.cps[4][0],
            base.cps[5][0],
            base.cps[6][0],
        ],
        SupportCellPatchSide::VMax => [
            base.cps[0][6],
            base.cps[1][6],
            base.cps[2][6],
            base.cps[3][6],
            base.cps[4][6],
            base.cps[5][6],
            base.cps[6][6],
        ],
    };
    let mut work = cps;
    for level in 1..=6 {
        for index in 0..=(6 - level) {
            let a = work[index];
            let b = work[index + 1];
            work[index] = a + (b - a) * t;
        }
    }
    work[0]
}

fn macro_piece_bbox_diagonal(base: &Kp20Bi6Patch) -> f64 {
    let mut min = base.cps[0][0];
    let mut max = base.cps[0][0];
    for row in 0..7 {
        for column in 0..7 {
            let point = base.cps[row][column];
            min = Point3::new(min.x.min(point.x), min.y.min(point.y), min.z.min(point.z));
            max = Point3::new(max.x.max(point.x), max.y.max(point.y), max.z.max(point.z));
        }
    }
    max.distance_to(&min)
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_cap_vertex_split_tests.rs"]
mod tests;
