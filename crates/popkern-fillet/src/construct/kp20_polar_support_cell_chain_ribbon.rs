//! Cross-tangent ribbon projection along chain curves.
//!
//! After the chain-curve network authors exact shared boundary
//! rows (G0), this pass makes the faces meet with a shared
//! tangent plane along every chain (G1): the pinned node normals
//! are interpolated into a normal field along the chain, and
//! each member row's row-1 controls are re-solved with minimal
//! motion so the cross-boundary derivative is perpendicular to
//! the field at every sample. The constraint is sign-agnostic
//! (perpendicularity), so no orientation bookkeeping is needed;
//! both faces' tangent planes coincide pointwise because they
//! share the curve tangent AND the plane normal.
//!
//! Row-1 endpoint controls (j ∈ {0, 9}) are frozen: they belong
//! to the transverse boundaries' row 0 (chain curves), and the
//! tangent-plane constraints applied to those curves at the
//! shared node already make them G1-compatible.

use std::collections::HashMap;

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::BoundarySide;
use super::chain_author::boundary_controls;
use super::chain_discovery::SubChain;
use super::chain_fit::ChainRowKey;
use super::chain_solver::{dbernstein9, solve_dense};
use super::row_fit::bernstein9;

/// Perpendicularity samples per row, strictly inside the row.
const RIBBON_SAMPLES_PER_ROW: usize = 12;

/// Tikhonov weight keeping the delta-field fit anchored to the
/// old controls where the samples leave slack.
const RIBBON_FIT_REGULARIZATION: f64 = 1.0e-9;

/// Weight of the second-pass perpendicularity terms against the
/// stay-near-rotation objective (both squared lengths, so the
/// weight is dimensionless).
const RIBBON_REFINE_WEIGHT: f64 = 1.0e6;

/// Parameter tolerance for node-field bracketing.
const FIELD_PARAM_TOL: f64 = 1.0e-9;

/// Per-sample frame of one row's projection: Bernstein weights,
/// the shared curve tangent, the in-plane transverse direction
/// `Ŵ = N_eff × T̂`, and the old delta's decomposition.
struct SampleFrame {
    basis: [f64; 10],
    unit_tangent: Vec3,
    binormal: Vec3,
    along: f64,
    magnitude: f64,
}

/// Outcome counts of one sub-chain's ribbon projection — kept
/// visible so silently-damped rows cannot masquerade as
/// projected ones.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct RibbonProjectionOutcome {
    pub(super) projected_row_count: usize,
    pub(super) damped_row_count: usize,
    pub(super) skipped_row_count: usize,
}

/// Project every member row's row-1 controls of one sub-chain
/// onto the shared tangent-plane field.
pub(super) fn project_sub_chain_rows_to_shared_tangent_planes(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub: &SubChain<ChainRowKey, VertexHandle>,
    node_normals: &HashMap<VertexHandle, Vec3>,
) -> Result<RibbonProjectionOutcome, &'static str> {
    let field = normal_field(sub, node_normals)?;
    let mut outcome = RibbonProjectionOutcome::default();
    for &(key, embedding) in &sub.rows {
        let patch = patches.get_mut(&key.face).ok_or("ribbon_patch_missing")?;
        match project_row(patch, key.side, |t| normal_at(&field, embedding.map(t)))? {
            RowProjection::Full => outcome.projected_row_count += 1,
            RowProjection::Damped => outcome.damped_row_count += 1,
            RowProjection::Skipped => outcome.skipped_row_count += 1,
        }
    }
    Ok(outcome)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RowProjection {
    Full,
    Damped,
    Skipped,
}

/// Sorted, sequentially sign-aligned `(chain param, normal)`
/// samples for interpolation along the chain.
fn normal_field(
    sub: &SubChain<ChainRowKey, VertexHandle>,
    node_normals: &HashMap<VertexHandle, Vec3>,
) -> Result<Vec<(f64, Vec3)>, &'static str> {
    let mut field = Vec::with_capacity(sub.nodes.len());
    for &(vertex, param) in &sub.nodes {
        if let Some(&normal) = node_normals.get(&vertex) {
            field.push((param, normal));
        }
    }
    field.sort_by(|a, b| a.0.total_cmp(&b.0));
    field.dedup_by(|a, b| (a.0 - b.0).abs() <= FIELD_PARAM_TOL);
    if field.is_empty() {
        return Err("ribbon_no_node_normals");
    }
    for index in 1..field.len() {
        if field[index - 1].1.dot(&field[index].1) < 0.0 {
            field[index].1 = field[index].1 * -1.0;
        }
    }
    Ok(field)
}

/// Clamped piecewise-linear interpolation of the normal field,
/// renormalized at the sample.
fn normal_at(field: &[(f64, Vec3)], param: f64) -> Vec3 {
    let first = field[0];
    let last = field[field.len() - 1];
    if param <= first.0 {
        return first.1;
    }
    if param >= last.0 {
        return last.1;
    }
    for window in field.windows(2) {
        let (p0, n0) = window[0];
        let (p1, n1) = window[1];
        if param <= p1 + FIELD_PARAM_TOL {
            let span = p1 - p0;
            if span <= FIELD_PARAM_TOL {
                return n0;
            }
            let fraction = (param - p0) / span;
            let blended = n0 * (1.0 - fraction) + n1 * fraction;
            return blended.normalized().unwrap_or(n0);
        }
    }
    last.1
}

/// Re-solve row-1 controls `j ∈ 1..=8` so the cross-boundary
/// derivative field lies in the shared tangent plane: at every
/// interior sample the old cross-derivative is decomposed
/// against the frame `(T̂, Ŵ = N_eff × T̂, N_eff)` and its
/// off-tangent component is rotated rigidly about the curve
/// tangent into the plane (`b' = sign(b)·√(b² + c²)`) —
/// length-preserving, so the derivative never collapses and the
/// curvature measurement never degenerates. The free controls
/// are then least-squares fit to the rotated delta field
/// (axis-separable; frozen endpoint contributions moved to the
/// right-hand side).
fn project_row(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    normal_for_t: impl Fn(f64) -> Vec3,
) -> Result<RowProjection, &'static str> {
    let row0 = boundary_controls(patch, side);
    let row1 = second_row_controls(patch, side);
    // First pass: per-sample frames and decompositions. The
    // rotated in-plane component takes ONE consistent sign for
    // the whole row (energy-weighted majority of the old
    // in-plane components): near-orthogonal sheets have b ≈ 0
    // with large c, where a per-sample sign flips arbitrarily
    // and makes the rotated field discontinuous — the fit then
    // oscillates and the collapse guard rejects the row.
    let mut frames = Vec::with_capacity(RIBBON_SAMPLES_PER_ROW);
    let mut sign_vote = 0.0_f64;
    for sample in 0..RIBBON_SAMPLES_PER_ROW {
        let t = (sample as f64 + 1.0) / (RIBBON_SAMPLES_PER_ROW as f64 + 1.0);
        let basis = bernstein9(t);
        // G1 needs the two faces' tangent planes to COINCIDE,
        // and both planes contain the shared curve tangent, so
        // the field normal is orthogonalized against the curve
        // tangent pointwise: the shared plane is span(T̂, Ŵ) on
        // both sides. (At nodes the chain curves are already
        // tangent-plane-constrained, so the correction vanishes
        // there.)
        let Some(unit_tangent) = row_tangent(&row0, t).normalized() else {
            continue;
        };
        let raw_normal = normal_for_t(t);
        let Some(normal) = (raw_normal - unit_tangent * raw_normal.dot(&unit_tangent)).normalized()
        else {
            continue;
        };
        let binormal = normal.cross(&unit_tangent);
        let mut delta = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..10 {
            let difference = row1[i] - row0[i];
            delta = delta + difference * basis[i];
        }
        let along = delta.dot(&unit_tangent);
        let in_plane = delta.dot(&binormal);
        let out_of_plane = delta.dot(&normal);
        let magnitude = (in_plane * in_plane + out_of_plane * out_of_plane).sqrt();
        sign_vote += in_plane * magnitude;
        frames.push(SampleFrame {
            basis,
            unit_tangent,
            binormal,
            along,
            magnitude,
        });
    }
    if frames.is_empty() {
        return Ok(RowProjection::Skipped);
    }
    let row_sign = if sign_vote >= 0.0 { 1.0 } else { -1.0 };
    let mut gram = vec![vec![0.0_f64; 8]; 8];
    let mut rhs = [[0.0_f64; 8], [0.0_f64; 8], [0.0_f64; 8]];
    for frame in &frames {
        let target =
            frame.unit_tangent * frame.along + frame.binormal * (row_sign * frame.magnitude);
        // Move the frozen endpoint deltas across so the fit sees
        // only the free controls' share of the field.
        let frozen = (row1[0] - row0[0]) * frame.basis[0] + (row1[9] - row0[9]) * frame.basis[9];
        let residual = target - frozen;
        for i in 0..8 {
            for j in 0..8 {
                gram[i][j] += frame.basis[i + 1] * frame.basis[j + 1];
            }
            rhs[0][i] += frame.basis[i + 1] * residual.x;
            rhs[1][i] += frame.basis[i + 1] * residual.y;
            rhs[2][i] += frame.basis[i + 1] * residual.z;
        }
    }
    // Prior: stay near the old free deltas where the samples
    // leave slack.
    for i in 0..8 {
        gram[i][i] += RIBBON_FIT_REGULARIZATION;
        let old_delta = row1[i + 1] - row0[i + 1];
        rhs[0][i] += RIBBON_FIT_REGULARIZATION * old_delta.x;
        rhs[1][i] += RIBBON_FIT_REGULARIZATION * old_delta.y;
        rhs[2][i] += RIBBON_FIT_REGULARIZATION * old_delta.z;
    }
    let xs = solve_dense(gram.clone(), rhs[0].to_vec()).ok_or("ribbon_projection_singular")?;
    let ys = solve_dense(gram.clone(), rhs[1].to_vec()).ok_or("ribbon_projection_singular")?;
    let zs = solve_dense(gram, rhs[2].to_vec()).ok_or("ribbon_projection_singular")?;
    let mut rotated = row1;
    for i in 0..8 {
        rotated[i + 1] = Point3::new(
            row0[i + 1].x + xs[i],
            row0[i + 1].y + ys[i],
            row0[i + 1].z + zs[i],
        );
    }
    // The rotated target field is not polynomial (its magnitude
    // carries a square root), so where the old sheet geometry is
    // wild the 8-control fit leaves an out-of-plane residual of
    // the field's own scale. Second pass: starting FROM the
    // rotated solution, drive the residual at the samples to
    // zero with heavily weighted perpendicularity terms — the
    // correction is small now, so it cannot collapse the
    // derivative the way a single-pass perpendicularity solve
    // can.
    let refined = perpendicularity_refine(&row0, &rotated, &frames);
    // Accept the strongest candidate whose derivative field
    // stays non-collapsed everywhere; an unchanged row is always
    // acceptable (counted as skipped, never silent).
    let mut candidates: Vec<([Point3; 10], RowProjection)> = Vec::new();
    if let Some(refined) = refined {
        candidates.push((refined, RowProjection::Full));
    }
    candidates.push((rotated, RowProjection::Full));
    for blend in [0.5, 0.25] {
        let mut damped = row1;
        for index in 1..9 {
            damped[index] = row1[index] + (rotated[index] - row1[index]) * blend;
        }
        candidates.push((damped, RowProjection::Damped));
    }
    for (candidate, projection) in candidates {
        if cross_derivative_stays_resolved(&row0, &candidate) {
            for (index, point) in candidate.iter().enumerate().take(9).skip(1) {
                set_second_row_control(patch, side, index, *point);
            }
            return Ok(projection);
        }
    }
    Ok(RowProjection::Skipped)
}

/// Second-pass refinement: minimize motion from the rotated
/// solution subject to heavily weighted perpendicularity of the
/// cross-derivative to the (tangent-orthogonalized) plane normal
/// at every sample. Solved per the 24 free coordinates with
/// normal equations `(I + W·CᵀC) q = q_prior + W·Cᵀd`.
fn perpendicularity_refine(
    row0: &[Point3; 10],
    rotated: &[Point3; 10],
    frames: &[SampleFrame],
) -> Option<[Point3; 10]> {
    let unknowns = 24;
    let mut matrix = vec![vec![0.0_f64; unknowns]; unknowns];
    let mut rhs = vec![0.0_f64; unknowns];
    for i in 0..8 {
        let prior = rotated[i + 1] - row0[i + 1];
        matrix[i][i] = 1.0;
        matrix[8 + i][8 + i] = 1.0;
        matrix[16 + i][16 + i] = 1.0;
        rhs[i] = prior.x;
        rhs[8 + i] = prior.y;
        rhs[16 + i] = prior.z;
    }
    for frame in frames {
        // The plane normal of this sample's shared tangent plane.
        let normal = frame.unit_tangent.cross(&frame.binormal);
        let mut coefficients = [0.0_f64; 24];
        for i in 0..8 {
            coefficients[i] = frame.basis[i + 1] * normal.x;
            coefficients[8 + i] = frame.basis[i + 1] * normal.y;
            coefficients[16 + i] = frame.basis[i + 1] * normal.z;
        }
        // Frozen endpoint deltas keep their share of the field.
        let frozen =
            (rotated[0] - row0[0]) * frame.basis[0] + (rotated[9] - row0[9]) * frame.basis[9];
        let datum = -frozen.dot(&normal);
        for i in 0..unknowns {
            for j in 0..unknowns {
                matrix[i][j] += RIBBON_REFINE_WEIGHT * coefficients[i] * coefficients[j];
            }
            rhs[i] += RIBBON_REFINE_WEIGHT * coefficients[i] * datum;
        }
    }
    let solution = solve_dense(matrix, rhs)?;
    let mut out = *rotated;
    for i in 0..8 {
        out[i + 1] = Point3::new(
            row0[i + 1].x + solution[i],
            row0[i + 1].y + solution[8 + i],
            row0[i + 1].z + solution[16 + i],
        );
    }
    Some(out)
}

/// Whether the row-1 candidate keeps a resolved (non-collapsed)
/// cross-derivative along the whole row, relative to its own
/// largest magnitude.
fn cross_derivative_stays_resolved(row0: &[Point3; 10], row1: &[Point3; 10]) -> bool {
    let mut min_magnitude = f64::MAX;
    let mut max_magnitude = 0.0_f64;
    for sample in 0..=32 {
        let t = sample as f64 / 32.0;
        let basis = bernstein9(t);
        let mut delta = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..10 {
            delta = delta + (row1[i] - row0[i]) * basis[i];
        }
        let magnitude = delta.dot(&delta).sqrt();
        min_magnitude = min_magnitude.min(magnitude);
        max_magnitude = max_magnitude.max(magnitude);
    }
    min_magnitude > 1.0e-4 * max_magnitude
}

/// Chain-curve tangent at `t`, from the shared row-0 controls.
fn row_tangent(row0: &[Point3; 10], t: f64) -> Vec3 {
    let dbasis = dbernstein9(t);
    let mut tangent = Vec3::new(0.0, 0.0, 0.0);
    for (index, weight) in dbasis.iter().enumerate() {
        tangent = Vec3::new(
            tangent.x + row0[index].x * weight,
            tangent.y + row0[index].y * weight,
            tangent.z + row0[index].z * weight,
        );
    }
    tangent
}

fn second_row_controls(patch: &RectangularGregoryBiNonic, side: BoundarySide) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = match side {
            BoundarySide::UMin => patch.grid[1][index],
            BoundarySide::UMax => patch.grid[8][index],
            BoundarySide::VMin => patch.grid[index][1],
            BoundarySide::VMax => patch.grid[index][8],
        };
    }
    out
}

fn set_second_row_control(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    index: usize,
    point: Point3,
) {
    match side {
        BoundarySide::UMin => patch.grid[1][index] = point,
        BoundarySide::UMax => patch.grid[8][index] = point,
        BoundarySide::VMin => patch.grid[index][1] = point,
        BoundarySide::VMax => patch.grid[index][8] = point,
    }
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_chain_ribbon_tests.rs"]
mod tests;
