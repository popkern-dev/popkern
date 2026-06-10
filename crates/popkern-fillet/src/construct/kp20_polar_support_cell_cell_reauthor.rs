//! Constructive support-cell re-authoring from the chain network.
//!
//! Instead of projecting the old (mutually disagreeing) sheet
//! geometry toward shared planes, every cell's near-boundary
//! band is DERIVED from shared objects: per sub-chain, the
//! canonical curve `C` and a degree-8 transverse ribbon field
//! `R`. Each member row's cross-derivative becomes
//! `D(t) = α(t)·C'(t) + β(t)·R(t)` with LINEAR `α, β` solved
//! from the exact end legs (the transverse chain curves' first
//! control legs) — a degree-9 polynomial, exactly representable
//! in row 1, lying in `span(C'(s), R(s))` pointwise. Both sides
//! of every trim therefore share a tangent plane by
//! construction, with no fit residual. Corner solves are exact
//! because the legs, `C'`, and `R` all lie in the pinned node
//! plane (the chain curves are tangent-plane-constrained there).
//!
//! Twist controls (a grid corner's diagonal neighbour) are
//! demanded by two boundaries; the writes are averaged and
//! counted — the one approximation of the construction. Cell
//! interiors are refilled by a discrete Coons blend of the four
//! derived bands.

use std::collections::HashMap;

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::BoundarySide;
use super::chain_fit::{ChainRowKey, boundary_controls};
use super::chain_ribbon_field::{
    eval_ribbon_field, fit_ribbon_field, ribbon_field_segment, row_derivative_coefficients,
};

/// Sample count for the per-sub-chain ribbon-field fit.
const FIELD_FIT_SAMPLES: usize = 17;

/// Determinant floor for the 2x2 end-leg solves (relative to the
/// basis vectors' scale).
const LEG_SOLVE_DET_TOL: f64 = 1.0e-10;

/// One authored sub-chain, retained by the curve-authoring pass
/// for the constructive re-author: the canonical curve, every
/// member row's local interval, and the node params/identities
/// (all in the sub-chain's local `[0, 1]` domain).
#[derive(Debug, Clone)]
pub(super) struct AuthoredSubChain {
    pub(super) chain_cps: [Point3; 10],
    pub(super) rows: Vec<(ChainRowKey, (f64, f64))>,
    pub(super) nodes: Vec<(VertexHandle, f64)>,
}

/// Outcome counts of the constructive re-author.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct CellReauthorOutcome {
    pub(super) derived_row_count: usize,
    pub(super) leg_solve_declined_count: usize,
    pub(super) field_declined_sub_chain_count: usize,
    pub(super) twist_conflict_count: usize,
    pub(super) refilled_cell_count: usize,
}

/// Derive every chained row's row-1 band from the shared
/// `(C, R)` fields, then refill the touched cells' interiors.
pub(super) fn reauthor_cells_from_chain_network(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub_chains: &[AuthoredSubChain],
    node_normals: &HashMap<VertexHandle, Vec3>,
) -> CellReauthorOutcome {
    let mut outcome = CellReauthorOutcome::default();
    // Collect every derived control write; twist controls are
    // demanded by two boundaries of the same face and averaged.
    let mut writes: HashMap<(FaceHandle, usize, usize), (Point3, usize)> = HashMap::new();
    for sub in sub_chains {
        let Some(field) = author_ribbon_field(sub, node_normals) else {
            outcome.field_declined_sub_chain_count += 1;
            continue;
        };
        for &(key, (start, end)) in &sub.rows {
            let Some(patch) = patches.get(&key.face) else {
                continue;
            };
            let row0 = boundary_controls(patch, key.side);
            let row1 = second_row_controls(patch, key.side);
            let field_segment = ribbon_field_segment(&field, start, end);
            match derive_row1(&row0, &row1, &field_segment) {
                Some(derived) => {
                    for (k, point) in derived.iter().enumerate().take(9).skip(1) {
                        let (i, j) = second_row_grid_index(key.side, k);
                        writes
                            .entry((key.face, i, j))
                            .and_modify(|(sum, count)| {
                                *sum =
                                    Point3::new(sum.x + point.x, sum.y + point.y, sum.z + point.z);
                                *count += 1;
                            })
                            .or_insert((*point, 1));
                    }
                    outcome.derived_row_count += 1;
                }
                None => outcome.leg_solve_declined_count += 1,
            }
        }
    }
    let mut touched_faces = Vec::new();
    for (&(face, i, j), &(sum, count)) in &writes {
        if count > 1 {
            outcome.twist_conflict_count += 1;
        }
        let inverse = 1.0 / count as f64;
        if let Some(patch) = patches.get_mut(&face) {
            patch.grid[i][j] = Point3::new(sum.x * inverse, sum.y * inverse, sum.z * inverse);
            if !touched_faces.contains(&face) {
                touched_faces.push(face);
            }
        }
    }
    for face in touched_faces {
        if let Some(patch) = patches.get_mut(&face) {
            refill_interior_coons(patch);
            outcome.refilled_cell_count += 1;
        }
    }
    outcome
}

/// Author the sub-chain's transverse ribbon field: direction
/// targets `W(s) = N(s) × T̂(s)` (the in-plane transverse of the
/// interpolated node-normal field), interpolated exactly at the
/// nodes so the corner solves close.
fn author_ribbon_field(
    sub: &AuthoredSubChain,
    node_normals: &HashMap<VertexHandle, Vec3>,
) -> Option<[Vec3; 9]> {
    let mut normal_field: Vec<(f64, Vec3)> = Vec::with_capacity(sub.nodes.len());
    for &(vertex, param) in &sub.nodes {
        if let Some(&normal) = node_normals.get(&vertex) {
            normal_field.push((param, normal));
        }
    }
    normal_field.sort_by(|a, b| a.0.total_cmp(&b.0));
    normal_field.dedup_by(|a, b| (a.0 - b.0).abs() <= 1.0e-9);
    if normal_field.is_empty() {
        return None;
    }
    for index in 1..normal_field.len() {
        if normal_field[index - 1].1.dot(&normal_field[index].1) < 0.0 {
            normal_field[index].1 = normal_field[index].1 * -1.0;
        }
    }
    let derivative = row_derivative_coefficients(&sub.chain_cps);
    let transverse_at = |t: f64| -> Option<Vec3> {
        let tangent = eval_ribbon_field(&derivative, t).normalized()?;
        let normal = interpolate_normal(&normal_field, t);
        (normal.cross(&tangent)).normalized()
    };
    let mut samples = Vec::with_capacity(FIELD_FIT_SAMPLES);
    let mut previous: Option<Vec3> = None;
    for index in 0..FIELD_FIT_SAMPLES {
        let t = index as f64 / (FIELD_FIT_SAMPLES - 1) as f64;
        let Some(mut direction) = transverse_at(t) else {
            continue;
        };
        if let Some(previous) = previous
            && previous.dot(&direction) < 0.0
        {
            direction = direction * -1.0;
        }
        previous = Some(direction);
        samples.push((t, direction));
    }
    let mut node_constraints = Vec::with_capacity(normal_field.len());
    for &(param, _) in &normal_field {
        let Some(mut direction) = transverse_at(param) else {
            continue;
        };
        // Align with the fitted sample stream's orientation.
        if let Some((_, nearest)) = samples
            .iter()
            .min_by(|a, b| (a.0 - param).abs().total_cmp(&(b.0 - param).abs()))
            && nearest.dot(&direction) < 0.0
        {
            direction = direction * -1.0;
        }
        node_constraints.push((param, direction));
    }
    fit_ribbon_field(&samples, &node_constraints).ok()
}

/// Piecewise-linear interpolation of the sign-aligned node
/// normals, renormalized.
fn interpolate_normal(field: &[(f64, Vec3)], param: f64) -> Vec3 {
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
        if param <= p1 + 1.0e-9 {
            let span = p1 - p0;
            if span <= 1.0e-9 {
                return n0;
            }
            let fraction = (param - p0) / span;
            let blended = n0 * (1.0 - fraction) + n1 * fraction;
            return blended.normalized().unwrap_or(n0);
        }
    }
    last.1
}

/// Derive the full row-1 band from `(C', R)` and the exact end
/// legs: `D(t) = α(t)·C'(t) + β(t)·R(t)` with linear `α, β`
/// matched to `D(0) = E₀`, `D(1) = E₁`, expanded into degree-9
/// Bezier coefficients, then `row1ₖ = row0ₖ + Dₖ / 9`.
fn derive_row1(
    row0: &[Point3; 10],
    row1: &[Point3; 10],
    field_segment: &[Vec3; 9],
) -> Option<[Point3; 10]> {
    let derivative = row_derivative_coefficients(row0);
    let leg_start = (row1[0] - row0[0]) * 9.0;
    let leg_end = (row1[9] - row0[9]) * 9.0;
    let (alpha_start, beta_start) = solve_leg(leg_start, derivative[0], field_segment[0])?;
    let (alpha_end, beta_end) = solve_leg(leg_end, derivative[8], field_segment[8])?;
    let mut out = *row1;
    for k in 1..9 {
        // Product of a linear blend with a degree-8 field:
        // coeffₖ = (9−k)/9·(α₀ d_k + β₀ r_k)
        //        +    k/9·(α₁ d_{k−1} + β₁ r_{k−1}).
        let head = derivative[k] * alpha_start + field_segment[k] * beta_start;
        let tail = derivative[k - 1] * alpha_end + field_segment[k - 1] * beta_end;
        let d_k = head * ((9 - k) as f64 / 9.0) + tail * (k as f64 / 9.0);
        out[k] = row0[k] + d_k * (1.0 / 9.0);
    }
    Some(out)
}

/// Solve `leg = α·u + β·v` in the least-squares sense (exact
/// when the leg is coplanar with the basis, which the node
/// tangent-plane constraints guarantee).
fn solve_leg(leg: Vec3, u: Vec3, v: Vec3) -> Option<(f64, f64)> {
    let uu = u.dot(&u);
    let uv = u.dot(&v);
    let vv = v.dot(&v);
    let det = uu * vv - uv * uv;
    if det.abs() <= LEG_SOLVE_DET_TOL * uu.max(vv).max(1.0e-300) {
        return None;
    }
    let bu = u.dot(&leg);
    let bv = v.dot(&leg);
    Some(((vv * bu - uv * bv) / det, (uu * bv - uv * bu) / det))
}

/// Discrete Coons fill of the interior controls (indices 2..=7)
/// from the four derived band rows/columns at indices 1 and 8.
fn refill_interior_coons(patch: &mut RectangularGregoryBiNonic) {
    let grid = patch.grid;
    let weighted = |terms: &[(f64, Point3)]| -> Point3 {
        let mut out = Point3::new(0.0, 0.0, 0.0);
        for &(weight, point) in terms {
            out = Point3::new(
                out.x + weight * point.x,
                out.y + weight * point.y,
                out.z + weight * point.z,
            );
        }
        out
    };
    for i in 2..8 {
        let a = (i as f64 - 1.0) / 7.0;
        for j in 2..8 {
            let b = (j as f64 - 1.0) / 7.0;
            patch.grid[i][j] = weighted(&[
                (1.0 - a, grid[1][j]),
                (a, grid[8][j]),
                (1.0 - b, grid[i][1]),
                (b, grid[i][8]),
                (-(1.0 - a) * (1.0 - b), grid[1][1]),
                (-(1.0 - a) * b, grid[1][8]),
                (-a * (1.0 - b), grid[8][1]),
                (-a * b, grid[8][8]),
            ]);
        }
    }
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

fn second_row_grid_index(side: BoundarySide, k: usize) -> (usize, usize) {
    match side {
        BoundarySide::UMin => (1, k),
        BoundarySide::UMax => (8, k),
        BoundarySide::VMin => (k, 1),
        BoundarySide::VMax => (k, 8),
    }
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_cell_reauthor_tests.rs"]
mod tests;
