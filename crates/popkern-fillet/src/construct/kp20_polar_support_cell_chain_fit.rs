//! Canonical chain-curve authoring for KP20 support-cell rows.
//!
//! For each discovered sub-chain, fit ONE degree-9 canonical
//! curve over the chain domain — interpolating the pinned node
//! positions exactly, least-squares against the current row
//! geometry elsewhere — then assign every member row the exact
//! de Casteljau restriction of that curve over its interval.
//! Every trim between rows of the chain becomes a control-point
//! identity: full-span, one-sided, and both-partial T-junction
//! pairs alike.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::chain_discovery::{ChainTrim, RowChainEmbedding, SubChain, discover_sub_chains};
use super::row_fit::bernstein9;
use super::row_segment::bezier_row_segment;
use super::{
    BoundarySide, full_boundary_trims, rectangular_binonic_on_face, set_rectangular_binonic_on_face,
};

/// Samples per member row feeding the least-squares part of the
/// canonical fit.
const CHAIN_FIT_SAMPLES_PER_ROW: usize = 13;
/// Maximum interpolated nodes per sub-chain: a degree-9 curve
/// has 10 control points; pinning more than 9 leaves no
/// least-squares freedom and degrades conditioning.
const MAX_INTERPOLATED_NODES: usize = 9;
/// Parameter tolerance for matching a row endpoint to a node
/// (chain-parameter units — one row spans 1.0).
const NODE_PARAM_TOL: f64 = 1.0e-6;
/// Tikhonov weight tying the fit to the sample cloud scale.
const FIT_REGULARIZATION: f64 = 1.0e-9;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ChainRowKey {
    face: FaceHandle,
    side: BoundarySide,
}

/// Author canonical curves for every discoverable chain and
/// assign all member rows as exact restrictions. The status
/// records authored/declined chain counts for the candidate
/// score fact.
pub(super) fn author_chain_curve_network(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    let trims = full_boundary_trims(brep, slots)?;
    let mut chain_trims = Vec::with_capacity(trims.len());
    for trim in &trims {
        chain_trims.push(ChainTrim {
            row_a: ChainRowKey {
                face: trim.face_a,
                side: trim.side_a,
            },
            interval_a: (trim.start_a, trim.end_a),
            row_b: ChainRowKey {
                face: trim.face_b,
                side: trim.side_b,
            },
            interval_b: (trim.start_b, trim.end_b),
            // start_a corresponds to half-edge A's origin and
            // end_a to its destination (= half-edge B's origin),
            // so vertex handles give exact node identity.
            node_start: brep.half_edges[trim.half_edge_a].origin,
            node_end: brep.half_edges[trim.half_edge_b].origin,
        });
    }
    let discovery = discover_sub_chains(&chain_trims);
    let mut patches: HashMap<FaceHandle, RectangularGregoryBiNonic> = HashMap::new();
    for sub in &discovery.sub_chains {
        for (key, _) in &sub.rows {
            if !patches.contains_key(&key.face) {
                patches.insert(key.face, rectangular_binonic_on_face(brep, key.face)?);
            }
        }
    }
    // Pin every node position ONCE, globally per vertex, before
    // any chain is fitted: a face corner control is shared by
    // its U-row and V-row, which belong to transverse chains —
    // if each chain pinned the node from its own members' views,
    // the two values would differ by the current gap scale and
    // the later writer would break the earlier chain's trims at
    // the corner.
    let node_positions = global_node_positions(&discovery.sub_chains, &patches);
    let mut authored_count = 0usize;
    let mut fit_declined_count = 0usize;
    let mut fit_decline_reason = "none";
    for sub in &discovery.sub_chains {
        match author_sub_chain(&mut patches, sub, &node_positions) {
            Ok(()) => authored_count += 1,
            Err(reason) => {
                fit_declined_count += 1;
                if fit_decline_reason == "none" {
                    fit_decline_reason = reason;
                }
            }
        }
    }
    for (face, patch) in patches {
        set_rectangular_binonic_on_face(brep, face, patch)?;
    }
    let label = if authored_count > 0 {
        "chain_curve_network_authored"
    } else {
        "chain_curve_network_declined"
    };
    Ok(format!(
        "{label};sub_chain_count={};authored_sub_chain_count={authored_count};fit_declined_sub_chain_count={fit_declined_count};fit_decline_reason={fit_decline_reason};discovery_declined_chain_count={};discovery_decline_reason={}",
        discovery.sub_chains.len(),
        discovery.declined_chain_count,
        discovery.decline_reason
    ))
}

fn author_sub_chain(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub: &SubChain<ChainRowKey, VertexHandle>,
    node_positions: &HashMap<VertexHandle, Point3>,
) -> Result<(), &'static str> {
    let domain = sub
        .rows
        .iter()
        .fold((f64::MAX, f64::MIN), |acc, (_, embedding)| {
            let (lo, hi) = embedding.chain_interval();
            (acc.0.min(lo), acc.1.max(hi))
        });
    let length = domain.1 - domain.0;
    if length <= NODE_PARAM_TOL {
        return Err("sub_chain_domain_degenerate");
    }
    let local = |param: f64| (param - domain.0) / length;
    let rows: Vec<(ChainRowKey, RowChainEmbedding, [Point3; 10])> = sub
        .rows
        .iter()
        .map(|&(key, embedding)| {
            let controls = boundary_controls(&patches[&key.face], key.side);
            (key, embedding, controls)
        })
        .collect();
    let nodes = sub_chain_node_constraints(sub, node_positions);
    if nodes.len() > MAX_INTERPOLATED_NODES {
        return Err("sub_chain_node_count_exceeds_degree");
    }
    let mut samples = Vec::with_capacity(rows.len() * CHAIN_FIT_SAMPLES_PER_ROW);
    for (_, embedding, controls) in &rows {
        for index in 0..CHAIN_FIT_SAMPLES_PER_ROW {
            let t = index as f64 / (CHAIN_FIT_SAMPLES_PER_ROW - 1) as f64;
            samples.push((local(embedding.map(t)), eval_row(controls, t)));
        }
    }
    let constraints: Vec<(f64, Point3)> = nodes
        .iter()
        .map(|&(param, position)| (local(param), position))
        .collect();
    let chain_cps = fit_constrained_degree9(&samples, &constraints)?;
    for (key, embedding, _) in &rows {
        let start = local(embedding.map(0.0));
        let end = local(embedding.map(1.0));
        let row = bezier_row_segment(&chain_cps, start, end);
        let Some(patch) = patches.get_mut(&key.face) else {
            continue;
        };
        for (index, point) in row.iter().enumerate() {
            set_boundary_control(patch, key.side, index, *point);
        }
    }
    Ok(())
}

/// Pinned 3D position per node, computed ONCE across every
/// chain: the average of every member row's view of the node —
/// endpoint control points for rows ending there, curve
/// evaluations for rows spanning it. Keyed by vertex handle so
/// transverse chains meeting at a face corner interpolate the
/// SAME point.
fn global_node_positions(
    sub_chains: &[SubChain<ChainRowKey, VertexHandle>],
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
) -> HashMap<VertexHandle, Point3> {
    let mut sums: HashMap<VertexHandle, (Point3, usize)> = HashMap::new();
    for sub in sub_chains {
        for &(vertex, param) in &sub.nodes {
            for &(key, embedding) in &sub.rows {
                let (lo, hi) = embedding.chain_interval();
                if param < lo - NODE_PARAM_TOL || param > hi + NODE_PARAM_TOL {
                    continue;
                }
                let controls = boundary_controls(&patches[&key.face], key.side);
                let t = (param - embedding.offset) / embedding.scale;
                let position = if t.abs() <= NODE_PARAM_TOL {
                    controls[0]
                } else if (t - 1.0).abs() <= NODE_PARAM_TOL {
                    controls[9]
                } else {
                    eval_row(&controls, t)
                };
                sums.entry(vertex)
                    .and_modify(|(sum, count)| {
                        *sum =
                            Point3::new(sum.x + position.x, sum.y + position.y, sum.z + position.z);
                        *count += 1;
                    })
                    .or_insert((position, 1));
            }
        }
    }
    sums.into_iter()
        .map(|(vertex, (sum, count))| {
            let inverse = 1.0 / count as f64;
            (
                vertex,
                Point3::new(sum.x * inverse, sum.y * inverse, sum.z * inverse),
            )
        })
        .collect()
}

/// Interpolation constraints for one sub-chain: its node
/// parameters paired with the globally pinned vertex positions.
fn sub_chain_node_constraints(
    sub: &SubChain<ChainRowKey, VertexHandle>,
    node_positions: &HashMap<VertexHandle, Point3>,
) -> Vec<(f64, Point3)> {
    let mut out = Vec::with_capacity(sub.nodes.len());
    for &(vertex, param) in &sub.nodes {
        if let Some(&position) = node_positions.get(&vertex) {
            out.push((param, position));
        }
    }
    // Deduplicate parameters (a rebased cut node appears at both
    // window ends with distinct params, which is correct; true
    // duplicates within tolerance collapse to the first).
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out.dedup_by(|a, b| (a.0 - b.0).abs() <= NODE_PARAM_TOL);
    out
}

/// Constrained least squares for one degree-9 Bezier curve:
/// minimize the sample residual subject to exact interpolation
/// of the pinned nodes, via the KKT system
/// `[[AᵀA + λI, Cᵀ], [C, 0]] [p; μ] = [Aᵀb + λ·b̄; d]` solved
/// per axis.
fn fit_constrained_degree9(
    samples: &[(f64, Point3)],
    constraints: &[(f64, Point3)],
) -> Result<[Point3; 10], &'static str> {
    if samples.len() < 10 {
        return Err("sub_chain_sample_count_insufficient");
    }
    let k = constraints.len();
    let n = 10 + k;
    let mut matrix = vec![vec![0.0_f64; n]; n];
    let mut rhs = [vec![0.0_f64; n], vec![0.0_f64; n], vec![0.0_f64; n]];
    for &(t, target) in samples {
        let basis = bernstein9(t);
        for i in 0..10 {
            for j in 0..10 {
                matrix[i][j] += basis[i] * basis[j];
            }
            rhs[0][i] += basis[i] * target.x;
            rhs[1][i] += basis[i] * target.y;
            rhs[2][i] += basis[i] * target.z;
        }
    }
    for (i, row) in matrix.iter_mut().enumerate().take(10) {
        row[i] += FIT_REGULARIZATION;
    }
    for (constraint_index, &(t, position)) in constraints.iter().enumerate() {
        let basis = bernstein9(t);
        let row = 10 + constraint_index;
        for i in 0..10 {
            matrix[row][i] = basis[i];
            matrix[i][row] = basis[i];
        }
        rhs[0][row] = position.x;
        rhs[1][row] = position.y;
        rhs[2][row] = position.z;
    }
    let xs = solve_dense(matrix.clone(), rhs[0].clone()).ok_or("sub_chain_fit_singular")?;
    let ys = solve_dense(matrix.clone(), rhs[1].clone()).ok_or("sub_chain_fit_singular")?;
    let zs = solve_dense(matrix, rhs[2].clone()).ok_or("sub_chain_fit_singular")?;
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = Point3::new(xs[index], ys[index], zs[index]);
    }
    Ok(out)
}

fn solve_dense(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Option<Vec<f64>> {
    let n = rhs.len();
    for pivot in 0..n {
        let mut best = pivot;
        for row in pivot + 1..n {
            if matrix[row][pivot].abs() > matrix[best][pivot].abs() {
                best = row;
            }
        }
        matrix.swap(pivot, best);
        rhs.swap(pivot, best);
        let scale = matrix[pivot][pivot];
        if scale.abs() <= 1.0e-13 {
            return None;
        }
        for column in pivot..n {
            matrix[pivot][column] /= scale;
        }
        rhs[pivot] /= scale;
        for row in 0..n {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            if factor == 0.0 {
                continue;
            }
            for column in pivot..n {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    Some(rhs)
}

fn eval_row(row: &[Point3; 10], t: f64) -> Point3 {
    let basis = bernstein9(t);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for (index, weight) in basis.iter().enumerate() {
        out = Point3::new(
            out.x + row[index].x * weight,
            out.y + row[index].y * weight,
            out.z + row[index].z * weight,
        );
    }
    out
}

fn boundary_controls(patch: &RectangularGregoryBiNonic, side: BoundarySide) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = match side {
            BoundarySide::UMin => patch.grid[0][index],
            BoundarySide::UMax => patch.grid[9][index],
            BoundarySide::VMin => patch.grid[index][0],
            BoundarySide::VMax => patch.grid[index][9],
        };
    }
    out
}

fn set_boundary_control(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    index: usize,
    point: Point3,
) {
    match side {
        BoundarySide::UMin => patch.grid[0][index] = point,
        BoundarySide::UMax => patch.grid[9][index] = point,
        BoundarySide::VMin => patch.grid[index][0] = point,
        BoundarySide::VMax => patch.grid[index][9] = point,
    }
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_chain_fit_tests.rs"]
mod tests;
