//! Globally pinned node registries for the chain network.
//!
//! Both registries follow one principle: every chain incident to
//! a node must agree on the node's data, so it is computed ONCE
//! per vertex from all incident rows and shared by every fit.

use std::collections::HashMap;

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::BoundarySide;
use super::chain_author::boundary_controls;
use super::chain_discovery::SubChain;
use super::chain_fit::ChainRowKey;
use super::row_fit::eval_row;

/// Parameter tolerance for matching a row endpoint to a node.
const NODE_PARAM_TOL: f64 = 1.0e-6;

/// Pinned 3D position per node, computed ONCE across every
/// chain: the average of every member row's view of the node —
/// endpoint control points for rows ending there, curve
/// evaluations for rows spanning it. Keyed by vertex handle so
/// transverse chains meeting at a face corner interpolate the
/// SAME point.
pub(super) fn global_node_positions(
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

/// Globally pinned tangent-plane normal per node, fit to the
/// incident DIRECTION data: every row ending at the node
/// contributes its boundary tangent leg and its transverse
/// (cross-derivative) leg — exactly the vectors the band
/// derivation must place in one plane. The normal is the
/// least-squares plane normal of those directions (sign-aligned
/// cross-product accumulation), so legs from unchained rim rows
/// help DEFINE the plane instead of violating it. Surface-normal
/// sampling is the fallback when a node has too few directions.
pub(super) fn global_node_normals(
    sub_chains: &[SubChain<ChainRowKey, VertexHandle>],
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
) -> HashMap<VertexHandle, Vec3> {
    let mut directions: HashMap<VertexHandle, Vec<Vec3>> = HashMap::new();
    for sub in sub_chains {
        for &(vertex, param) in &sub.nodes {
            for &(key, embedding) in &sub.rows {
                let (lo, hi) = embedding.chain_interval();
                if param < lo - NODE_PARAM_TOL || param > hi + NODE_PARAM_TOL {
                    continue;
                }
                let Some(patch) = patches.get(&key.face) else {
                    continue;
                };
                let row0 = boundary_controls(patch, key.side);
                let row1 = second_row_controls(patch, key.side);
                let t = (param - embedding.offset) / embedding.scale;
                let entry = directions.entry(vertex).or_default();
                if t.abs() <= NODE_PARAM_TOL {
                    push_direction(entry, row1[0] - row0[0]);
                    push_direction(entry, row0[1] - row0[0]);
                } else if (t - 1.0).abs() <= NODE_PARAM_TOL {
                    push_direction(entry, row1[9] - row0[9]);
                    push_direction(entry, row0[8] - row0[9]);
                } else {
                    // The node is interior to this row (a T-stem):
                    // contribute the row tangent at the node.
                    let derivative = eval_row_derivative(&row0, t.clamp(0.0, 1.0));
                    push_direction(entry, derivative);
                }
            }
        }
    }
    let mut normals = HashMap::new();
    for (vertex, vectors) in directions {
        if let Some(normal) = plane_normal_of_directions(&vectors) {
            normals.insert(vertex, normal);
            continue;
        }
        // Fallback: surface normal at the node from any incident
        // member row.
        'fallback: for sub in sub_chains {
            for &(node, param) in &sub.nodes {
                if node != vertex {
                    continue;
                }
                for &(key, embedding) in &sub.rows {
                    let (lo, hi) = embedding.chain_interval();
                    if param < lo - NODE_PARAM_TOL || param > hi + NODE_PARAM_TOL {
                        continue;
                    }
                    let t = ((param - embedding.offset) / embedding.scale).clamp(0.0, 1.0);
                    if let Some(patch) = patches.get(&key.face)
                        && let Some(normal) = boundary_surface_normal(patch, key.side, t)
                    {
                        normals.insert(vertex, normal);
                        break 'fallback;
                    }
                }
            }
        }
    }
    normals
}

fn push_direction(directions: &mut Vec<Vec3>, vector: Vec3) {
    if let Some(unit) = vector.normalized() {
        directions.push(unit);
    }
}

/// Least-squares plane normal of a set of (unit) directions via
/// sign-aligned cross-product accumulation, weighted by pair
/// independence. Returns None when fewer than two independent
/// directions exist.
fn plane_normal_of_directions(directions: &[Vec3]) -> Option<Vec3> {
    // Deterministic accumulation: sort the directions so the
    // sign alignment and averaging order do not depend on map
    // iteration order (which made downstream geometry vary
    // run-to-run).
    let mut directions = directions.to_vec();
    directions.sort_by(|a, b| {
        a.x.total_cmp(&b.x)
            .then(a.y.total_cmp(&b.y))
            .then(a.z.total_cmp(&b.z))
    });
    let directions = &directions;
    let mut accumulated = Vec3::new(0.0, 0.0, 0.0);
    let mut accumulated_any = false;
    for first in 0..directions.len() {
        for second in first + 1..directions.len() {
            let cross = directions[first].cross(&directions[second]);
            let magnitude = cross.dot(&cross).sqrt();
            if magnitude <= 1.0e-6 {
                continue;
            }
            let aligned = if accumulated_any && accumulated.dot(&cross) < 0.0 {
                cross * -1.0
            } else {
                cross
            };
            accumulated = accumulated + aligned;
            accumulated_any = true;
        }
    }
    if !accumulated_any {
        return None;
    }
    accumulated.normalized()
}

/// Derivative of a degree-9 row at `t`.
fn eval_row_derivative(row: &[Point3; 10], t: f64) -> Vec3 {
    let mut delta = [Vec3::new(0.0, 0.0, 0.0); 9];
    for (index, vector) in delta.iter_mut().enumerate() {
        *vector = (row[index + 1] - row[index]) * 9.0;
    }
    // de Casteljau on the derivative control vectors.
    let mut working = delta;
    for level in 1..9 {
        for index in 0..(9 - level) {
            working[index] = working[index] + (working[index + 1] - working[index]) * t;
        }
    }
    working[0]
}

pub(super) fn second_row_controls(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
) -> [Point3; 10] {
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

/// Unit surface normal of the patch at a boundary-row point.
fn boundary_surface_normal(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
    t: f64,
) -> Option<Vec3> {
    let (u, v) = match side {
        BoundarySide::UMin => (0.0, t),
        BoundarySide::UMax => (1.0, t),
        BoundarySide::VMin => (t, 0.0),
        BoundarySide::VMax => (t, 1.0),
    };
    patch
        .partial_u(u, v)
        .cross(&patch.partial_v(u, v))
        .normalized()
}

/// Node identities with their parameters in the sub-chain's
/// local `[0, 1]` domain.
pub(super) fn local_node_params(
    sub: &SubChain<ChainRowKey, VertexHandle>,
    local: &impl Fn(f64) -> f64,
) -> Vec<(VertexHandle, f64)> {
    let mut out: Vec<(VertexHandle, f64)> = sub
        .nodes
        .iter()
        .map(|&(vertex, param)| (vertex, local(param)))
        .collect();
    out.sort_by(|a, b| a.1.total_cmp(&b.1));
    out.dedup_by(|a, b| (a.1 - b.1).abs() <= NODE_PARAM_TOL);
    out
}

/// Interpolation constraints for one sub-chain: its node
/// parameters paired with the globally pinned vertex positions.
pub(super) fn sub_chain_node_constraints(
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

/// Tangent-plane constraints for one sub-chain: at every node
/// with a pinned normal, the curve tangent must lie in the plane
/// perpendicular to it, so all chains incident to the node share
/// one tangent plane and corner cross-derivatives become
/// G1-compatible by construction.
pub(super) fn sub_chain_tangent_constraints(
    sub: &SubChain<ChainRowKey, VertexHandle>,
    node_normals: &HashMap<VertexHandle, Vec3>,
) -> Vec<(f64, Vec3)> {
    let mut out = Vec::with_capacity(sub.nodes.len());
    for &(vertex, param) in &sub.nodes {
        if let Some(&normal) = node_normals.get(&vertex) {
            out.push((param, normal));
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out.dedup_by(|a, b| (a.0 - b.0).abs() <= NODE_PARAM_TOL);
    out
}
