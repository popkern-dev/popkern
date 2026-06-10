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
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::chain_discovery::{ChainTrim, RowChainEmbedding, SubChain, discover_sub_chains};
use super::chain_ribbon::project_sub_chain_rows_to_shared_tangent_planes;
use super::chain_solver::{ChainCurveConstraints, fit_constrained_degree9};
use super::row_fit::bernstein9;
use super::row_segment::bezier_row_segment;
use super::{
    BoundarySide, full_boundary_trims, rectangular_binonic_on_face, set_rectangular_binonic_on_face,
};

/// Samples per member row feeding the least-squares part of the
/// canonical fit.
const CHAIN_FIT_SAMPLES_PER_ROW: usize = 13;
/// Maximum interpolated nodes per sub-chain: a degree-9 curve
/// has 10 control points per axis; pinning more than 9 leaves
/// no least-squares freedom.
const MAX_INTERPOLATED_NODES: usize = 9;
/// Hard-constraint row budget for one sub-chain fit (30 control
/// coordinates; a few rows of slack keep the KKT system well
/// conditioned). Position rows always win; tangent-plane rows
/// are dropped beyond the budget — losing localized G1 at those
/// nodes instead of declining the whole sub-chain (which would
/// cost exact G0).
const MAX_CONSTRAINT_ROWS: usize = 28;
/// Parameter tolerance for matching a row endpoint to a node
/// (chain-parameter units — one row spans 1.0).
const NODE_PARAM_TOL: f64 = 1.0e-6;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) struct ChainRowKey {
    pub(super) face: FaceHandle,
    pub(super) side: BoundarySide,
}

/// Whether the network authors curves only (exact G0) or also
/// projects row-1 controls onto shared tangent planes (G1).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChainNetworkMode {
    Curves,
    CurvesAndCrossTangentRibbons,
}

/// Author canonical curves for every discoverable chain and
/// assign all member rows as exact restrictions. The status
/// records authored/declined chain counts for the candidate
/// score fact.
pub(super) fn author_chain_curve_network(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(brep, slots, ChainNetworkMode::Curves)
}

/// Author the chain-curve network AND project every member
/// row's row-1 controls onto shared cross-tangent planes: chain
/// curves are re-fit with tangent-plane constraints at every
/// pinned node normal, then both sides of every chain receive
/// cross-derivative fields perpendicular to the interpolated
/// normal field — exact G0 (shared curves) plus measured G1
/// (shared tangent planes).
pub(super) fn author_chain_network_with_cross_tangent_ribbons(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(brep, slots, ChainNetworkMode::CurvesAndCrossTangentRibbons)
}

fn author_chain_network_internal(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
    mode: ChainNetworkMode,
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
    // Node normals follow the same global-pinning principle so
    // every chain incident to a node constrains its tangent to
    // the SAME plane. Computed from the pre-network geometry.
    let node_normals = if mode == ChainNetworkMode::CurvesAndCrossTangentRibbons {
        global_node_normals(&discovery.sub_chains, &patches)
    } else {
        HashMap::new()
    };
    let mut authored_count = 0usize;
    let mut fit_declined_count = 0usize;
    let mut fit_decline_reason = "none";
    let mut projected_row_count = 0usize;
    let mut damped_row_count = 0usize;
    let mut skipped_row_count = 0usize;
    let mut ribbon_declined_count = 0usize;
    for sub in &discovery.sub_chains {
        match author_sub_chain(&mut patches, sub, &node_positions, &node_normals) {
            Ok(()) => authored_count += 1,
            Err(reason) => {
                fit_declined_count += 1;
                if fit_decline_reason == "none" {
                    fit_decline_reason = reason;
                }
                continue;
            }
        }
        if mode == ChainNetworkMode::CurvesAndCrossTangentRibbons {
            match project_sub_chain_rows_to_shared_tangent_planes(&mut patches, sub, &node_normals)
            {
                Ok(outcome) => {
                    projected_row_count += outcome.projected_row_count;
                    damped_row_count += outcome.damped_row_count;
                    skipped_row_count += outcome.skipped_row_count;
                }
                Err(_) => ribbon_declined_count += 1,
            }
        }
    }
    for (face, patch) in patches {
        set_rectangular_binonic_on_face(brep, face, patch)?;
    }
    let label = match (mode, authored_count > 0) {
        (ChainNetworkMode::Curves, true) => "chain_curve_network_authored",
        (ChainNetworkMode::Curves, false) => "chain_curve_network_declined",
        (ChainNetworkMode::CurvesAndCrossTangentRibbons, true) => "chain_ribbon_network_authored",
        (ChainNetworkMode::CurvesAndCrossTangentRibbons, false) => "chain_ribbon_network_declined",
    };
    Ok(format!(
        "{label};sub_chain_count={};authored_sub_chain_count={authored_count};fit_declined_sub_chain_count={fit_declined_count};fit_decline_reason={fit_decline_reason};discovery_declined_chain_count={};discovery_decline_reason={};projected_row_count={projected_row_count};damped_row_count={damped_row_count};skipped_row_count={skipped_row_count};ribbon_declined_sub_chain_count={ribbon_declined_count}",
        discovery.sub_chains.len(),
        discovery.declined_chain_count,
        discovery.decline_reason
    ))
}

fn author_sub_chain(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub: &SubChain<ChainRowKey, VertexHandle>,
    node_positions: &HashMap<VertexHandle, Point3>,
    node_normals: &HashMap<VertexHandle, Vec3>,
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
    let mut constraints = ChainCurveConstraints {
        positions: nodes
            .iter()
            .map(|&(param, position)| (local(param), position))
            .collect(),
        tangent_planes: sub_chain_tangent_constraints(sub, node_normals)
            .into_iter()
            .map(|(param, normal)| (local(param), normal))
            .collect(),
    };
    let tangent_budget = MAX_CONSTRAINT_ROWS.saturating_sub(3 * constraints.positions.len());
    constraints.tangent_planes.truncate(tangent_budget);
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

/// Globally pinned tangent-plane normal per node, computed like
/// node positions: every member face contributes its surface
/// normal at the node's boundary point, sign-aligned to the
/// running sum (only the normal's LINE matters — the
/// perpendicularity constraints are sign-agnostic).
fn global_node_normals(
    sub_chains: &[SubChain<ChainRowKey, VertexHandle>],
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
) -> HashMap<VertexHandle, Vec3> {
    let mut sums: HashMap<VertexHandle, Vec3> = HashMap::new();
    for sub in sub_chains {
        for &(vertex, param) in &sub.nodes {
            for &(key, embedding) in &sub.rows {
                let (lo, hi) = embedding.chain_interval();
                if param < lo - NODE_PARAM_TOL || param > hi + NODE_PARAM_TOL {
                    continue;
                }
                let t = ((param - embedding.offset) / embedding.scale).clamp(0.0, 1.0);
                let Some(normal) = boundary_surface_normal(&patches[&key.face], key.side, t) else {
                    continue;
                };
                let entry = sums.entry(vertex).or_insert(Vec3::new(0.0, 0.0, 0.0));
                let aligned = if entry.dot(&normal) < 0.0 {
                    normal * -1.0
                } else {
                    normal
                };
                *entry = *entry + aligned;
            }
        }
    }
    sums.into_iter()
        .filter_map(|(vertex, sum)| sum.normalized().map(|normal| (vertex, normal)))
        .collect()
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

/// Tangent-plane constraints for one sub-chain: at every node
/// with a pinned normal, the curve tangent must lie in the plane
/// perpendicular to it, so all chains incident to the node share
/// one tangent plane and corner cross-derivatives become
/// G1-compatible by construction.
fn sub_chain_tangent_constraints(
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

pub(super) fn eval_row(row: &[Point3; 10], t: f64) -> Point3 {
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

pub(super) fn boundary_controls(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
) -> [Point3; 10] {
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
