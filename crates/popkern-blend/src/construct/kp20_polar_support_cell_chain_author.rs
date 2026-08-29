//! Single sub-chain curve authoring for the KP20 chain network.

use std::collections::HashMap;

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::BoundarySide;
use super::cell_reauthor::AuthoredSubChain;
use super::chain_discovery::{RowChainEmbedding, SubChain};
use super::chain_fit::ChainRowKey;
use super::chain_nodes::{
    local_node_params, sub_chain_node_constraints, sub_chain_tangent_constraints,
};
use super::chain_solver::{ChainCurveConstraints, fit_constrained_degree9};
use super::chain_source::{ChainCurveSource, ChainCurveSourcePreference};
use super::chain_spacing::{end_decompression_targets, node_speed_decompression_targets};
use super::row_fit::eval_row;
use super::row_segment::bezier_row_segment;

const CHAIN_FIT_SAMPLES_PER_ROW: usize = 13;
const MAX_INTERPOLATED_NODES: usize = 9;
const NODE_PARAM_TOL: f64 = 1.0e-6;

pub(super) fn author_sub_chain(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub: &SubChain<ChainRowKey, VertexHandle>,
    node_positions: &HashMap<VertexHandle, Point3>,
    node_normals: &HashMap<VertexHandle, Vec3>,
    leg_class_targets: &[(f64, Vec3, f64)],
    source_preference: ChainCurveSourcePreference,
) -> Result<AuthoredSubChain, &'static str> {
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
    let nodes_local = local_node_params(sub, &local);
    if source_preference.prefers_spanning_owner_row()
        && let Some(chain_cps) = spanning_owner_row_controls(&rows, domain)
    {
        return Ok(authored_sub_chain_from_controls(
            patches,
            sub,
            &local,
            chain_cps,
            nodes_local,
            ChainCurveSource::SpanningOwnerRow,
        ));
    }
    let mut samples = Vec::with_capacity(rows.len() * CHAIN_FIT_SAMPLES_PER_ROW);
    for (_, embedding, controls) in &rows {
        for index in 0..CHAIN_FIT_SAMPLES_PER_ROW {
            let t = index as f64 / (CHAIN_FIT_SAMPLES_PER_ROW - 1) as f64;
            samples.push((local(embedding.map(t)), eval_row(controls, t)));
        }
    }
    let constraints = ChainCurveConstraints {
        positions: nodes
            .iter()
            .map(|&(param, position)| (local(param), position))
            .collect(),
        tangent_planes: sub_chain_tangent_constraints(sub, node_normals)
            .into_iter()
            .map(|(param, normal)| (local(param), normal))
            .collect(),
        spacing_targets: Vec::new(),
        derivative_targets: Vec::new(),
    };
    let mut chain_cps = fit_constrained_degree9(&samples, &constraints)?;
    let spacing = end_decompression_targets(&chain_cps);
    let mut derivatives = node_speed_decompression_targets(&chain_cps, &nodes_local);
    derivatives.extend_from_slice(leg_class_targets);
    if spacing.is_some() || !derivatives.is_empty() {
        let mut retargeted = constraints.clone();
        retargeted.spacing_targets = spacing.unwrap_or_default();
        retargeted.derivative_targets = derivatives;
        chain_cps = fit_constrained_degree9(&samples, &retargeted)?;
    }
    let authored = authored_sub_chain_from_controls(
        patches,
        sub,
        &local,
        chain_cps,
        nodes_local,
        ChainCurveSource::LeastSquaresFit,
    );
    Ok(authored)
}

pub(super) fn authored_sub_chain_from_controls(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub: &SubChain<ChainRowKey, VertexHandle>,
    local: &impl Fn(f64) -> f64,
    chain_cps: [Point3; 10],
    nodes: Vec<(VertexHandle, f64)>,
    source: ChainCurveSource,
) -> AuthoredSubChain {
    let mut authored = AuthoredSubChain {
        chain_cps,
        rows: Vec::with_capacity(sub.rows.len()),
        nodes,
        source,
    };
    for &(key, embedding) in &sub.rows {
        let start = local(embedding.map(0.0));
        let end = local(embedding.map(1.0));
        let row = bezier_row_segment(&chain_cps, start, end);
        let Some(patch) = patches.get_mut(&key.face) else {
            continue;
        };
        for (index, point) in row.iter().enumerate() {
            set_boundary_control(patch, key.side, index, *point);
        }
        authored.rows.push((key, (start, end)));
    }
    authored
}

fn spanning_owner_row_controls(
    rows: &[(ChainRowKey, RowChainEmbedding, [Point3; 10])],
    domain: (f64, f64),
) -> Option<[Point3; 10]> {
    rows.iter().find_map(|(_, embedding, controls)| {
        let (lo, hi) = embedding.chain_interval();
        if lo > domain.0 + NODE_PARAM_TOL || hi < domain.1 - NODE_PARAM_TOL {
            return None;
        }
        let mut source = *controls;
        if embedding.scale < 0.0 {
            source.reverse();
        }
        Some(source)
    })
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
