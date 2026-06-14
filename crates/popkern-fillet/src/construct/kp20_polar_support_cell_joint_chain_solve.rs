//! Coupled chain solve for shared endpoint leg consistency.
//!
//! This keeps the accepted single-chain invariants — exact node
//! interpolation, soft-heavy tangent planes, end-spacing and
//! node-speed decompression — but assembles one system for every
//! authored sub-chain so endpoint-leg equations can connect two
//! chain derivatives directly.

use std::collections::HashMap;

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::cell_reauthor::AuthoredSubChain;
use super::chain_author::{authored_sub_chain_from_controls, boundary_controls};
use super::chain_discovery::{RowChainEmbedding, SubChain};
use super::chain_fit::ChainRowKey;
use super::chain_nodes::{
    local_node_params, sub_chain_node_constraints, sub_chain_tangent_constraints,
};
use super::chain_solver::{
    FIT_REGULARIZATION, FIT_SMOOTHING_WEIGHT, SOFT_TANGENT_PLANE_WEIGHT,
    add_third_difference_penalty, dbernstein9, solve_dense,
};
use super::chain_spacing::{end_decompression_targets, node_speed_decompression_targets};
use super::leg_vector::{
    JointEndpointLeg, JointEndpointLegEquation, JointEndpointPlaneConstraint, chain_derivative,
};
use super::row_fit::{bernstein9, eval_row};

const CHAIN_FIT_SAMPLES_PER_ROW: usize = 13;
const MAX_INTERPOLATED_NODES: usize = 9;
const NODE_PARAM_TOL: f64 = 1.0e-6;
const JOINT_ENDPOINT_LEG_WEIGHT: f64 = 50.0;
const JOINT_ENDPOINT_PLANE_WEIGHT: f64 = 1.0e4;

#[derive(Debug, Clone)]
pub(super) struct JointChainSolveOutcome {
    pub(super) authored_sub_chains: Vec<AuthoredSubChain>,
    pub(super) plane_constraint_count: usize,
    pub(super) max_plane_residual_before: f64,
    pub(super) max_plane_residual_after: f64,
}

struct PreparedSubChain {
    samples: Vec<(f64, Point3)>,
    positions: Vec<(f64, Point3)>,
    tangent_planes: Vec<(f64, Vec3)>,
    spacing_targets: Vec<(usize, Point3, f64)>,
    derivative_targets: Vec<(f64, Vec3, f64)>,
    nodes: Vec<(VertexHandle, f64)>,
    domain_start: f64,
    domain_length: f64,
}

pub(super) fn solve_joint_chain_network(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub_chains: &[SubChain<ChainRowKey, VertexHandle>],
    initial_authored: &[AuthoredSubChain],
    node_positions: &HashMap<VertexHandle, Point3>,
    node_normals: &HashMap<VertexHandle, Vec3>,
    derivative_targets: &HashMap<usize, Vec<(f64, Vec3, f64)>>,
    endpoint_equations: &[JointEndpointLegEquation],
    plane_constraints: &[JointEndpointPlaneConstraint],
) -> Result<JointChainSolveOutcome, &'static str> {
    if endpoint_equations.is_empty() && plane_constraints.is_empty() {
        return Err("joint_chain_solve_has_no_joint_constraints");
    }
    if sub_chains.len() != initial_authored.len() {
        return Err("joint_chain_solve_requires_all_sub_chains_authored");
    }
    let mut prepared = Vec::with_capacity(sub_chains.len());
    for (index, sub) in sub_chains.iter().enumerate() {
        prepared.push(prepare_sub_chain(
            patches,
            sub,
            &initial_authored[index],
            node_positions,
            node_normals,
            derivative_targets
                .get(&index)
                .map(Vec::as_slice)
                .unwrap_or(&[]),
        )?);
    }
    let solution = solve_prepared_network(&prepared, endpoint_equations, plane_constraints)?;
    let mut authored_sub_chains = Vec::with_capacity(sub_chains.len());
    for (chain_index, sub) in sub_chains.iter().enumerate() {
        let local = |param: f64| {
            (param - prepared[chain_index].domain_start) / prepared[chain_index].domain_length
        };
        let nodes = prepared[chain_index].nodes.clone();
        authored_sub_chains.push(authored_sub_chain_from_controls(
            patches,
            sub,
            &local,
            solution[chain_index],
            nodes,
        ));
    }
    let max_plane_before = plane_constraints
        .iter()
        .map(|constraint| constraint.residual_before)
        .fold(0.0_f64, f64::max);
    let max_plane_after = plane_constraints
        .iter()
        .map(|constraint| plane_residual(&authored_sub_chains, constraint))
        .fold(0.0_f64, f64::max);
    Ok(JointChainSolveOutcome {
        authored_sub_chains,
        plane_constraint_count: plane_constraints.len(),
        max_plane_residual_before: max_plane_before,
        max_plane_residual_after: max_plane_after,
    })
}

fn prepare_sub_chain(
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub: &SubChain<ChainRowKey, VertexHandle>,
    initial_authored: &AuthoredSubChain,
    node_positions: &HashMap<VertexHandle, Point3>,
    node_normals: &HashMap<VertexHandle, Vec3>,
    derivative_targets: &[(f64, Vec3, f64)],
) -> Result<PreparedSubChain, &'static str> {
    let domain = sub
        .rows
        .iter()
        .fold((f64::MAX, f64::MIN), |acc, (_, embedding)| {
            let (lo, hi) = embedding.chain_interval();
            (acc.0.min(lo), acc.1.max(hi))
        });
    let length = domain.1 - domain.0;
    if length <= NODE_PARAM_TOL {
        return Err("joint_chain_domain_degenerate");
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
        return Err("joint_chain_node_count_exceeds_degree");
    }
    let mut samples = Vec::with_capacity(rows.len() * CHAIN_FIT_SAMPLES_PER_ROW);
    for (_, embedding, controls) in &rows {
        for index in 0..CHAIN_FIT_SAMPLES_PER_ROW {
            let t = index as f64 / (CHAIN_FIT_SAMPLES_PER_ROW - 1) as f64;
            samples.push((local(embedding.map(t)), eval_row(controls, t)));
        }
    }
    if samples.len() < 10 {
        return Err("joint_chain_sample_count_insufficient");
    }
    let mut derivative_rows =
        node_speed_decompression_targets(&initial_authored.chain_cps, &initial_authored.nodes);
    derivative_rows.extend_from_slice(derivative_targets);
    Ok(PreparedSubChain {
        samples,
        positions: nodes
            .iter()
            .map(|&(param, position)| (local(param), position))
            .collect(),
        tangent_planes: sub_chain_tangent_constraints(sub, node_normals)
            .into_iter()
            .map(|(param, normal)| (local(param), normal))
            .collect(),
        spacing_targets: end_decompression_targets(&initial_authored.chain_cps).unwrap_or_default(),
        derivative_targets: derivative_rows,
        nodes: local_node_params(sub, &local),
        domain_start: domain.0,
        domain_length: length,
    })
}

fn solve_prepared_network(
    chains: &[PreparedSubChain],
    endpoint_equations: &[JointEndpointLegEquation],
    plane_constraints: &[JointEndpointPlaneConstraint],
) -> Result<Vec<[Point3; 10]>, &'static str> {
    let hard_rows = chains
        .iter()
        .map(|chain| 3 * chain.positions.len())
        .sum::<usize>();
    let unknowns = 30 * chains.len();
    let n = unknowns + hard_rows;
    let mut matrix = vec![vec![0.0_f64; n]; n];
    let mut rhs = vec![0.0_f64; n];
    for (chain_index, chain) in chains.iter().enumerate() {
        add_chain_objective(&mut matrix, &mut rhs, chain_index, chain);
    }
    for equation in endpoint_equations {
        add_joint_endpoint_equation(&mut matrix, chains.len(), equation)?;
    }
    for constraint in plane_constraints {
        add_joint_endpoint_plane_constraint(&mut matrix, chains.len(), constraint)?;
    }
    let mut row = unknowns;
    for (chain_index, chain) in chains.iter().enumerate() {
        for &(t, position) in &chain.positions {
            let basis = bernstein9(t);
            for (axis, value) in [(0, position.x), (1, position.y), (2, position.z)] {
                for (cp, weight) in basis.iter().enumerate() {
                    let column = variable_index(chain_index, axis, cp);
                    matrix[row][column] = *weight;
                    matrix[column][row] = *weight;
                }
                rhs[row] = value;
                row += 1;
            }
        }
    }
    debug_assert_eq!(row, n);
    let solution = solve_dense(matrix, rhs).ok_or("joint_chain_solve_singular")?;
    let mut out = Vec::with_capacity(chains.len());
    for chain_index in 0..chains.len() {
        let mut controls = [Point3::new(0.0, 0.0, 0.0); 10];
        for (cp, point) in controls.iter_mut().enumerate() {
            *point = Point3::new(
                solution[variable_index(chain_index, 0, cp)],
                solution[variable_index(chain_index, 1, cp)],
                solution[variable_index(chain_index, 2, cp)],
            );
        }
        out.push(controls);
    }
    Ok(out)
}

fn add_chain_objective(
    matrix: &mut [Vec<f64>],
    rhs: &mut [f64],
    chain_index: usize,
    chain: &PreparedSubChain,
) {
    for &(t, target) in &chain.samples {
        let basis = bernstein9(t);
        for axis in 0..3 {
            let value = match axis {
                0 => target.x,
                1 => target.y,
                _ => target.z,
            };
            for i in 0..10 {
                let row = variable_index(chain_index, axis, i);
                rhs[row] += basis[i] * value;
                for j in 0..10 {
                    matrix[row][variable_index(chain_index, axis, j)] += basis[i] * basis[j];
                }
            }
        }
    }
    for axis in 0..3 {
        for cp in 0..10 {
            matrix[variable_index(chain_index, axis, cp)][variable_index(chain_index, axis, cp)] +=
                FIT_REGULARIZATION;
        }
        add_third_difference_penalty(
            matrix,
            variable_index(chain_index, axis, 0),
            10,
            FIT_SMOOTHING_WEIGHT,
        );
    }
    for &(cp, target, weight) in &chain.spacing_targets {
        for (axis, value) in [(0, target.x), (1, target.y), (2, target.z)] {
            let index = variable_index(chain_index, axis, cp);
            matrix[index][index] += weight;
            rhs[index] += weight * value;
        }
    }
    for &(t, target, weight) in &chain.derivative_targets {
        add_derivative_target(matrix, rhs, chain_index, t, target, weight);
    }
    for &(t, normal) in &chain.tangent_planes {
        add_tangent_plane(matrix, chain_index, t, normal);
    }
}

fn add_derivative_target(
    matrix: &mut [Vec<f64>],
    rhs: &mut [f64],
    chain_index: usize,
    t: f64,
    target: Vec3,
    weight: f64,
) {
    let dbasis = dbernstein9(t);
    for (axis, value) in [(0, target.x), (1, target.y), (2, target.z)] {
        for i in 0..10 {
            let row = variable_index(chain_index, axis, i);
            rhs[row] += weight * dbasis[i] * value;
            for j in 0..10 {
                matrix[row][variable_index(chain_index, axis, j)] += weight * dbasis[i] * dbasis[j];
            }
        }
    }
}

fn add_tangent_plane(matrix: &mut [Vec<f64>], chain_index: usize, t: f64, normal: Vec3) {
    let dbasis = dbernstein9(t);
    let mut terms = Vec::with_capacity(30);
    for cp in 0..10 {
        terms.push((variable_index(chain_index, 0, cp), dbasis[cp] * normal.x));
        terms.push((variable_index(chain_index, 1, cp), dbasis[cp] * normal.y));
        terms.push((variable_index(chain_index, 2, cp), dbasis[cp] * normal.z));
    }
    add_zero_target_row(matrix, &terms, SOFT_TANGENT_PLANE_WEIGHT);
}

fn add_joint_endpoint_equation(
    matrix: &mut [Vec<f64>],
    chain_count: usize,
    equation: &JointEndpointLegEquation,
) -> Result<(), &'static str> {
    if equation.first.chain_index >= chain_count || equation.second.chain_index >= chain_count {
        return Err("joint_chain_endpoint_index_out_of_range");
    }
    if equation.first.derivative_scale.abs() <= 1.0e-12
        || equation.second.derivative_scale.abs() <= 1.0e-12
    {
        return Err("joint_chain_endpoint_scale_degenerate");
    }
    let first = endpoint_terms(equation.first, 1.0);
    let second = endpoint_terms(equation.second, -equation.second_orientation);
    for axis in 0..3 {
        let mut terms = Vec::with_capacity(20);
        terms.extend(first[axis].iter().copied());
        terms.extend(second[axis].iter().copied());
        add_zero_target_row(matrix, &terms, JOINT_ENDPOINT_LEG_WEIGHT);
    }
    Ok(())
}

fn add_joint_endpoint_plane_constraint(
    matrix: &mut [Vec<f64>],
    chain_count: usize,
    constraint: &JointEndpointPlaneConstraint,
) -> Result<(), &'static str> {
    if constraint.endpoint.chain_index >= chain_count {
        return Err("joint_chain_endpoint_index_out_of_range");
    }
    if constraint.endpoint.derivative_scale.abs() <= 1.0e-12 {
        return Err("joint_chain_endpoint_scale_degenerate");
    }
    let terms_by_axis = endpoint_terms(constraint.endpoint, 1.0);
    let mut terms = Vec::with_capacity(30);
    for (axis, normal_component) in [
        (0, constraint.plane_normal.x),
        (1, constraint.plane_normal.y),
        (2, constraint.plane_normal.z),
    ] {
        for &(index, coefficient) in &terms_by_axis[axis] {
            terms.push((index, coefficient * normal_component));
        }
    }
    add_zero_target_row(matrix, &terms, JOINT_ENDPOINT_PLANE_WEIGHT);
    Ok(())
}

fn endpoint_terms(endpoint: JointEndpointLeg, sign: f64) -> [Vec<(usize, f64)>; 3] {
    let dbasis = dbernstein9(endpoint.station);
    std::array::from_fn(|axis| {
        (0..10)
            .map(|cp| {
                (
                    variable_index(endpoint.chain_index, axis, cp),
                    sign * dbasis[cp] / endpoint.derivative_scale,
                )
            })
            .collect()
    })
}

fn add_zero_target_row(matrix: &mut [Vec<f64>], terms: &[(usize, f64)], weight: f64) {
    for &(i, coeff_i) in terms {
        if coeff_i == 0.0 {
            continue;
        }
        for &(j, coeff_j) in terms {
            if coeff_j != 0.0 {
                matrix[i][j] += weight * coeff_i * coeff_j;
            }
        }
    }
}

fn plane_residual(
    authored_sub_chains: &[AuthoredSubChain],
    constraint: &JointEndpointPlaneConstraint,
) -> f64 {
    authored_endpoint_leg(authored_sub_chains, constraint.endpoint)
        .dot(&constraint.plane_normal)
        .abs()
}

fn authored_endpoint_leg(
    authored_sub_chains: &[AuthoredSubChain],
    endpoint: JointEndpointLeg,
) -> Vec3 {
    chain_derivative(
        &authored_sub_chains[endpoint.chain_index].chain_cps,
        endpoint.station,
    ) * (1.0 / endpoint.derivative_scale)
}

fn variable_index(chain: usize, axis: usize, cp: usize) -> usize {
    chain * 30 + axis * 10 + cp
}
