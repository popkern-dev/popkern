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

use super::cell_reauthor::{AuthoredSubChain, reauthor_cells_from_chain_network};
use super::chain_discovery::{ChainTrim, RowChainEmbedding, SubChain, discover_sub_chains};
use super::chain_nodes::{
    global_node_normals, global_node_positions, local_node_params, sub_chain_node_constraints,
    sub_chain_tangent_constraints,
};
use super::chain_ribbon::project_sub_chain_rows_to_shared_tangent_planes;
use super::chain_solver::{ChainCurveConstraints, fit_constrained_degree9};
use super::chain_spacing::{
    end_decompression_targets, end_spacing_ratio, node_speed_decompression_targets,
};
use super::leg_class::leg_class_targets;
use super::row_fit::bernstein9;
use super::row_segment::bezier_row_segment;
use super::{
    BoundarySide, count_skipped_boundary_trims, full_boundary_trims, rectangular_binonic_on_face,
    set_rectangular_binonic_dual_twist_on_face, set_rectangular_binonic_on_face,
};

/// Samples per member row feeding the least-squares part of the
/// canonical fit.
const CHAIN_FIT_SAMPLES_PER_ROW: usize = 13;
/// Maximum interpolated nodes per sub-chain: a degree-9 curve
/// has 10 control points per axis; pinning more than 9 leaves
/// no least-squares freedom.
const MAX_INTERPOLATED_NODES: usize = 9;
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
    CurvesAndConstructiveCells,
    /// Constructive cells emitted as dual-twist patches: each
    /// boundary's derived corner band coefficient lands in its
    /// own Chiyokura-Kimura candidate slot — twist conflicts
    /// cease to exist instead of being compromised.
    CurvesAndConstructiveDualTwistCells,
}

impl ChainNetworkMode {
    fn constructive(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveCells
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCells
        )
    }
}

/// Author canonical curves for every discoverable chain and
/// assign all member rows as exact restrictions.
pub(super) fn author_chain_curve_network(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(brep, slots, ChainNetworkMode::Curves)
}

/// Chain network plus row-1 projection onto shared
/// cross-tangent planes (exact G0, measured G1).
pub(super) fn author_chain_network_with_cross_tangent_ribbons(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(brep, slots, ChainNetworkMode::CurvesAndCrossTangentRibbons)
}

/// Chain network plus constructive cell re-authoring: per-chain
/// ribbon fields, derived row-1 bands, Coons interiors.
pub(super) fn author_chain_network_with_constructive_cells(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(brep, slots, ChainNetworkMode::CurvesAndConstructiveCells)
}

/// Constructive cells on the dual-twist patch substrate: the
/// per-boundary corner candidates eliminate the contested
/// single twist control entirely (design:
/// `docs/decisions/2026-06-11-kp20-binonic-dual-twist-corners.md`).
pub(super) fn author_chain_network_with_constructive_dual_twist_cells(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(
        brep,
        slots,
        ChainNetworkMode::CurvesAndConstructiveDualTwistCells,
    )
}

fn author_chain_network_internal(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
    mode: ChainNetworkMode,
) -> Result<String, CoreError> {
    let trims = full_boundary_trims(brep, slots)?;
    let skipped_trims = count_skipped_boundary_trims(brep, slots)?;
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
    // Row-conservation observability: every row glued by a trim
    // must be assigned to exactly one sub-chain; a row silently
    // dropped by the clean-node splitter keeps its old geometry
    // and surfaces as a localized G0/G1 hole.
    let mut glued_rows: std::collections::HashSet<ChainRowKey> = std::collections::HashSet::new();
    for trim in &chain_trims {
        glued_rows.insert(trim.row_a);
        glued_rows.insert(trim.row_b);
    }
    let mut assigned_rows: std::collections::HashSet<ChainRowKey> =
        std::collections::HashSet::new();
    let mut multiply_assigned_row_count = 0usize;
    for sub in &discovery.sub_chains {
        for (key, _) in &sub.rows {
            if !assigned_rows.insert(*key) {
                multiply_assigned_row_count += 1;
            }
        }
    }
    let dropped_rows: Vec<ChainRowKey> = glued_rows
        .iter()
        .filter(|key| !assigned_rows.contains(key))
        .copied()
        .collect();
    // Shared-field invariant: a trim's two rows must land in the
    // SAME sub-chain, else they are derived from different
    // ribbon fields and their tangent planes cannot coincide.
    let mut row_sub_chain: HashMap<ChainRowKey, usize> = HashMap::new();
    for (index, sub) in discovery.sub_chains.iter().enumerate() {
        for (key, _) in &sub.rows {
            row_sub_chain.insert(*key, index);
        }
    }
    let cross_sub_chain_trim_count = chain_trims
        .iter()
        .filter(|trim| {
            match (
                row_sub_chain.get(&trim.row_a),
                row_sub_chain.get(&trim.row_b),
            ) {
                (Some(a), Some(b)) => a != b,
                _ => true,
            }
        })
        .count();
    let dropped_row_list = dropped_rows
        .iter()
        .take(6)
        .map(|key| format!("{:?}:{:?}", key.face, key.side))
        .collect::<Vec<_>>()
        .join(",");
    let mut patches: HashMap<FaceHandle, RectangularGregoryBiNonic> = HashMap::new();
    for sub in &discovery.sub_chains {
        for (key, _) in &sub.rows {
            if !patches.contains_key(&key.face) {
                patches.insert(key.face, rectangular_binonic_on_face(brep, key.face)?);
            }
        }
    }
    // Coverage observability: a slot face with no chained row is
    // invisible to the entire network (its bands stay underived),
    // which surfaces as a localized G1 hole — name the faces so
    // a coverage regression cannot hide inside an aggregate.
    let uncovered_faces: Vec<FaceHandle> = slots
        .iter()
        .map(|slot| slot.face)
        .filter(|face| !patches.contains_key(face))
        .collect();
    let uncovered_face_list = uncovered_faces
        .iter()
        .take(8)
        .map(|face| format!("{face:?}"))
        .collect::<Vec<_>>()
        .join(",");
    let chained_face_count = patches.len();
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
    let node_normals = if mode == ChainNetworkMode::Curves {
        HashMap::new()
    } else {
        global_node_normals(&discovery.sub_chains, &patches)
    };
    let mut authored_count = 0usize;
    let mut fit_declined_count = 0usize;
    let mut fit_decline_reason = "none";
    let mut projected_row_count = 0usize;
    let mut damped_row_count = 0usize;
    let mut skipped_row_count = 0usize;
    let mut ribbon_declined_count = 0usize;
    let mut authored_sub_chains = Vec::new();
    let mut min_end_spacing_ratio = f64::INFINITY;
    for sub in &discovery.sub_chains {
        match author_sub_chain(&mut patches, sub, &node_positions, &node_normals, &[]) {
            Ok(authored) => {
                authored_count += 1;
                min_end_spacing_ratio =
                    min_end_spacing_ratio.min(end_spacing_ratio(&authored.chain_cps));
                if mode.constructive() {
                    authored_sub_chains.push(authored);
                }
            }
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
    let mut leg_class_thin_count = 0usize;
    let mut leg_class_chain_count = 0usize;
    if mode.constructive() {
        // Per-vertex leg-class solve — see the leg_class module.
        let (mut class_targets, thin_count) = leg_class_targets(&authored_sub_chains);
        leg_class_thin_count += thin_count;
        let mut chain_order: Vec<usize> = class_targets.keys().copied().collect();
        chain_order.sort_unstable();
        for chain_index in chain_order {
            let mut targets = class_targets.remove(&chain_index).unwrap_or_default();
            targets.sort_by(|a, b| a.0.total_cmp(&b.0));
            let sub = &discovery.sub_chains[chain_index];
            if let Ok(authored) =
                author_sub_chain(&mut patches, sub, &node_positions, &node_normals, &targets)
            {
                authored_sub_chains[chain_index] = authored;
                leg_class_chain_count += 1;
            }
        }
    }
    let mut reauthor_detail = String::new();
    let mut dual_twist_pairs = HashMap::new();
    if mode.constructive() {
        let outcome = reauthor_cells_from_chain_network(
            &mut patches,
            &authored_sub_chains,
            &node_normals,
            mode == ChainNetworkMode::CurvesAndConstructiveDualTwistCells,
        );
        dual_twist_pairs = outcome.dual_twist_pairs.clone();
        reauthor_detail = format!(
            ";derived_row_count={};leg_solve_declined_count={};field_declined_sub_chain_count={};twist_conflict_count={};max_twist_conflict_distance={:.6e};max_derived_span_residual={:.6e};worst_span_residual_row={};repaired_leg_count={};max_leg_out_of_plane_ratio={:.6e};max_twist_conflict_normal_component={:.6e};max_relative_twist_displacement={:.6e};refilled_cell_count={};predicted_worst_trim_g1={:.6e};clamped_leg_beta_count={};one_sided_band_count={};dual_multi_demand_side_count={};dual_twist_pair_count={};max_dual_twist_pair_separation={:.6e};predicted_worst_detail=[{}]",
            outcome.derived_row_count,
            outcome.leg_solve_declined_count,
            outcome.field_declined_sub_chain_count,
            outcome.twist_conflict_count,
            outcome.max_twist_conflict_distance,
            outcome.max_derived_span_residual,
            outcome.worst_span_residual_row.as_deref().unwrap_or("none"),
            outcome.repaired_leg_count,
            outcome.max_leg_out_of_plane_ratio,
            outcome.max_twist_conflict_normal_component,
            outcome.max_relative_twist_displacement,
            outcome.refilled_cell_count,
            outcome.predicted_worst_trim_g1,
            outcome.clamped_leg_beta_count,
            outcome.one_sided_band_count,
            outcome.dual_multi_demand_side_count,
            outcome.dual_twist_pair_count,
            outcome.max_dual_twist_pair_separation,
            outcome.predicted_worst_detail.as_deref().unwrap_or("none")
        );
    }
    for (face, patch) in patches {
        match dual_twist_pairs.get(&face) {
            Some(corners) if corners.iter().any(Option::is_some) => {
                set_rectangular_binonic_dual_twist_on_face(brep, face, patch, corners)?;
            }
            _ => set_rectangular_binonic_on_face(brep, face, patch)?,
        }
    }
    let label = match (mode, authored_count > 0) {
        (ChainNetworkMode::Curves, true) => "chain_curve_network_authored",
        (ChainNetworkMode::Curves, false) => "chain_curve_network_declined",
        (ChainNetworkMode::CurvesAndCrossTangentRibbons, true) => "chain_ribbon_network_authored",
        (ChainNetworkMode::CurvesAndCrossTangentRibbons, false) => "chain_ribbon_network_declined",
        (ChainNetworkMode::CurvesAndConstructiveCells, true) => "chain_constructive_cells_authored",
        (ChainNetworkMode::CurvesAndConstructiveCells, false) => {
            "chain_constructive_cells_declined"
        }
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCells, true) => {
            "chain_constructive_dual_twist_cells_authored"
        }
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCells, false) => {
            "chain_constructive_dual_twist_cells_declined"
        }
    };
    Ok(format!(
        "{label};sub_chain_count={};authored_sub_chain_count={authored_count};fit_declined_sub_chain_count={fit_declined_count};fit_decline_reason={fit_decline_reason};discovery_declined_chain_count={};discovery_decline_reason={};chained_face_count={};uncovered_face_count={};uncovered_faces=[{uncovered_face_list}];glued_trim_count={};skipped_trim_missing_pcurve_count={};skipped_trim_unrecognized_span_count={};glued_row_count={};assigned_row_count={};dropped_row_count={};dropped_rows=[{dropped_row_list}];multiply_assigned_row_count={multiply_assigned_row_count};cross_sub_chain_trim_count={cross_sub_chain_trim_count};projected_row_count={projected_row_count};damped_row_count={damped_row_count};skipped_row_count={skipped_row_count};ribbon_declined_sub_chain_count={ribbon_declined_count};min_end_spacing_ratio={min_end_spacing_ratio:.4e};leg_class_thin_count={leg_class_thin_count};leg_class_chain_count={leg_class_chain_count}{reauthor_detail}",
        discovery.sub_chains.len(),
        discovery.declined_chain_count,
        discovery.decline_reason,
        chained_face_count,
        uncovered_faces.len(),
        trims.len(),
        skipped_trims.missing_pcurve,
        skipped_trims.unrecognized_boundary_span,
        glued_rows.len(),
        assigned_rows.len(),
        dropped_rows.len()
    ))
}

fn author_sub_chain(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub: &SubChain<ChainRowKey, VertexHandle>,
    node_positions: &HashMap<VertexHandle, Point3>,
    node_normals: &HashMap<VertexHandle, Vec3>,
    leg_class_targets: &[(f64, Vec3, f64)],
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
    let nodes_local = local_node_params(sub, &local);
    let spacing = end_decompression_targets(&chain_cps);
    let mut derivatives = node_speed_decompression_targets(&chain_cps, &nodes_local);
    derivatives.extend_from_slice(leg_class_targets);
    if spacing.is_some() || !derivatives.is_empty() {
        let mut retargeted = constraints.clone();
        retargeted.spacing_targets = spacing.unwrap_or_default();
        retargeted.derivative_targets = derivatives;
        chain_cps = fit_constrained_degree9(&samples, &retargeted)?;
    }
    let mut authored = AuthoredSubChain {
        chain_cps,
        rows: Vec::with_capacity(rows.len()),
        nodes: nodes_local,
    };
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
        authored.rows.push((*key, (start, end)));
    }
    Ok(authored)
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
