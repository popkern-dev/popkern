//! Canonical chain-curve authoring for KP20 support-cell rows:
//! one constrained degree-9 curve per sub-chain (hard node
//! pins, soft tangent planes), every member row its exact de
//! Casteljau restriction — all trims of a chain become CP
//! identities.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::{BRep, FaceHandle};

use super::cell_reauthor::reauthor_cells_from_chain_network;
#[path = "kp20_polar_support_cell_chain_status.rs"]
mod chain_status;
use super::chain_author::author_sub_chain;
use super::chain_discovery::{ChainTrim, discover_sub_chains};
use super::chain_nodes::{global_node_normals, global_node_positions};
use super::chain_ribbon::project_sub_chain_rows_to_shared_tangent_planes;
use super::chain_spacing::end_spacing_ratio;
use super::joint_chain_solve::solve_joint_chain_network;
use super::leg_class::leg_class_targets;
use super::leg_direction::leg_direction_targets;
use super::leg_vector::{
    joint_endpoint_leg_equations, joint_endpoint_plane_constraints, leg_vector_targets,
};
use super::{
    BoundarySide, count_skipped_boundary_trims, full_boundary_trims, rectangular_binonic_on_face,
    set_rectangular_binonic_dual_twist_on_face, set_rectangular_binonic_on_face,
};

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
    CurvesAndConstructiveDualTwistCellsWithJointLegSolve,
    CurvesAndConstructiveDualTwistCellsWithBetaEnvelope,
    CurvesAndConstructiveDualTwistCellsWithPairNormalSolve,
    CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve,
}

impl ChainNetworkMode {
    fn constructive(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveCells
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCells
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve
        )
    }

    fn dual_twist(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCells
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve
        )
    }

    fn joint_leg_solve(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve
        )
    }

    fn beta_envelope(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope
        )
    }

    fn pair_normal_solve(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve
        )
    }

    fn active_normal_solve(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve
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

/// Chain network plus shared cross-tangent row-1 projection.
pub(super) fn author_chain_network_with_cross_tangent_ribbons(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(brep, slots, ChainNetworkMode::CurvesAndCrossTangentRibbons)
}

/// Chain network plus constructive cell re-authoring.
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

pub(super) fn author_chain_network_with_constructive_dual_twist_joint_leg_solve(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(
        brep,
        slots,
        ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve,
    )
}

pub(super) fn author_chain_network_with_constructive_dual_twist_beta_envelope(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(
        brep,
        slots,
        ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope,
    )
}

pub(super) fn author_chain_network_with_constructive_dual_twist_pair_normal_solve(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(
        brep,
        slots,
        ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve,
    )
}

pub(super) fn author_chain_network_with_constructive_dual_twist_active_normal_solve(
    brep: &mut BRep,
    slots: &[crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot],
) -> Result<String, CoreError> {
    author_chain_network_internal(
        brep,
        slots,
        ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve,
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
    let mut rotated_leg_count = 0usize;
    let mut joint_leg_vector_candidate_count = 0usize;
    let mut max_joint_leg_vector_candidate_delta = 0.0_f64;
    let mut joint_endpoint_equation_candidate_count = 0usize;
    let mut max_joint_endpoint_residual_before = 0.0_f64;
    let mut joint_plane_constraint_count = 0usize;
    let mut max_joint_plane_residual_before = 0.0_f64;
    let mut max_joint_plane_residual_after = 0.0_f64;
    let mut joint_leg_solve_declined_count = 0usize;
    if mode.constructive() {
        let (mut class_targets, thin_count) = leg_class_targets(&authored_sub_chains);
        leg_class_thin_count += thin_count;
        let (direction_targets, rotated) =
            leg_direction_targets(&authored_sub_chains, &node_normals);
        rotated_leg_count = rotated;
        for (chain_index, rows) in direction_targets {
            class_targets.entry(chain_index).or_default().extend(rows);
        }
        let (_vector_targets, vector_count, max_vector_delta) =
            leg_vector_targets(&authored_sub_chains, &chain_trims, &patches);
        joint_leg_vector_candidate_count = vector_count;
        max_joint_leg_vector_candidate_delta = max_vector_delta;
        if mode.joint_leg_solve() {
            let equations = joint_endpoint_leg_equations(&authored_sub_chains, &chain_trims);
            let plane_constraints =
                joint_endpoint_plane_constraints(&authored_sub_chains, &chain_trims, &patches);
            match solve_joint_chain_network(
                &mut patches,
                &discovery.sub_chains,
                &authored_sub_chains,
                &node_positions,
                &node_normals,
                &class_targets,
                &[],
                &plane_constraints,
            ) {
                Ok(outcome) => {
                    joint_endpoint_equation_candidate_count = equations.len();
                    max_joint_endpoint_residual_before = equations
                        .iter()
                        .map(|equation| equation.residual_before)
                        .fold(0.0_f64, f64::max);
                    joint_plane_constraint_count = outcome.plane_constraint_count;
                    max_joint_plane_residual_before = outcome.max_plane_residual_before;
                    max_joint_plane_residual_after = outcome.max_plane_residual_after;
                    leg_class_chain_count = class_targets.len();
                    authored_sub_chains = outcome.authored_sub_chains;
                }
                Err(_) => joint_leg_solve_declined_count = 1,
            }
        } else {
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
    }
    let mut reauthor_detail = String::new();
    let mut dual_twist_pairs = HashMap::new();
    if mode.constructive() {
        let outcome = reauthor_cells_from_chain_network(
            &mut patches,
            &authored_sub_chains,
            &node_normals,
            mode.dual_twist(),
            mode.beta_envelope(),
            mode.pair_normal_solve(),
            mode.active_normal_solve(),
        );
        dual_twist_pairs = outcome.dual_twist_pairs.clone();
        reauthor_detail = chain_status::reauthor_detail(&outcome);
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
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve, true) => {
            "chain_constructive_dual_twist_joint_leg_solve_authored"
        }
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve, false) => {
            "chain_constructive_dual_twist_joint_leg_solve_declined"
        }
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope, true) => {
            "chain_constructive_dual_twist_beta_envelope_authored"
        }
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope, false) => {
            "chain_constructive_dual_twist_beta_envelope_declined"
        }
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve, true) => {
            "chain_constructive_dual_twist_pair_normal_solve_authored"
        }
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve, false) => {
            "chain_constructive_dual_twist_pair_normal_solve_declined"
        }
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve, true) => {
            "chain_constructive_dual_twist_active_normal_solve_authored"
        }
        (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve, false) => {
            "chain_constructive_dual_twist_active_normal_solve_declined"
        }
    };
    Ok(format!(
        "{label};sub_chain_count={};authored_sub_chain_count={authored_count};fit_declined_sub_chain_count={fit_declined_count};fit_decline_reason={fit_decline_reason};discovery_declined_chain_count={};discovery_decline_reason={};chained_face_count={};uncovered_face_count={};uncovered_faces=[{uncovered_face_list}];glued_trim_count={};skipped_trim_missing_pcurve_count={};skipped_trim_unrecognized_span_count={};glued_row_count={};assigned_row_count={};dropped_row_count={};dropped_rows=[{dropped_row_list}];multiply_assigned_row_count={multiply_assigned_row_count};cross_sub_chain_trim_count={cross_sub_chain_trim_count};projected_row_count={projected_row_count};damped_row_count={damped_row_count};skipped_row_count={skipped_row_count};ribbon_declined_sub_chain_count={ribbon_declined_count};min_end_spacing_ratio={min_end_spacing_ratio:.4e};leg_class_thin_count={leg_class_thin_count};leg_class_chain_count={leg_class_chain_count};rotated_leg_count={rotated_leg_count};joint_leg_vector_candidate_count={joint_leg_vector_candidate_count};max_joint_leg_vector_candidate_delta={max_joint_leg_vector_candidate_delta:.6e};joint_endpoint_equation_candidate_count={joint_endpoint_equation_candidate_count};max_joint_endpoint_residual_before={max_joint_endpoint_residual_before:.6e};joint_plane_constraint_count={joint_plane_constraint_count};max_joint_plane_residual_before={max_joint_plane_residual_before:.6e};max_joint_plane_residual_after={max_joint_plane_residual_after:.6e};joint_leg_solve_declined_count={joint_leg_solve_declined_count}{reauthor_detail}",
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

#[cfg(test)]
#[path = "kp20_polar_support_cell_chain_fit_tests.rs"]
mod tests;
