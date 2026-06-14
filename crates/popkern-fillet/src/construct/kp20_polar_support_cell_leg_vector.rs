//! Per-vertex leg-vector consistency targets for constructive
//! support cells.
//!
//! The row-1 derivation keeps each boundary row's endpoint
//! controls frozen because those controls belong to the
//! perpendicular boundary chains. At T-junctions and shared trim
//! endpoints, that means one side can carry a CP1 leg that no
//! longer agrees with the other side's row-1 polynomial at the
//! same vertex. This registry compares the cross-derivative
//! vectors at T-junction shared-trim endpoints and emits soft
//! derivative targets for the chain-authored perpendicular row
//! that owns the endpoint CP1 control. Endpoint-endpoint joins
//! are left to the existing magnitude/direction registries; this
//! pass is for the asymmetric case where one side's endpoint
//! must agree with the other side's row-1 polynomial.

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::FaceHandle;

use super::BoundarySide;
use super::cell_reauthor::{AuthoredSubChain, second_row_controls};
use super::chain_author::boundary_controls;
use super::chain_discovery::ChainTrim;
use super::chain_fit::ChainRowKey;
use super::chain_solver::dbernstein9;
use super::row_fit::bernstein9;

const ROW_ENDPOINT_TOL: f64 = 1.0e-6;
const LEG_VECTOR_TARGET_WEIGHT: f64 = 5.0;

#[derive(Debug, Clone, Copy)]
struct AuthoredRow {
    chain_index: usize,
    start: f64,
    end: f64,
}

#[derive(Debug, Clone, Copy)]
struct EndpointOwner {
    key: ChainRowKey,
    row_param: f64,
    sign: f64,
}

#[derive(Debug, Clone, Copy)]
struct RetargetableEndpoint {
    chain_index: usize,
    station: f64,
    scale: f64,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct JointEndpointLeg {
    pub(super) chain_index: usize,
    pub(super) station: f64,
    pub(super) derivative_scale: f64,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct JointEndpointLegEquation {
    pub(super) first: JointEndpointLeg,
    pub(super) second: JointEndpointLeg,
    pub(super) second_orientation: f64,
    pub(super) residual_before: f64,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct JointEndpointPlaneConstraint {
    pub(super) endpoint: JointEndpointLeg,
    pub(super) plane_normal: Vec3,
    pub(super) residual_before: f64,
}

/// Build soft derivative targets that make endpoint CP1 legs at
/// T-junctions agree with the opposing trim side's interior
/// cross-derivative vector, up to orientation. Returns per-chain
/// `(station, target derivative, weight)` rows, the emitted
/// endpoint-target count, and the largest requested leg-vector
/// displacement.
pub(super) fn leg_vector_targets(
    authored_sub_chains: &[AuthoredSubChain],
    chain_trims: &[ChainTrim<ChainRowKey, popkern_core::topology::VertexHandle>],
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
) -> (HashMap<usize, Vec<(f64, Vec3, f64)>>, usize, f64) {
    let authored_rows = authored_row_lookup(authored_sub_chains);
    let mut targets: HashMap<usize, Vec<(f64, Vec3, f64)>> = HashMap::new();
    let mut target_count = 0usize;
    let mut max_leg_delta = 0.0_f64;
    for trim in chain_trims {
        for (row_a, param_a, row_b, param_b) in [
            (trim.row_a, trim.interval_a.0, trim.row_b, trim.interval_b.1),
            (trim.row_a, trim.interval_a.1, trim.row_b, trim.interval_b.0),
        ] {
            let Some(vector_a) = row_cross_derivative(patches, row_a, param_a) else {
                continue;
            };
            let Some(vector_b) = row_cross_derivative(patches, row_b, param_b) else {
                continue;
            };
            match (row_is_endpoint(param_a), row_is_endpoint(param_b)) {
                (true, false) => {
                    if let Some(endpoint_a) = retargetable_endpoint(&authored_rows, row_a, param_a)
                    {
                        let target_a = aligned_to(vector_a, vector_b);
                        push_target(
                            &mut targets,
                            endpoint_a,
                            target_a,
                            vector_a,
                            &mut target_count,
                            &mut max_leg_delta,
                        );
                    }
                }
                (false, true) => {
                    if let Some(endpoint_b) = retargetable_endpoint(&authored_rows, row_b, param_b)
                    {
                        let target_b = aligned_to(vector_b, vector_a);
                        push_target(
                            &mut targets,
                            endpoint_b,
                            target_b,
                            vector_b,
                            &mut target_count,
                            &mut max_leg_delta,
                        );
                    }
                }
                _ => {}
            }
        }
    }
    (targets, target_count, max_leg_delta)
}

/// Build coupled endpoint equations where both cross-derivative
/// legs are owned by chain-curve endpoint derivatives. Unlike
/// [`leg_vector_targets`], these rows do not freeze either side:
/// the global chain solve links both derivative unknowns in one
/// system and can move the pair coherently.
pub(super) fn joint_endpoint_leg_equations(
    authored_sub_chains: &[AuthoredSubChain],
    chain_trims: &[ChainTrim<ChainRowKey, popkern_core::topology::VertexHandle>],
) -> Vec<JointEndpointLegEquation> {
    let authored_rows = authored_row_lookup(authored_sub_chains);
    let mut equations = Vec::new();
    let mut seen = HashSet::new();
    for trim in chain_trims {
        for (row_a, param_a, row_b, param_b) in [
            (trim.row_a, trim.interval_a.0, trim.row_b, trim.interval_b.1),
            (trim.row_a, trim.interval_a.1, trim.row_b, trim.interval_b.0),
        ] {
            if !row_is_endpoint(param_a) || !row_is_endpoint(param_b) {
                continue;
            }
            let Some(first) = retargetable_endpoint(&authored_rows, row_a, param_a) else {
                continue;
            };
            let Some(second) = retargetable_endpoint(&authored_rows, row_b, param_b) else {
                continue;
            };
            if same_endpoint(first, second) {
                continue;
            }
            let key = endpoint_pair_key(first, second);
            if !seen.insert(key) {
                continue;
            }
            let first_leg = endpoint_leg(authored_sub_chains, first);
            let second_leg = endpoint_leg(authored_sub_chains, second);
            let orientation = if first_leg.dot(&second_leg) < 0.0 {
                -1.0
            } else {
                1.0
            };
            let residual = first_leg - second_leg * orientation;
            equations.push(JointEndpointLegEquation {
                first: JointEndpointLeg {
                    chain_index: first.chain_index,
                    station: first.station,
                    derivative_scale: first.scale,
                },
                second: JointEndpointLeg {
                    chain_index: second.chain_index,
                    station: second.station,
                    derivative_scale: second.scale,
                },
                second_orientation: orientation,
                residual_before: residual.dot(&residual).sqrt(),
            });
        }
    }
    equations
}

/// Build asymmetric endpoint-to-interior plane constraints for
/// T-junction trim endpoints. The endpoint side contributes the
/// chain-owned transverse leg unknown; the interior side
/// contributes the tangent plane measured from its current row-1
/// field. This enforces the G1 invariant that the endpoint leg
/// lies in the peer tangent plane without forcing vector
/// equality to a frozen peer derivative.
pub(super) fn joint_endpoint_plane_constraints(
    authored_sub_chains: &[AuthoredSubChain],
    chain_trims: &[ChainTrim<ChainRowKey, popkern_core::topology::VertexHandle>],
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
) -> Vec<JointEndpointPlaneConstraint> {
    let authored_rows = authored_row_lookup(authored_sub_chains);
    let mut constraints = Vec::new();
    let mut seen = HashSet::new();
    for trim in chain_trims {
        for (row_a, param_a, row_b, param_b) in [
            (trim.row_a, trim.interval_a.0, trim.row_b, trim.interval_b.1),
            (trim.row_a, trim.interval_a.1, trim.row_b, trim.interval_b.0),
        ] {
            match (row_is_endpoint(param_a), row_is_endpoint(param_b)) {
                (true, false) => push_plane_constraint(
                    &authored_rows,
                    authored_sub_chains,
                    patches,
                    row_a,
                    param_a,
                    row_b,
                    param_b,
                    &mut constraints,
                    &mut seen,
                ),
                (false, true) => push_plane_constraint(
                    &authored_rows,
                    authored_sub_chains,
                    patches,
                    row_b,
                    param_b,
                    row_a,
                    param_a,
                    &mut constraints,
                    &mut seen,
                ),
                _ => {}
            }
        }
    }
    constraints
}

fn authored_row_lookup(
    authored_sub_chains: &[AuthoredSubChain],
) -> HashMap<ChainRowKey, AuthoredRow> {
    let mut rows = HashMap::new();
    for (chain_index, authored) in authored_sub_chains.iter().enumerate() {
        for &(key, (start, end)) in &authored.rows {
            rows.entry(key).or_insert(AuthoredRow {
                chain_index,
                start,
                end,
            });
        }
    }
    rows
}

fn row_cross_derivative(
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
    key: ChainRowKey,
    row_param: f64,
) -> Option<Vec3> {
    let patch = patches.get(&key.face)?;
    let row0 = boundary_controls(patch, key.side);
    let row1 = second_row_controls(patch, key.side);
    let basis = bernstein9(row_param.clamp(0.0, 1.0));
    let mut vector = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..10 {
        vector = vector + (row1[i] - row0[i]) * basis[i];
    }
    Some(vector * 9.0)
}

fn row_tangent(
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
    key: ChainRowKey,
    row_param: f64,
) -> Option<Vec3> {
    let patch = patches.get(&key.face)?;
    let row0 = boundary_controls(patch, key.side);
    Some(chain_derivative(&row0, row_param.clamp(0.0, 1.0)))
}

fn push_plane_constraint(
    authored_rows: &HashMap<ChainRowKey, AuthoredRow>,
    authored_sub_chains: &[AuthoredSubChain],
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
    endpoint_row: ChainRowKey,
    endpoint_param: f64,
    interior_row: ChainRowKey,
    interior_param: f64,
    constraints: &mut Vec<JointEndpointPlaneConstraint>,
    seen: &mut HashSet<String>,
) {
    let Some(endpoint) = retargetable_endpoint(authored_rows, endpoint_row, endpoint_param) else {
        return;
    };
    let Some(tangent) = row_tangent(patches, interior_row, interior_param) else {
        return;
    };
    let Some(cross) = row_cross_derivative(patches, interior_row, interior_param) else {
        return;
    };
    let Some(plane_normal) = tangent.cross(&cross).normalized() else {
        return;
    };
    let key = format!(
        "{}:{:.12}:{:.12}:{:?}:{:?}:{:.12}",
        endpoint.chain_index,
        endpoint.station,
        endpoint.scale,
        interior_row.face,
        interior_row.side,
        interior_param
    );
    if !seen.insert(key) {
        return;
    }
    let leg = endpoint_leg(authored_sub_chains, endpoint);
    constraints.push(JointEndpointPlaneConstraint {
        endpoint: JointEndpointLeg {
            chain_index: endpoint.chain_index,
            station: endpoint.station,
            derivative_scale: endpoint.scale,
        },
        plane_normal,
        residual_before: leg.dot(&plane_normal).abs(),
    });
}

fn retargetable_endpoint(
    authored_rows: &HashMap<ChainRowKey, AuthoredRow>,
    row: ChainRowKey,
    row_param: f64,
) -> Option<RetargetableEndpoint> {
    let owner = endpoint_owner(row, row_param)?;
    let authored = authored_rows.get(&owner.key)?;
    let station = if owner.row_param <= ROW_ENDPOINT_TOL {
        authored.start
    } else {
        authored.end
    };
    let width = authored.end - authored.start;
    if width.abs() <= 1.0e-12 {
        return None;
    }
    Some(RetargetableEndpoint {
        chain_index: authored.chain_index,
        station,
        scale: owner.sign / width,
    })
}

fn same_endpoint(first: RetargetableEndpoint, second: RetargetableEndpoint) -> bool {
    first.chain_index == second.chain_index && (first.station - second.station).abs() <= 1.0e-9
}

fn endpoint_pair_key(first: RetargetableEndpoint, second: RetargetableEndpoint) -> String {
    let a = format!(
        "{}:{:.12}:{:.12}",
        first.chain_index, first.station, first.scale
    );
    let b = format!(
        "{}:{:.12}:{:.12}",
        second.chain_index, second.station, second.scale
    );
    if a <= b {
        format!("{a}|{b}")
    } else {
        format!("{b}|{a}")
    }
}

fn endpoint_leg(authored_sub_chains: &[AuthoredSubChain], endpoint: RetargetableEndpoint) -> Vec3 {
    let derivative = chain_derivative(
        &authored_sub_chains[endpoint.chain_index].chain_cps,
        endpoint.station,
    );
    derivative * (1.0 / endpoint.scale)
}

pub(super) fn chain_derivative(chain_cps: &[popkern_core::geometry::Point3; 10], t: f64) -> Vec3 {
    let basis = dbernstein9(t);
    let mut derivative = Vec3::new(0.0, 0.0, 0.0);
    for (index, weight) in basis.iter().enumerate() {
        derivative = Vec3::new(
            derivative.x + chain_cps[index].x * weight,
            derivative.y + chain_cps[index].y * weight,
            derivative.z + chain_cps[index].z * weight,
        );
    }
    derivative
}

fn endpoint_owner(row: ChainRowKey, row_param: f64) -> Option<EndpointOwner> {
    if row_is_start(row_param) {
        return Some(endpoint_owner_for_corner(row, 0));
    }
    if row_is_end(row_param) {
        return Some(endpoint_owner_for_corner(row, 9));
    }
    None
}

fn row_is_endpoint(row_param: f64) -> bool {
    row_is_start(row_param) || row_is_end(row_param)
}

fn row_is_start(row_param: f64) -> bool {
    row_param.abs() <= ROW_ENDPOINT_TOL
}

fn row_is_end(row_param: f64) -> bool {
    (row_param - 1.0).abs() <= ROW_ENDPOINT_TOL
}

fn endpoint_owner_for_corner(row: ChainRowKey, corner_k: usize) -> EndpointOwner {
    let (side, row_param, sign) = match (row.side, corner_k) {
        (BoundarySide::UMin, 0) => (BoundarySide::VMin, 0.0, 1.0),
        (BoundarySide::UMin, _) => (BoundarySide::VMax, 0.0, 1.0),
        (BoundarySide::UMax, 0) => (BoundarySide::VMin, 1.0, -1.0),
        (BoundarySide::UMax, _) => (BoundarySide::VMax, 1.0, -1.0),
        (BoundarySide::VMin, 0) => (BoundarySide::UMin, 0.0, 1.0),
        (BoundarySide::VMin, _) => (BoundarySide::UMax, 0.0, 1.0),
        (BoundarySide::VMax, 0) => (BoundarySide::UMin, 1.0, -1.0),
        (BoundarySide::VMax, _) => (BoundarySide::UMax, 1.0, -1.0),
    };
    EndpointOwner {
        key: ChainRowKey {
            face: row.face,
            side,
        },
        row_param,
        sign,
    }
}

fn aligned_to(current: Vec3, authority: Vec3) -> Vec3 {
    if current.dot(&authority) < 0.0 {
        authority * -1.0
    } else {
        authority
    }
}

fn push_target(
    targets: &mut HashMap<usize, Vec<(f64, Vec3, f64)>>,
    endpoint: RetargetableEndpoint,
    target_leg: Vec3,
    current_leg: Vec3,
    target_count: &mut usize,
    max_leg_delta: &mut f64,
) {
    let target_derivative = target_leg * endpoint.scale;
    targets.entry(endpoint.chain_index).or_default().push((
        endpoint.station,
        target_derivative,
        LEG_VECTOR_TARGET_WEIGHT,
    ));
    *target_count += 1;
    let delta = (target_leg - current_leg)
        .dot(&(target_leg - current_leg))
        .sqrt();
    *max_leg_delta = max_leg_delta.max(delta);
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_leg_vector_tests.rs"]
mod tests;
