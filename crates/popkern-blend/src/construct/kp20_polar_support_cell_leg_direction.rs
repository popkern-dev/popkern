//! Per-vertex joint leg-direction solve.
//!
//! The residual G1 mechanism (PK-085 oos-decomposition probe):
//! at vertices where two incident chains OSCULATE (T-junction
//! stem vertices meet their through-chain at 180 degrees by
//! construction; sliver tips near 0), the crossing rows' frozen
//! corner legs are chain-collinear, the band derivation must
//! deviate from them, and the per-side deviations make the trim
//! sides disagree. Speed retargeting cannot help — the legs are
//! wrong in DIRECTION. This pass builds a per-vertex registry of
//! incident chain tangents and rotates each rotatable chain END
//! that osculates another incident direction AWAY from it,
//! inside the node tangent plane, magnitude preserved. Through-
//! chain tangents at mid-chain stations are direction-fixed
//! (their continuity owns them) and act as obstacles only.
//! Healthy leg directions let the bands interpolate their frozen
//! ends exactly, restoring the in-span product identity on both
//! trim sides.

use std::collections::HashMap;

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::VertexHandle;

use super::cell_reauthor::AuthoredSubChain;
use super::chain_solver::dbernstein9;

type ChainDerivativeTargets = HashMap<usize, Vec<(f64, Vec3, f64)>>;

/// Sine bound under which two incident directions count as
/// osculating (matches the band derivation's degenerate-corner
/// tolerance with headroom).
const OSCULATION_SIN_BOUND: f64 = 0.05;

/// Rotated separation target (sine): far enough that the leg
/// solve is well-conditioned, small enough to bound the shape
/// change at the corner.
const ROTATION_SIN_TARGET: f64 = 0.25;

/// Soft weight for the rotation derivative targets (the scale
/// measured optimal for station retargets).
const ROTATION_TARGET_WEIGHT: f64 = 5.0;

/// One incident chain direction at a vertex.
struct IncidentDirection {
    chain_index: usize,
    station: f64,
    direction: Vec3,
    speed: f64,
    /// Chain ends may rotate; mid-chain (through) stations are
    /// fixed obstacles.
    rotatable: bool,
}

/// Build per-vertex direction registries and emit derivative
/// retargets rotating osculating chain ends apart. Returns
/// per-chain `(station, target derivative, weight)` rows plus
/// the rotated-end count.
pub(super) fn leg_direction_targets(
    authored_sub_chains: &[AuthoredSubChain],
    node_normals: &HashMap<VertexHandle, Vec3>,
) -> (ChainDerivativeTargets, usize) {
    let mut vertices: HashMap<VertexHandle, Vec<IncidentDirection>> = HashMap::new();
    for (chain_index, authored) in authored_sub_chains.iter().enumerate() {
        for &(vertex, station) in &authored.nodes {
            let derivative = chain_derivative(&authored.chain_cps, station);
            let speed = derivative.dot(&derivative).sqrt();
            if speed <= 1.0e-300 {
                continue;
            }
            vertices.entry(vertex).or_default().push(IncidentDirection {
                chain_index,
                station,
                direction: derivative * (1.0 / speed),
                speed,
                rotatable: station <= 1.0e-6 || station >= 1.0 - 1.0e-6,
            });
        }
    }
    let mut targets: ChainDerivativeTargets = HashMap::new();
    let mut rotated = 0usize;
    let mut vertex_order: Vec<&VertexHandle> = vertices.keys().collect();
    vertex_order.sort_by_key(|vertex| format!("{vertex:?}"));
    for vertex in vertex_order {
        let incidents = &vertices[vertex];
        let Some(&normal) = node_normals.get(vertex) else {
            continue;
        };
        for (index, incident) in incidents.iter().enumerate() {
            if !incident.rotatable {
                continue;
            }
            // The nearest other incident direction (by sine).
            let mut worst: Option<(f64, Vec3)> = None;
            for (other_index, other) in incidents.iter().enumerate() {
                if other_index == index {
                    continue;
                }
                let sine = incident
                    .direction
                    .cross(&other.direction)
                    .dot(&incident.direction.cross(&other.direction))
                    .sqrt();
                if worst.is_none_or(|(s, _)| sine < s) {
                    worst = Some((sine, other.direction));
                }
            }
            let Some((sine, obstacle)) = worst else {
                continue;
            };
            if sine >= OSCULATION_SIN_BOUND {
                continue;
            }
            // Rotate inside the node tangent plane: keep the
            // obstacle-parallel component's sign, set the
            // transverse component to the separation target.
            let along = obstacle * incident.direction.dot(&obstacle);
            let across = incident.direction - along;
            let across_unit = match across.normalized() {
                Some(unit) => unit,
                // Exactly collinear: take the in-plane
                // perpendicular to the obstacle.
                None => match normal.cross(&obstacle).normalized() {
                    Some(unit) => unit,
                    None => continue,
                },
            };
            let along_sign = if incident.direction.dot(&obstacle) >= 0.0 {
                1.0
            } else {
                -1.0
            };
            let rotated_direction = obstacle
                * (along_sign * (1.0 - ROTATION_SIN_TARGET.powi(2)).sqrt())
                + across_unit * ROTATION_SIN_TARGET;
            rotated += 1;
            targets.entry(incident.chain_index).or_default().push((
                incident.station,
                rotated_direction * incident.speed,
                ROTATION_TARGET_WEIGHT,
            ));
        }
    }
    (targets, rotated)
}

fn chain_derivative(cps: &[Point3; 10], param: f64) -> Vec3 {
    let dbasis = dbernstein9(param);
    let mut derivative = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..10 {
        derivative = derivative + Vec3::new(cps[i].x, cps[i].y, cps[i].z) * dbasis[i];
    }
    derivative
}
