//! Per-vertex leg-class registry for the constructive
//! support-cell pass.
//!
//! Every row corner leg at a vertex inherits its magnitude from
//! its chain's CP spacing at that station, and the chains' free
//! speed profiles leave some legs an order of magnitude thinner
//! than their siblings — those corners' bands turn tangential
//! and their trim normals are noise. The healthiest incident leg
//! defines the vertex's scale (a SHARED registry, like the
//! pinned node positions and normals); stations whose legs fall
//! below the class floor get direction-preserving derivative
//! retargets raised to the class scale. SINGLE round by design:
//! iterating re-measures against a raised class scale (a moving
//! target); the thin set GROWS round over round and occasionally
//! destabilizes the whole candidate — measured divergent.

use std::collections::HashMap;

use popkern_core::geometry::Vec3;
use popkern_core::topology::VertexHandle;

use super::cell_reauthor::AuthoredSubChain;
use super::chain_solver::dbernstein9;

/// A corner leg below this fraction of its vertex's class scale
/// (the healthiest incident leg) is raised to the class target.
const LEG_CLASS_FLOOR_FRACTION: f64 = 0.45;
/// Demanded leg fraction of the class scale after the refit.
const LEG_CLASS_TARGET_FRACTION: f64 = 0.6;
/// Soft weight for leg-class derivative targets (the gentle
/// scale measured optimal for station retargets).
const LEG_CLASS_TARGET_WEIGHT: f64 = 5.0;

/// Build the per-vertex leg registry from the authored chains
/// and return per-chain derivative retargets for the stations
/// whose legs fall below the class floor, plus the thin-leg
/// count. Iteration is fully deterministic (sorted vertices,
/// chains, and stations): the refit order changes the float
/// accumulation path enough to flip downstream measured-gate
/// decisions run to run otherwise.
pub(super) fn leg_class_targets(
    authored_sub_chains: &[AuthoredSubChain],
) -> (HashMap<usize, Vec<(f64, Vec3, f64)>>, usize) {
    let mut thin_count = 0usize;
    let mut vertex_legs: HashMap<VertexHandle, Vec<(usize, f64, f64, f64)>> = HashMap::new();
    for (chain_index, authored) in authored_sub_chains.iter().enumerate() {
        for &(_, (start, end)) in &authored.rows {
            let width = (end - start).abs().max(1.0e-12);
            for station in [start, end] {
                let Some(&(vertex, _)) = authored
                    .nodes
                    .iter()
                    .find(|&&(_, param)| (param - station).abs() <= 1.0e-6)
                else {
                    continue;
                };
                let speed = chain_speed(&authored.chain_cps, station);
                let leg = speed.dot(&speed).sqrt() * width / 9.0;
                vertex_legs
                    .entry(vertex)
                    .or_default()
                    .push((chain_index, station, width, leg));
            }
        }
    }
    let mut targets: HashMap<usize, Vec<(f64, Vec3, f64)>> = HashMap::new();
    let mut vertex_order: Vec<&VertexHandle> = vertex_legs.keys().collect();
    vertex_order.sort_by_key(|vertex| format!("{vertex:?}"));
    for vertex in vertex_order {
        let legs = &vertex_legs[vertex];
        let scale = legs
            .iter()
            .fold(0.0_f64, |acc, &(_, _, _, leg)| acc.max(leg));
        if scale <= 1.0e-300 {
            continue;
        }
        for &(chain_index, station, width, leg) in legs {
            if leg >= LEG_CLASS_FLOOR_FRACTION * scale || leg <= 1.0e-300 {
                continue;
            }
            thin_count += 1;
            let target_speed = LEG_CLASS_TARGET_FRACTION * scale * 9.0 / width;
            let derivative = chain_speed(&authored_sub_chains[chain_index].chain_cps, station);
            let speed = derivative.dot(&derivative).sqrt();
            if speed <= 1.0e-300 {
                continue;
            }
            targets.entry(chain_index).or_default().push((
                station,
                derivative * (target_speed / speed),
                LEG_CLASS_TARGET_WEIGHT,
            ));
        }
    }
    (targets, thin_count)
}

fn chain_speed(cps: &[popkern_core::geometry::Point3; 10], param: f64) -> Vec3 {
    let dbasis = dbernstein9(param);
    let mut derivative = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..10 {
        derivative = derivative + Vec3::new(cps[i].x, cps[i].y, cps[i].z) * dbasis[i];
    }
    derivative
}
