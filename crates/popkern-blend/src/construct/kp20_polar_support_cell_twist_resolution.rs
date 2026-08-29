//! Twist-control resolution for the constructive support-cell
//! pass: a grid corner's diagonal neighbour is demanded by two
//! boundaries of the same patch. The plain bi-nonic route
//! resolves the conflict by inverse-square field-magnitude
//! weighting (the thin boundary is protected); the dual-twist
//! route gives each boundary its own Chiyokura-Kimura candidate
//! slot, so the conflict ceases to exist.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::FaceHandle;

use super::cell_reauthor::CellReauthorOutcome;

/// Map a grid index to the dual-twist corner slot of the
/// `RectangularGregoryBiNonicDualTwist` convention (CCW from the
/// parameter origin).
pub(super) fn dual_twist_corner_index(i: usize, j: usize) -> Option<usize> {
    match (i, j) {
        (1, 1) => Some(0),
        (1, 8) => Some(1),
        (8, 8) => Some(2),
        (8, 1) => Some(3),
        _ => None,
    }
}

/// Dual-twist resolution of a corner control: each demanding
/// boundary keeps its own candidate (u-constant boundaries own
/// `P_u`, v-constant own `P_v`); a missing side copies the
/// other, degenerating to the plain control. The grid slot
/// receives the pair midpoint — the corner-limit convention of
/// the rational blend, used by the interior Coons refill; the
/// surface's boundary cross-derivatives read the candidates.
pub(super) fn resolve_dual_twist_corner(
    face: FaceHandle,
    i: usize,
    j: usize,
    demands: &[TwistDemand],
    outcome: &mut CellReauthorOutcome,
) -> Point3 {
    let corner = dual_twist_corner_index(i, j).unwrap_or(0);
    let mean_of = |u_constant: bool| -> Option<Point3> {
        let mut sum = Point3::new(0.0, 0.0, 0.0);
        let mut count = 0usize;
        for demand in demands {
            if demand.u_constant_boundary == u_constant {
                sum = Point3::new(
                    sum.x + demand.position.x,
                    sum.y + demand.position.y,
                    sum.z + demand.position.z,
                );
                count += 1;
            }
        }
        if count == 0 {
            return None;
        }
        let inverse = 1.0 / count as f64;
        Some(Point3::new(
            sum.x * inverse,
            sum.y * inverse,
            sum.z * inverse,
        ))
    };
    let u_count = demands.iter().filter(|d| d.u_constant_boundary).count();
    let v_count = demands.len() - u_count;
    if u_count > 1 || v_count > 1 {
        outcome.dual_multi_demand_side_count += 1;
    }
    let (u_side, v_side) = (mean_of(true), mean_of(false));
    let p_u = u_side.or(v_side).unwrap_or(Point3::new(0.0, 0.0, 0.0));
    let p_v = v_side.or(u_side).unwrap_or(p_u);
    let separation = p_u.distance_to(&p_v);
    if separation > 0.0 {
        outcome.dual_twist_pair_count += 1;
        outcome.max_dual_twist_pair_separation =
            outcome.max_dual_twist_pair_separation.max(separation);
    }
    outcome.dual_twist_pairs.entry(face).or_default()[corner] = Some((p_u, p_v));
    Point3::new(
        0.5 * (p_u.x + p_v.x),
        0.5 * (p_u.y + p_v.y),
        0.5 * (p_u.z + p_v.z),
    )
}

/// One boundary's demand for a shared band control: the desired
/// position, the corner tangent-plane normal (twist controls
/// only), the demanding boundary's local field magnitude, and
/// whether the demanding boundary is u-constant (`UMin`/`UMax`)
/// — the dual-twist candidate slot it owns.
#[derive(Debug, Clone, Copy)]
pub(super) struct TwistDemand {
    pub(super) position: Point3,
    pub(super) plane_normal: Option<Vec3>,
    pub(super) field_magnitude: f64,
    pub(super) u_constant_boundary: bool,
}

/// Resolve possibly-conflicting demands for one band control.
/// Conflicts take the inverse-square field-magnitude weighted
/// mean — the closed-form minimizer of
/// `Σ_b (|q − demand_b| / |Δ_b|)²` — so the displacement lands
/// on the boundary that can absorb it relative to its own field
/// strength. Field magnitudes are clamped below at a fraction of
/// the largest demand magnitude to keep the weights finite.
pub(super) fn resolve_twist_demands(
    demands: &[TwistDemand],
    outcome: &mut CellReauthorOutcome,
) -> Point3 {
    let count = demands.len();
    if count == 1 {
        return demands[0].position;
    }
    outcome.twist_conflict_count += 1;
    for first in 0..count {
        for second in first + 1..count {
            let difference = demands[second].position - demands[first].position;
            outcome.max_twist_conflict_distance = outcome
                .max_twist_conflict_distance
                .max(difference.dot(&difference).sqrt());
            for normal in [demands[first].plane_normal, demands[second].plane_normal]
                .into_iter()
                .flatten()
            {
                outcome.max_twist_conflict_normal_component = outcome
                    .max_twist_conflict_normal_component
                    .max(difference.dot(&normal).abs());
            }
        }
    }
    let magnitude_floor = demands
        .iter()
        .map(|demand| demand.field_magnitude)
        .fold(0.0_f64, f64::max)
        * 1.0e-3;
    let mut weight_sum = 0.0_f64;
    let mut weighted = Point3::new(0.0, 0.0, 0.0);
    for demand in demands {
        let magnitude = demand.field_magnitude.max(magnitude_floor).max(1.0e-300);
        let weight = 1.0 / (magnitude * magnitude);
        weight_sum += weight;
        weighted = Point3::new(
            weighted.x + demand.position.x * weight,
            weighted.y + demand.position.y * weight,
            weighted.z + demand.position.z * weight,
        );
    }
    let resolved = Point3::new(
        weighted.x / weight_sum,
        weighted.y / weight_sum,
        weighted.z / weight_sum,
    );
    for demand in demands {
        let displacement = resolved.distance_to(&demand.position);
        let magnitude = demand.field_magnitude.max(magnitude_floor).max(1.0e-300);
        outcome.max_relative_twist_displacement = outcome
            .max_relative_twist_displacement
            .max(displacement / magnitude);
    }
    resolved
}
