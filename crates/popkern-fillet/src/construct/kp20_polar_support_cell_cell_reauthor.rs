//! Constructive support-cell re-authoring from the chain network.
//!
//! Instead of projecting the old (mutually disagreeing) sheet
//! geometry toward shared planes, every cell's near-boundary
//! band is DERIVED from shared objects: per sub-chain, the
//! canonical curve `C` and a degree-8 transverse ribbon field
//! `R`. Each member row's cross-derivative becomes
//! `D(t) = α(t)·C'(t) + β(t)·R(t)` with LINEAR `α, β` solved
//! from the exact end legs (the transverse chain curves' first
//! control legs) — a degree-9 polynomial, exactly representable
//! in row 1, lying in `span(C'(s), R(s))` pointwise. Both sides
//! of every trim therefore share a tangent plane by
//! construction, with no fit residual. Corner solves are exact
//! because the legs, `C'`, and `R` all lie in the pinned node
//! plane (the chain curves are tangent-plane-constrained there).
//!
//! Twist controls (a grid corner's diagonal neighbour) are
//! demanded by two boundaries; the writes are averaged and
//! counted — the one approximation of the construction. Cell
//! interiors are refilled by a discrete Coons blend of the four
//! derived bands.

use std::collections::HashMap;

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::BoundarySide;
use super::band_field::author_ribbon_field;
use super::band_predictor::{derived_span_residual, predict_worst_trim_g1};
use super::chain_fit::{ChainRowKey, boundary_controls};
use super::chain_ribbon_field::{ribbon_field_segment, row_derivative_coefficients};

/// Determinant floor for the 2x2 end-leg solves (relative to the
/// basis vectors' scale).
const LEG_SOLVE_DET_TOL: f64 = 1.0e-10;

/// One authored sub-chain, retained by the curve-authoring pass
/// for the constructive re-author: the canonical curve, every
/// member row's local interval, and the node params/identities
/// (all in the sub-chain's local `[0, 1]` domain).
#[derive(Debug, Clone)]
pub(super) struct AuthoredSubChain {
    pub(super) chain_cps: [Point3; 10],
    pub(super) rows: Vec<(ChainRowKey, (f64, f64))>,
    pub(super) nodes: Vec<(VertexHandle, f64)>,
}

/// Outcome counts of the constructive re-author.
#[derive(Debug, Clone, Default)]
pub(super) struct CellReauthorOutcome {
    pub(super) derived_row_count: usize,
    pub(super) leg_solve_declined_count: usize,
    pub(super) field_declined_sub_chain_count: usize,
    pub(super) twist_conflict_count: usize,
    /// Largest distance between two boundaries' demands for one
    /// twist control. Averaging perturbs each boundary's derived
    /// field by `9·B₁(t)·conflict/2`, so this magnitude bounds
    /// the construction's G1 residual near corners — the bi-9
    /// twist-compatibility limit.
    pub(super) max_twist_conflict_distance: f64,
    /// Construction self-check: worst out-of-span component of
    /// the derived cross-derivative fields, relative to the
    /// field magnitude, sampled BEFORE the twist writes are
    /// averaged. The derivation is a polynomial identity in
    /// `span(C', R)`, so anything above roundoff here is an
    /// implementation defect, not a twist trade-off.
    pub(super) max_derived_span_residual: f64,
    /// Legs rotated into their node planes by the repair
    /// pre-pass, and the worst out-of-plane ratio found before
    /// repair. Chained rows' legs are already in-plane (their
    /// chains are tangent-constrained), so repairs discipline
    /// the unchained rim rows' transverse controls.
    pub(super) repaired_leg_count: usize,
    pub(super) max_leg_out_of_plane_ratio: f64,
    /// Worst NORMAL component of a twist-demand disagreement —
    /// the only part that damages G1 (in-plane disagreement
    /// only redistributes the tangential field).
    pub(super) max_twist_conflict_normal_component: f64,
    /// Worst displacement of a resolved twist from one of its
    /// demands, RELATIVE to that boundary's own field magnitude
    /// — the quantity that perturbs measured G1. The
    /// magnitude-aware resolution keeps it bounded where the
    /// unweighted average let it reach ~8x on thin rows.
    pub(super) max_relative_twist_displacement: f64,
    pub(super) refilled_cell_count: usize,
    /// The construction's own post-hoc G1 prediction: worst
    /// `1 − |n̂_a · n̂_b|` over every glued trim, sampled from the
    /// final control nets exactly the way the gate measures —
    /// splits construction-internal damage from later pipeline
    /// damage decisively.
    pub(super) predicted_worst_trim_g1: f64,
    /// Anatomy of the worst predicted sample (sides' own-span
    /// residuals, field values, intervals) for localization.
    pub(super) predicted_worst_detail: Option<String>,
    /// End legs whose transverse component was floor-clamped to
    /// prevent a mid-row parametric degeneracy.
    pub(super) clamped_leg_beta_count: usize,
}

/// Derive every chained row's row-1 band from the shared
/// `(C, R)` fields, then refill the touched cells' interiors.
pub(super) fn reauthor_cells_from_chain_network(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub_chains: &[AuthoredSubChain],
    node_normals: &HashMap<VertexHandle, Vec3>,
) -> CellReauthorOutcome {
    let mut outcome = CellReauthorOutcome::default();
    let chained_rows: std::collections::HashSet<ChainRowKey> = sub_chains
        .iter()
        .flat_map(|sub| sub.rows.iter().map(|&(key, _)| key))
        .collect();
    // Collect every derived control write; twist controls are
    // demanded by two boundaries of the same face and averaged,
    // and the demand spread is retained for observability.
    let mut writes: HashMap<(FaceHandle, usize, usize), Vec<TwistDemand>> = HashMap::new();
    for sub in sub_chains {
        let Some(field) = author_ribbon_field(sub, node_normals, patches, &chained_rows) else {
            outcome.field_declined_sub_chain_count += 1;
            continue;
        };
        for &(key, (start, end)) in &sub.rows {
            let Some(patch) = patches.get(&key.face) else {
                continue;
            };
            let row0 = boundary_controls(patch, key.side);
            let field_segment = ribbon_field_segment(&field, start, end);
            // Leg repair: the band derivation closes exactly only
            // if its end legs lie in the LOCAL derivation span
            // span(C'(end), R(end)). Chained transverse rows'
            // legs are already in-plane at nodes (tangent
            // constraints), so the rotation is a no-op for them;
            // unchained rim rows' junk legs get rotated in
            // (magnitude- and in-plane-direction-preserving) and
            // the transverse control is written back. Defined at
            // chain extremes too, where no node exists.
            repair_row_end_legs(
                patches,
                key,
                &row0,
                &field_segment,
                &chained_rows,
                &mut outcome,
            );
            let Some(patch) = patches.get(&key.face) else {
                continue;
            };
            let row1 = second_row_controls(patch, key.side);
            match derive_row1(
                &row0,
                &row1,
                &field_segment,
                &mut outcome.clamped_leg_beta_count,
            ) {
                Some(derived) => {
                    outcome.max_derived_span_residual = outcome
                        .max_derived_span_residual
                        .max(derived_span_residual(&row0, &derived, &field_segment));
                    let derivative = row_derivative_coefficients(&row0);
                    let corner_normal_start = derivative[0].cross(&field_segment[0]).normalized();
                    let corner_normal_end = derivative[8].cross(&field_segment[8]).normalized();
                    for (k, point) in derived.iter().enumerate().take(9).skip(1) {
                        let (i, j) = second_row_grid_index(key.side, k);
                        let plane_normal = match k {
                            1 => corner_normal_start,
                            8 => corner_normal_end,
                            _ => None,
                        };
                        let band_offset = *point - row0[k];
                        writes
                            .entry((key.face, i, j))
                            .or_default()
                            .push(TwistDemand {
                                position: *point,
                                plane_normal,
                                field_magnitude: band_offset.dot(&band_offset).sqrt(),
                            });
                    }
                    outcome.derived_row_count += 1;
                }
                None => outcome.leg_solve_declined_count += 1,
            }
        }
    }
    let mut touched_faces = Vec::new();
    for (&(face, i, j), demands) in &writes {
        let resolved = resolve_twist_demands(demands, &mut outcome);
        if let Some(patch) = patches.get_mut(&face) {
            patch.grid[i][j] = resolved;
            if !touched_faces.contains(&face) {
                touched_faces.push(face);
            }
        }
    }
    for face in touched_faces {
        if let Some(patch) = patches.get_mut(&face) {
            refill_interior_coons(patch);
            outcome.refilled_cell_count += 1;
        }
    }
    predict_worst_trim_g1(patches, sub_chains, node_normals, &mut outcome);
    outcome
}

/// One boundary's demand for a shared band control: the desired
/// position, the corner tangent-plane normal (twist controls
/// only), and the demanding boundary's local field magnitude.
#[derive(Debug, Clone, Copy)]
struct TwistDemand {
    position: Point3,
    plane_normal: Option<Vec3>,
    field_magnitude: f64,
}

/// Resolve possibly-conflicting demands for one band control.
/// Conflicts take the inverse-square field-magnitude weighted
/// mean — the closed-form minimizer of
/// `Σ_b (|q − demand_b| / |Δ_b|)²` — so the displacement lands
/// on the boundary that can absorb it relative to its own field
/// strength. Field magnitudes are clamped below at a fraction of
/// the largest demand magnitude to keep the weights finite.
fn resolve_twist_demands(demands: &[TwistDemand], outcome: &mut CellReauthorOutcome) -> Point3 {
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

/// Rotate a row's two end legs into the local derivation span
/// `span(C'(end), R(end))`, preserving magnitude and in-plane
/// direction, and write the repaired transverse controls back.
fn repair_row_end_legs(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    key: ChainRowKey,
    row0: &[Point3; 10],
    field_segment: &[Vec3; 9],
    chained_rows: &std::collections::HashSet<ChainRowKey>,
    outcome: &mut CellReauthorOutcome,
) {
    let derivative = row_derivative_coefficients(row0);
    for (corner_k, tangent, transverse) in [
        (0usize, derivative[0], field_segment[0]),
        (9, derivative[8], field_segment[8]),
    ] {
        // A chained transverse row's control is chain-authored
        // (exact G0); it must never move. Its leg is made
        // span-compatible on the field side instead (the ribbon
        // fit is constrained to contain it).
        if chained_rows.contains(&ChainRowKey {
            face: key.face,
            side: transverse_row_side(key.side, corner_k),
        }) {
            continue;
        }
        let Some(plane_normal) = tangent.cross(&transverse).normalized() else {
            continue;
        };
        let Some(patch) = patches.get(&key.face) else {
            continue;
        };
        let row1 = second_row_controls(patch, key.side);
        let leg = row1[corner_k] - row0[corner_k];
        let magnitude = leg.dot(&leg).sqrt();
        if magnitude <= 1.0e-300 {
            continue;
        }
        let out_of_plane = leg.dot(&plane_normal);
        let ratio = out_of_plane.abs() / magnitude;
        if ratio <= 1.0e-9 {
            continue;
        }
        outcome.max_leg_out_of_plane_ratio = outcome.max_leg_out_of_plane_ratio.max(ratio);
        let Some(direction) = (leg - plane_normal * out_of_plane).normalized() else {
            continue;
        };
        let (i, j) = second_row_grid_index(key.side, corner_k);
        if let Some(patch) = patches.get_mut(&key.face) {
            patch.grid[i][j] = row0[corner_k] + direction * magnitude;
            outcome.repaired_leg_count += 1;
        }
    }
}

/// Relative floor on the band's transverse (β·R) component: the
/// weak end's transverse magnitude is clamped to this fraction
/// of the strong end's, with matching sign. Without it, sliver
/// tips whose legs tilt almost tangent to the chain make β(t)
/// cross zero mid-row — ∂u ∥ ∂v there, a parametrically
/// degenerate stretch whose normals are noise on both sides of
/// the trim. The clamp keeps everything in-span (only the
/// in-plane shape near the corner changes), so G1 is unaffected
/// by construction.
const TRANSVERSE_FLOOR: f64 = 0.5;

/// Derive the full row-1 band from `(C', R)` and the exact end
/// legs: `D(t) = α(t)·C'(t) + β(t)·R(t)` with linear `α, β`
/// matched to `D(0) = E₀`, `D(1) = E₁`, expanded into degree-9
/// Bezier coefficients, then `row1ₖ = row0ₖ + Dₖ / 9`.
fn derive_row1(
    row0: &[Point3; 10],
    row1: &[Point3; 10],
    field_segment: &[Vec3; 9],
    clamped_count: &mut usize,
) -> Option<[Point3; 10]> {
    let derivative = row_derivative_coefficients(row0);
    let leg_start = (row1[0] - row0[0]) * 9.0;
    let leg_end = (row1[9] - row0[9]) * 9.0;
    let (alpha_start, mut beta_start) = solve_leg(leg_start, derivative[0], field_segment[0])?;
    let (alpha_end, mut beta_end) = solve_leg(leg_end, derivative[8], field_segment[8])?;
    // Transverse-floor clamp: each end's transverse magnitude
    // must reach a fraction of the LEG SCALE (not merely of the
    // other end — rows whose both legs are near-tangent would
    // otherwise stay collapsed), with one consistent sign taken
    // from the dominant end (or positive when both vanish).
    let leg_scale = leg_start
        .dot(&leg_start)
        .sqrt()
        .max(leg_end.dot(&leg_end).sqrt());
    let floor_magnitude = TRANSVERSE_FLOOR * leg_scale;
    let transverse_start = beta_start * field_segment[0].dot(&field_segment[0]).sqrt();
    let transverse_end = beta_end * field_segment[8].dot(&field_segment[8]).sqrt();
    let dominant_sign = if transverse_start.abs() >= transverse_end.abs() {
        if transverse_start >= 0.0 { 1.0 } else { -1.0 }
    } else if transverse_end >= 0.0 {
        1.0
    } else {
        -1.0
    };
    if transverse_start * dominant_sign < floor_magnitude {
        let magnitude = field_segment[0].dot(&field_segment[0]).sqrt().max(1.0e-300);
        beta_start = dominant_sign * floor_magnitude / magnitude;
        *clamped_count += 1;
    }
    if transverse_end * dominant_sign < floor_magnitude {
        let magnitude = field_segment[8].dot(&field_segment[8]).sqrt().max(1.0e-300);
        beta_end = dominant_sign * floor_magnitude / magnitude;
        *clamped_count += 1;
    }
    let mut out = *row1;
    for k in 1..9 {
        // Product of a linear blend with a degree-8 field:
        // coeffₖ = (9−k)/9·(α₀ d_k + β₀ r_k)
        //        +    k/9·(α₁ d_{k−1} + β₁ r_{k−1}).
        let head = derivative[k] * alpha_start + field_segment[k] * beta_start;
        let tail = derivative[k - 1] * alpha_end + field_segment[k - 1] * beta_end;
        let d_k = head * ((9 - k) as f64 / 9.0) + tail * (k as f64 / 9.0);
        out[k] = row0[k] + d_k * (1.0 / 9.0);
    }
    Some(out)
}

/// Solve `leg = α·u + β·v` in the least-squares sense (exact
/// when the leg is coplanar with the basis, which the node
/// tangent-plane constraints guarantee).
fn solve_leg(leg: Vec3, u: Vec3, v: Vec3) -> Option<(f64, f64)> {
    let uu = u.dot(&u);
    let uv = u.dot(&v);
    let vv = v.dot(&v);
    let det = uu * vv - uv * uv;
    if det.abs() <= LEG_SOLVE_DET_TOL * uu.max(vv).max(1.0e-300) {
        return None;
    }
    let bu = u.dot(&leg);
    let bv = v.dot(&leg);
    Some(((vv * bu - uv * bv) / det, (uu * bv - uv * bu) / det))
}

/// Discrete Coons fill of the interior controls (indices 2..=7)
/// from the four derived band rows/columns at indices 1 and 8.
fn refill_interior_coons(patch: &mut RectangularGregoryBiNonic) {
    let grid = patch.grid;
    let weighted = |terms: &[(f64, Point3)]| -> Point3 {
        let mut out = Point3::new(0.0, 0.0, 0.0);
        for &(weight, point) in terms {
            out = Point3::new(
                out.x + weight * point.x,
                out.y + weight * point.y,
                out.z + weight * point.z,
            );
        }
        out
    };
    for i in 2..8 {
        let a = (i as f64 - 1.0) / 7.0;
        for j in 2..8 {
            let b = (j as f64 - 1.0) / 7.0;
            patch.grid[i][j] = weighted(&[
                (1.0 - a, grid[1][j]),
                (a, grid[8][j]),
                (1.0 - b, grid[i][1]),
                (b, grid[i][8]),
                (-(1.0 - a) * (1.0 - b), grid[1][1]),
                (-(1.0 - a) * b, grid[1][8]),
                (-a * (1.0 - b), grid[8][1]),
                (-a * b, grid[8][8]),
            ]);
        }
    }
}

/// Relative out-of-span component of the derived band's
/// cross-derivative at interior samples: the residual of
/// `D(t) = 9·Σ (row1ᵢ − row0ᵢ) Bᵢ(t)` against the plane spanned
/// by the row tangent and the ribbon field.
/// The boundary side owning the transverse control at a row
/// corner: the perpendicular boundary through that corner.
pub(super) fn transverse_row_side(side: BoundarySide, corner_k: usize) -> BoundarySide {
    match (side, corner_k) {
        (BoundarySide::UMin | BoundarySide::UMax, 0) => BoundarySide::VMin,
        (BoundarySide::UMin | BoundarySide::UMax, _) => BoundarySide::VMax,
        (BoundarySide::VMin | BoundarySide::VMax, 0) => BoundarySide::UMin,
        (BoundarySide::VMin | BoundarySide::VMax, _) => BoundarySide::UMax,
    }
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

fn second_row_grid_index(side: BoundarySide, k: usize) -> (usize, usize) {
    match side {
        BoundarySide::UMin => (1, k),
        BoundarySide::UMax => (8, k),
        BoundarySide::VMin => (k, 1),
        BoundarySide::VMax => (k, 8),
    }
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_cell_reauthor_tests.rs"]
mod tests;
