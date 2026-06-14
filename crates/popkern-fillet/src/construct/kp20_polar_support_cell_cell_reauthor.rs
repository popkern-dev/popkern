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

#[path = "kp20_polar_support_cell_beta_envelope.rs"]
mod beta_envelope;
#[path = "kp20_polar_support_cell_pair_normal_solve.rs"]
mod pair_normal_solve;

use self::beta_envelope::repair_beta_envelope;
use self::pair_normal_solve::repair_pair_normal_solve;
use super::BoundarySide;
use super::band_derive::derive_row1;
use super::band_field::author_ribbon_field;
use super::band_predictor::{derived_span_residual, predict_worst_trim_g1};
use super::chain_author::boundary_controls;
use super::chain_fit::ChainRowKey;
use super::chain_ribbon_field::{ribbon_field_segment, row_derivative_coefficients};
use super::twist_resolution::{
    TwistDemand, dual_twist_corner_index, resolve_dual_twist_corner, resolve_twist_demands,
};

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

#[derive(Debug, Clone)]
pub(super) struct DerivedRowBand {
    pub(super) key: ChainRowKey,
    pub(super) interval: (f64, f64),
    pub(super) row0: [Point3; 10],
    pub(super) row1: [Point3; 10],
    pub(super) field_segment: [Vec3; 9],
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
    /// Row identity carrying the worst derived span residual.
    pub(super) worst_span_residual_row: Option<String>,
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
    /// Cross-field samples measured from final control nets
    /// against each row's local `(C', R)` span. These diagnostics
    /// separate tangent-dominated row-1 bands (tiny beta content)
    /// from twist/out-of-span damage before a repair candidate is
    /// allowed to move controls.
    pub(super) cross_field_sample_count: usize,
    pub(super) near_tangent_cross_count: usize,
    pub(super) beta_sign_flip_count: usize,
    pub(super) max_tangential_fraction: f64,
    pub(super) min_beta_magnitude_fraction: f64,
    pub(super) worst_cross_field_detail: Option<String>,
    /// Optional beta-envelope repair candidate: counts sampled
    /// row interiors whose realized beta content was lifted by a
    /// constrained solve over controls k=2..7 only.
    pub(super) beta_envelope_sample_count: usize,
    pub(super) beta_envelope_after_sample_count: usize,
    pub(super) beta_envelope_constraint_count: usize,
    pub(super) beta_envelope_repaired_row_count: usize,
    pub(super) max_beta_envelope_control_delta: f64,
    pub(super) min_beta_envelope_fraction_before: f64,
    pub(super) min_beta_envelope_fraction_after: f64,
    pub(super) pair_normal_sample_count: usize,
    pub(super) pair_normal_constraint_count: usize,
    pub(super) pair_normal_repaired_row_count: usize,
    pub(super) pair_normal_limited_row_count: usize,
    pub(super) max_pair_normal_control_delta: f64,
    pub(super) max_pair_normal_disagreement_before: f64,
    pub(super) max_pair_normal_disagreement_after: f64,
    pub(super) worst_pair_normal_detail: Option<String>,
    /// End legs whose transverse component was floor-clamped to
    /// prevent a mid-row parametric degeneracy.
    pub(super) clamped_leg_beta_count: usize,
    /// Bands derived one-sided because a corner leg was chain-
    /// collinear (T-junction stems, sliver tips).
    pub(super) one_sided_band_count: usize,
    /// Dual-twist slots where one DIRECTION side received more
    /// than one demand — the mean then differs from each demand
    /// and the realized candidate leaves that row's span.
    pub(super) dual_multi_demand_side_count: usize,
    /// Dual-twist mode: per-face corner candidate pairs `(P_u,
    /// P_v)` indexed by the patch-type corner order `[(0,0),
    /// (0,9), (9,9), (9,0)]`. `None` corners keep the plain
    /// grid control.
    pub(super) dual_twist_pairs: HashMap<FaceHandle, [Option<(Point3, Point3)>; 4]>,
    /// Corner pairs whose two candidates genuinely differ —
    /// twist conflicts that the dual-twist patch absorbs
    /// structurally instead of compromising.
    pub(super) dual_twist_pair_count: usize,
    /// Worst separation between a pair's two candidates — the
    /// twist incompatibility magnitude the old single-CP
    /// resolution had to destroy.
    pub(super) max_dual_twist_pair_separation: f64,
}

/// Derive every chained row's row-1 band from the shared
/// `(C, R)` fields, then refill the touched cells' interiors.
pub(super) fn reauthor_cells_from_chain_network(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub_chains: &[AuthoredSubChain],
    node_normals: &HashMap<VertexHandle, Vec3>,
    dual_twist: bool,
    beta_envelope: bool,
    pair_normal_solve: bool,
) -> CellReauthorOutcome {
    let mut outcome = CellReauthorOutcome::default();
    let chained_rows: std::collections::HashSet<ChainRowKey> = sub_chains
        .iter()
        .flat_map(|sub| sub.rows.iter().map(|&(key, _)| key))
        .collect();
    // Collect every derived control write; twist controls are
    // demanded by two boundaries of the same face and averaged,
    // and the demand spread is retained for observability.
    let mut row_bands = Vec::new();
    for sub in sub_chains {
        let Some(field) = author_ribbon_field(sub, node_normals, patches, &chained_rows) else {
            outcome.field_declined_sub_chain_count += 1;
            continue;
        };
        let first_row = row_bands.len();
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
                &mut outcome.one_sided_band_count,
            ) {
                Some(mut derived) => {
                    if beta_envelope {
                        repair_beta_envelope(&row0, &mut derived, &field_segment, &mut outcome);
                    }
                    row_bands.push(DerivedRowBand {
                        key,
                        interval: (start, end),
                        row0,
                        row1: derived,
                        field_segment,
                    });
                    outcome.derived_row_count += 1;
                }
                None => outcome.leg_solve_declined_count += 1,
            }
        }
        if pair_normal_solve {
            repair_pair_normal_solve(&mut row_bands[first_row..], &mut outcome);
        }
    }
    let mut writes: HashMap<(FaceHandle, usize, usize), Vec<TwistDemand>> = HashMap::new();
    for band in &row_bands {
        let (span_residual, residual_param) =
            derived_span_residual(&band.row0, &band.row1, &band.field_segment);
        if span_residual > outcome.max_derived_span_residual {
            outcome.max_derived_span_residual = span_residual;
            // Worst-row anatomy: end-leg angles vs the
            // chain (sine), leg magnitudes vs the local
            // field, and the band magnitude at the
            // worst sample — names the deviation source
            // (degenerate direction, thin magnitude, or
            // healthy-but-clamped).
            let derivative = row_derivative_coefficients(&band.row0);
            let leg_start = (band.row1[0] - band.row0[0]) * 9.0;
            let leg_end = (band.row1[9] - band.row0[9]) * 9.0;
            let sine = |leg: Vec3, chain: Vec3| -> f64 {
                let denominator = (leg.dot(&leg) * chain.dot(&chain)).sqrt();
                if denominator <= 1.0e-300 {
                    return f64::NAN;
                }
                let cross = leg.cross(&chain);
                cross.dot(&cross).sqrt() / denominator
            };
            outcome.worst_span_residual_row = Some(format!(
                "{:?}:{:?}@{residual_param:.2}(sin0={:.2},sin1={:.2},|E0|={:.2e},|E1|={:.2e},|r0|={:.2e},|r8|={:.2e})",
                band.key.face,
                band.key.side,
                sine(leg_start, derivative[0]),
                sine(leg_end, derivative[8]),
                leg_start.dot(&leg_start).sqrt(),
                leg_end.dot(&leg_end).sqrt(),
                band.field_segment[0].dot(&band.field_segment[0]).sqrt(),
                band.field_segment[8].dot(&band.field_segment[8]).sqrt()
            ));
        }
        let derivative = row_derivative_coefficients(&band.row0);
        let corner_normal_start = derivative[0].cross(&band.field_segment[0]).normalized();
        let corner_normal_end = derivative[8].cross(&band.field_segment[8]).normalized();
        for (k, point) in band.row1.iter().enumerate().take(9).skip(1) {
            let (i, j) = second_row_grid_index(band.key.side, k);
            let plane_normal = match k {
                1 => corner_normal_start,
                8 => corner_normal_end,
                _ => None,
            };
            let band_offset = *point - band.row0[k];
            writes
                .entry((band.key.face, i, j))
                .or_default()
                .push(TwistDemand {
                    position: *point,
                    plane_normal,
                    field_magnitude: band_offset.dot(&band_offset).sqrt(),
                    u_constant_boundary: matches!(
                        band.key.side,
                        BoundarySide::UMin | BoundarySide::UMax
                    ),
                });
        }
    }
    let mut touched_faces = Vec::new();
    for (&(face, i, j), demands) in &writes {
        let resolved = if dual_twist && dual_twist_corner_index(i, j).is_some() {
            resolve_dual_twist_corner(face, i, j, demands, &mut outcome)
        } else {
            resolve_twist_demands(demands, &mut outcome)
        };
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

pub(super) fn second_row_grid_index(side: BoundarySide, k: usize) -> (usize, usize) {
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
