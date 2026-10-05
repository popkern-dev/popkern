//! Construction-time G1 predictor for the chain constructive
//! cells: samples every glued trim's two cross-derivative fields
//! from the final control nets and reports the worst normal
//! disagreement, mirroring the gate's surface-level measurement
//! before the patches are committed. Also hosts the band span
//! self-check (`derived_span_residual`) and the cell-anisotropy
//! probe (`transverse_extent`) feeding the worst-sample anatomy.

use popkern_math::PortableMath;
use std::collections::HashMap;

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{FaceHandle, VertexHandle};

#[path = "kp20_polar_support_cell_cross_field_diagnostics.rs"]
mod cross_field_diagnostics;

use self::cross_field_diagnostics::{
    record_beta_sign_flip, record_cross_field_sample, sample_cross_field,
};
use super::BoundarySide;
use super::band_derive::derive_row1;
use super::band_field::author_ribbon_field;
use super::cell_reauthor::{
    AuthoredSubChain, CellReauthorOutcome, PredictedTopG1Sample, second_row_controls,
    second_row_grid_index,
};
use super::chain_author::boundary_controls;
use super::chain_discovery::ChainTrim;
use super::chain_fit::ChainRowKey;
use super::chain_ribbon_field::{
    eval_ribbon_field, ribbon_field_segment, row_derivative_coefficients,
};
use super::twist_resolution::dual_twist_corner_index;

/// Sample every glued trim's two cross-derivative fields from
/// the final control nets and report the worst normal
/// disagreement, mirroring the gate's measurement. The worst
/// sample's anatomy (per-side own-span residuals against the
/// live ribbon field, field values, intervals) is retained for
/// localization.
pub(super) fn predict_worst_trim_g1(
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub_chains: &[AuthoredSubChain],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
    node_normals: &HashMap<VertexHandle, Vec3>,
    outcome: &mut CellReauthorOutcome,
) {
    use super::row_fit::bernstein9;
    let chained_rows: std::collections::HashSet<ChainRowKey> = sub_chains
        .iter()
        .flat_map(|sub| sub.rows.iter().map(|&(key, _)| key))
        .collect();
    for sub in sub_chains {
        let Some(field) = author_ribbon_field(sub, node_normals, patches, &chained_rows) else {
            continue;
        };
        let row_intervals: HashMap<ChainRowKey, (f64, f64)> = sub.rows.iter().copied().collect();
        for trim in chain_trims {
            let (Some(&(start_a, end_a)), Some(&(start_b, end_b))) = (
                row_intervals.get(&trim.row_a),
                row_intervals.get(&trim.row_b),
            ) else {
                continue;
            };
            let (Some(patch_a), Some(patch_b)) =
                (patches.get(&trim.row_a.face), patches.get(&trim.row_b.face))
            else {
                continue;
            };
            let row0_a = boundary_controls(patch_a, trim.row_a.side);
            let mut row1_a = second_row_controls(patch_a, trim.row_a.side);
            let row0_b = boundary_controls(patch_b, trim.row_b.side);
            let mut row1_b = second_row_controls(patch_b, trim.row_b.side);
            realize_dual_twist_candidates(&mut row1_a, trim.row_a, outcome);
            realize_dual_twist_candidates(&mut row1_b, trim.row_b, outcome);
            let derivative_a = row_derivative_coefficients(&row0_a);
            let derivative_b = row_derivative_coefficients(&row0_b);
            let field_segment_a = ribbon_field_segment(&field, start_a, end_a);
            let field_segment_b = ribbon_field_segment(&field, start_b, end_b);
            let trim_span_a = (trim.interval_a.1 - trim.interval_a.0).abs();
            let trim_span_b = (trim.interval_b.1 - trim.interval_b.0).abs();
            if trim_span_a <= 1.0e-9 || trim_span_b <= 1.0e-9 {
                continue;
            }
            let mut previous_beta_sign_a = 0_i8;
            let mut previous_beta_sign_b = 0_i8;
            for sample in 0..9 {
                let fraction = sample as f64 / 8.0;
                let t_a = lerp(trim.interval_a, fraction).clamp(0.0, 1.0);
                let t_b = lerp(trim.interval_b, 1.0 - fraction).clamp(0.0, 1.0);
                let chain_param_a = row_to_chain((start_a, end_a), t_a);
                let chain_param_b = row_to_chain((start_b, end_b), t_b);
                let chain_param = 0.5 * (chain_param_a + chain_param_b);
                let basis_a = bernstein9(t_a);
                let basis_b = bernstein9(t_b);
                let tangent_a = eval_ribbon_field(&derivative_a, t_a);
                let tangent_b = eval_ribbon_field(&derivative_b, t_b);
                let mut cross_a = Vec3::new(0.0, 0.0, 0.0);
                let mut cross_b = Vec3::new(0.0, 0.0, 0.0);
                for i in 0..10 {
                    cross_a = cross_a + (row1_a[i] - row0_a[i]) * basis_a[i];
                    cross_b = cross_b + (row1_b[i] - row0_b[i]) * basis_b[i];
                }
                let field_a = eval_ribbon_field(&field_segment_a, t_a);
                let field_b = eval_ribbon_field(&field_segment_b, t_b);
                if let Some(sample_a) = sample_cross_field(cross_a * 9.0, tangent_a, field_a) {
                    record_cross_field_sample(outcome, sample_a, || {
                        format!(
                            "{:?}:{:?}@s={chain_param:.4},t={t_a:.4},interval=({start_a:.4},{end_a:.4}),tf={:.6e},beta={:.6e},beta_fraction={:.6e},|cross|={:.6e},|field|={:.6e}",
                            trim.row_a.face,
                            trim.row_a.side,
                            sample_a.tangential_fraction,
                            sample_a.beta,
                            sample_a.beta_magnitude_fraction,
                            sample_a.cross_magnitude,
                            sample_a.field_magnitude
                        )
                    });
                    record_beta_sign_flip(outcome, &mut previous_beta_sign_a, sample_a);
                }
                if let Some(sample_b) = sample_cross_field(cross_b * 9.0, tangent_b, field_b) {
                    record_cross_field_sample(outcome, sample_b, || {
                        format!(
                            "{:?}:{:?}@s={chain_param:.4},t={t_b:.4},interval=({start_b:.4},{end_b:.4}),tf={:.6e},beta={:.6e},beta_fraction={:.6e},|cross|={:.6e},|field|={:.6e}",
                            trim.row_b.face,
                            trim.row_b.side,
                            sample_b.tangential_fraction,
                            sample_b.beta,
                            sample_b.beta_magnitude_fraction,
                            sample_b.cross_magnitude,
                            sample_b.field_magnitude
                        )
                    });
                    record_beta_sign_flip(outcome, &mut previous_beta_sign_b, sample_b);
                }
                let (Some(normal_a), Some(normal_b)) = (
                    tangent_a.cross(&cross_a).normalized(),
                    tangent_b.cross(&cross_b).normalized(),
                ) else {
                    continue;
                };
                let disagreement = 1.0 - normal_a.dot(&normal_b).abs();
                let top_index = predicted_top_sample_index(outcome, disagreement);
                if disagreement > outcome.predicted_worst_trim_g1 || top_index.is_some() {
                    let field_value = eval_ribbon_field(&field, chain_param);
                    let span_residual = |cross: Vec3, tangent: Vec3| -> f64 {
                        let Some(plane) = tangent.cross(&field_value).normalized() else {
                            return -1.0;
                        };
                        let magnitude = cross.dot(&cross).sqrt();
                        if magnitude <= 1.0e-300 {
                            return -2.0;
                        }
                        cross.dot(&plane).abs() / magnitude
                    };
                    let polyline_length = |row: &[Point3; 10]| -> f64 {
                        (1..10).map(|i| row[i].distance_to(&row[i - 1])).sum()
                    };
                    // Tangential fraction |c·T̂|/|c| — near 1
                    // means the cross-derivative is almost
                    // parallel to the boundary and the normal
                    // direction is numerical noise.
                    let tangential_fraction = |cross: Vec3, tangent: Vec3| -> f64 {
                        let t_mag = tangent.dot(&tangent).sqrt();
                        let c_mag = cross.dot(&cross).sqrt();
                        if t_mag <= 1.0e-300 || c_mag <= 1.0e-300 {
                            return -1.0;
                        }
                        cross.dot(&tangent).abs() / (t_mag * c_mag)
                    };
                    // Corner angle: chain tangent vs the
                    // nearest frozen corner leg of each row
                    // (degrees). Near 0/180 = collapsed quad
                    // corner (parametric degeneracy).
                    let corner_leg_angle =
                        |row0: &[Point3; 10], row1: &[Point3; 10], t: f64| -> f64 {
                            let k = if t < 0.5 { 0 } else { 9 };
                            let leg = row1[k] - row0[k];
                            let chain_leg = if k == 0 {
                                row0[1] - row0[0]
                            } else {
                                row0[9] - row0[8]
                            };
                            let denom = (leg.dot(&leg) * chain_leg.dot(&chain_leg)).sqrt();
                            if denom <= 1.0e-300 {
                                return -1.0;
                            }
                            (leg.dot(&chain_leg) / denom).clamp(-1.0, 1.0).pacos() * 180.0
                                / std::f64::consts::PI
                        };
                    let transverse_a = transverse_extent(patch_a, trim.row_a.side);
                    let transverse_b = transverse_extent(patch_b, trim.row_b.side);
                    // Out-of-span decomposition of the cross
                    // by coefficient group: frozen end legs
                    // (k = 0, 9), the dual-twist-adjacent
                    // pair (k = 1, 8), and the derived
                    // interior (k = 2..=7). Each group's
                    // out-of-plane component (against the
                    // local span normal) and magnitude name
                    // WHICH writer's content leaves the span
                    // at the worst sample.
                    let decompose = |row0: &[Point3; 10],
                                     row1: &[Point3; 10],
                                     basis: &[f64; 10],
                                     tangent: Vec3|
                     -> [(f64, f64); 3] {
                        let span_normal = tangent.cross(&field_value).normalized();
                        let mut groups = [(Vec3::new(0.0, 0.0, 0.0)); 3];
                        for k in 0..10 {
                            let group = match k {
                                0 | 9 => 0,
                                1 | 8 => 1,
                                _ => 2,
                            };
                            groups[group] = groups[group] + (row1[k] - row0[k]) * (9.0 * basis[k]);
                        }
                        groups.map(|vector| {
                            let out_of_plane = span_normal
                                .map(|normal| vector.dot(&normal))
                                .unwrap_or(f64::NAN);
                            (out_of_plane, vector.dot(&vector).sqrt())
                        })
                    };
                    let decomp_a = decompose(&row0_a, &row1_a, &basis_a, tangent_a);
                    let decomp_b = decompose(&row0_b, &row1_b, &basis_b, tangent_b);
                    // Realized-vs-product deviation at the
                    // dual-twist-adjacent coefficients: the
                    // product is recomputed from this row's
                    // own (row0, frozen row1 ends, field
                    // segment); a nonzero k = 1/8 delta names
                    // the write path that replaced the
                    // in-span product value.
                    let product_deviation = |row0: &[Point3; 10],
                                             row1: &[Point3; 10],
                                             interval: (f64, f64)|
                     -> (f64, f64) {
                        let segment = ribbon_field_segment(&field, interval.0, interval.1);
                        let mut scratch_a = 0usize;
                        let mut scratch_b = 0usize;
                        match derive_row1(row0, row1, &segment, &mut scratch_a, &mut scratch_b) {
                            Some(product) => (
                                row1[1].distance_to(&product[1]),
                                row1[8].distance_to(&product[8]),
                            ),
                            None => (f64::NAN, f64::NAN),
                        }
                    };
                    let dev_a = product_deviation(&row0_a, &row1_a, (start_a, end_a));
                    let dev_b = product_deviation(&row0_b, &row1_b, (start_b, end_b));
                    // Frozen-leg out-of-CORNER-span residual:
                    // the product identity (cross in span
                    // pointwise) requires the frozen end legs
                    // to lie exactly in span(d, r) AT their
                    // corner; the residual against that
                    // corner's own plane is the identity's
                    // breakage budget.
                    let leg_residual = |row0: &[Point3; 10],
                                        row1: &[Point3; 10],
                                        interval: (f64, f64),
                                        t: f64|
                     -> f64 {
                        let segment = ribbon_field_segment(&field, interval.0, interval.1);
                        let derivative = row_derivative_coefficients(row0);
                        let (leg, d, r) = if t < 0.5 {
                            ((row1[0] - row0[0]) * 9.0, derivative[0], segment[0])
                        } else {
                            ((row1[9] - row0[9]) * 9.0, derivative[8], segment[8])
                        };
                        match d.cross(&r).normalized() {
                            Some(normal) => leg.dot(&normal).abs(),
                            None => f64::NAN,
                        }
                    };
                    let legres_a = leg_residual(&row0_a, &row1_a, (start_a, end_a), t_a);
                    let legres_b = leg_residual(&row0_b, &row1_b, (start_b, end_b), t_b);
                    let detail = format!(
                        "f={fraction:.4},s={chain_param:.4},ta={t_a:.4},tb={t_b:.4},ia=({start_a:.4},{end_a:.4}),ib=({start_b:.4},{end_b:.4}),oosA={:.3e},oosB={:.3e},|cA|={:.3e},|cB|={:.3e},fA={:?}:{:?},fB={:?}:{:?},alongA={:.3e},acrossA={:.3e},alongB={:.3e},acrossB={:.3e},tfA={:.3},tfB={:.3},angA={:.1},angB={:.1},dA=[frz={:.2e}/{:.2e} k18={:.2e}/{:.2e} mid={:.2e}/{:.2e}],dB=[frz={:.2e}/{:.2e} k18={:.2e}/{:.2e} mid={:.2e}/{:.2e}],devA=({:.2e};{:.2e}),devB=({:.2e};{:.2e}),legresA={:.2e},legresB={:.2e}",
                        span_residual(cross_a, tangent_a),
                        span_residual(cross_b, tangent_b),
                        cross_a.dot(&cross_a).sqrt(),
                        cross_b.dot(&cross_b).sqrt(),
                        trim.row_a.face,
                        trim.row_a.side,
                        trim.row_b.face,
                        trim.row_b.side,
                        polyline_length(&row0_a),
                        transverse_a,
                        polyline_length(&row0_b),
                        transverse_b,
                        tangential_fraction(cross_a, tangent_a),
                        tangential_fraction(cross_b, tangent_b),
                        corner_leg_angle(&row0_a, &row1_a, t_a),
                        corner_leg_angle(&row0_b, &row1_b, t_b),
                        decomp_a[0].0,
                        decomp_a[0].1,
                        decomp_a[1].0,
                        decomp_a[1].1,
                        decomp_a[2].0,
                        decomp_a[2].1,
                        decomp_b[0].0,
                        decomp_b[0].1,
                        decomp_b[1].0,
                        decomp_b[1].1,
                        decomp_b[2].0,
                        decomp_b[2].1,
                        dev_a.0,
                        dev_a.1,
                        dev_b.0,
                        dev_b.1,
                        legres_a,
                        legres_b
                    );
                    if let Some(index) = top_index {
                        insert_predicted_top_sample(outcome, index, disagreement, detail.clone());
                    }
                    if disagreement > outcome.predicted_worst_trim_g1 {
                        outcome.predicted_worst_trim_g1 = disagreement;
                        outcome.predicted_worst_detail = Some(detail);
                    }
                }
            }
        }
    }
}

fn lerp(interval: (f64, f64), fraction: f64) -> f64 {
    interval.0 + (interval.1 - interval.0) * fraction
}

fn row_to_chain(row_interval: (f64, f64), row_param: f64) -> f64 {
    row_interval.0 + (row_interval.1 - row_interval.0) * row_param
}

fn predicted_top_sample_index(outcome: &CellReauthorOutcome, disagreement: f64) -> Option<usize> {
    outcome.predicted_top_g1_samples.iter().position(|sample| {
        sample
            .as_ref()
            .is_none_or(|entry| disagreement > entry.disagreement)
    })
}

fn insert_predicted_top_sample(
    outcome: &mut CellReauthorOutcome,
    index: usize,
    disagreement: f64,
    detail: String,
) {
    for move_index in (index + 1..outcome.predicted_top_g1_samples.len()).rev() {
        outcome.predicted_top_g1_samples[move_index] =
            outcome.predicted_top_g1_samples[move_index - 1].take();
    }
    outcome.predicted_top_g1_samples[index] = Some(PredictedTopG1Sample {
        disagreement,
        detail,
    });
}

/// Replace a second row's corner-adjacent twist entries with the
/// boundary's OWN dual-twist candidates (when authored): the
/// rational blend reduces to that candidate on the boundary, so
/// this is exactly the cross-derivative the surface realizes —
/// the grid slot only holds the pair midpoint for the interior.
fn realize_dual_twist_candidates(
    row1: &mut [Point3; 10],
    key: ChainRowKey,
    outcome: &CellReauthorOutcome,
) {
    let Some(corners) = outcome.dual_twist_pairs.get(&key.face) else {
        return;
    };
    let u_constant = matches!(key.side, BoundarySide::UMin | BoundarySide::UMax);
    for k in [1usize, 8] {
        let (i, j) = second_row_grid_index(key.side, k);
        let Some(slot) = dual_twist_corner_index(i, j) else {
            continue;
        };
        if let Some((p_u, p_v)) = corners[slot] {
            row1[k] = if u_constant { p_u } else { p_v };
        }
    }
}

/// Arc-length of the boundary OPPOSITE-direction mid-line: the
/// patch's extent across the given side, measured along the
/// transverse mid-parameter polyline. Together with the side's
/// own arc-length this exposes cell anisotropy (sliver evidence).
fn transverse_extent(patch: &RectangularGregoryBiNonic, side: BoundarySide) -> f64 {
    let mut length = 0.0;
    for i in 1..10 {
        let (previous, current) = match side {
            BoundarySide::UMin | BoundarySide::UMax => (patch.grid[i - 1][4], patch.grid[i][4]),
            BoundarySide::VMin | BoundarySide::VMax => (patch.grid[4][i - 1], patch.grid[4][i]),
        };
        length += current.distance_to(&previous);
    }
    length
}

pub(super) fn derived_span_residual(
    row0: &[Point3; 10],
    derived: &[Point3; 10],
    field_segment: &[Vec3; 9],
) -> (f64, f64) {
    use super::row_fit::bernstein9;
    let derivative = row_derivative_coefficients(row0);
    let mut worst = (0.0_f64, 0.0_f64);
    for sample in [0.05_f64, 0.11, 0.2, 0.35, 0.5, 0.65, 0.8, 0.89, 0.95] {
        let tangent = eval_ribbon_field(&derivative, sample);
        let transverse = eval_ribbon_field(field_segment, sample);
        let Some(plane_normal) = tangent.cross(&transverse).normalized() else {
            continue;
        };
        let basis = bernstein9(sample);
        let mut delta = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..10 {
            delta = delta + (derived[i] - row0[i]) * basis[i];
        }
        let magnitude = delta.dot(&delta).sqrt();
        if magnitude <= 1.0e-300 {
            continue;
        }
        let residual = delta.dot(&plane_normal).abs() / magnitude;
        if residual > worst.0 {
            worst = (residual, sample);
        }
    }
    worst
}
