//! Construction-time G1 predictor for the chain constructive
//! cells: samples every glued trim's two cross-derivative fields
//! from the final control nets and reports the worst normal
//! disagreement, mirroring the gate's surface-level measurement
//! before the patches are committed. Also hosts the band span
//! self-check (`derived_span_residual`) and the cell-anisotropy
//! probe (`transverse_extent`) feeding the worst-sample anatomy.

use std::collections::HashMap;

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::BoundarySide;
use super::band_field::author_ribbon_field;
use super::cell_reauthor::{AuthoredSubChain, CellReauthorOutcome, second_row_controls};
use super::chain_fit::{ChainRowKey, boundary_controls};
use super::chain_ribbon_field::{eval_ribbon_field, row_derivative_coefficients};

/// Sample every glued trim's two cross-derivative fields from
/// the final control nets and report the worst normal
/// disagreement, mirroring the gate's measurement. The worst
/// sample's anatomy (per-side own-span residuals against the
/// live ribbon field, field values, intervals) is retained for
/// localization.
pub(super) fn predict_worst_trim_g1(
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
    sub_chains: &[AuthoredSubChain],
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
        for first in 0..sub.rows.len() {
            for second in first + 1..sub.rows.len() {
                let (key_a, (start_a, end_a)) = sub.rows[first];
                let (key_b, (start_b, end_b)) = sub.rows[second];
                let (lo_a, hi_a) = (start_a.min(end_a), start_a.max(end_a));
                let (lo_b, hi_b) = (start_b.min(end_b), start_b.max(end_b));
                let lo = lo_a.max(lo_b);
                let hi = hi_a.min(hi_b);
                if hi - lo <= 1.0e-6 {
                    continue;
                }
                let (Some(patch_a), Some(patch_b)) =
                    (patches.get(&key_a.face), patches.get(&key_b.face))
                else {
                    continue;
                };
                let row0_a = boundary_controls(patch_a, key_a.side);
                let row1_a = second_row_controls(patch_a, key_a.side);
                let row0_b = boundary_controls(patch_b, key_b.side);
                let row1_b = second_row_controls(patch_b, key_b.side);
                for sample in 1..16 {
                    let chain_param = lo + (hi - lo) * (sample as f64 / 16.0);
                    let t_a = ((chain_param - start_a) / (end_a - start_a)).clamp(0.0, 1.0);
                    let t_b = ((chain_param - start_b) / (end_b - start_b)).clamp(0.0, 1.0);
                    let basis_a = bernstein9(t_a);
                    let basis_b = bernstein9(t_b);
                    let tangent = eval_ribbon_field(&row_derivative_coefficients(&row0_a), t_a);
                    let mut cross_a = Vec3::new(0.0, 0.0, 0.0);
                    let mut cross_b = Vec3::new(0.0, 0.0, 0.0);
                    for i in 0..10 {
                        cross_a = cross_a + (row1_a[i] - row0_a[i]) * basis_a[i];
                        cross_b = cross_b + (row1_b[i] - row0_b[i]) * basis_b[i];
                    }
                    let (Some(normal_a), Some(normal_b)) = (
                        tangent.cross(&cross_a).normalized(),
                        tangent.cross(&cross_b).normalized(),
                    ) else {
                        continue;
                    };
                    let disagreement = 1.0 - normal_a.dot(&normal_b).abs();
                    if disagreement > outcome.predicted_worst_trim_g1 {
                        outcome.predicted_worst_trim_g1 = disagreement;
                        let field_value = eval_ribbon_field(&field, chain_param);
                        let span_residual = |cross: Vec3| -> f64 {
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
                        let transverse_a = transverse_extent(patch_a, key_a.side);
                        let transverse_b = transverse_extent(patch_b, key_b.side);
                        outcome.predicted_worst_detail = Some(format!(
                            "s={chain_param:.4},ta={t_a:.4},tb={t_b:.4},ia=({start_a:.4},{end_a:.4}),ib=({start_b:.4},{end_b:.4}),oosA={:.3e},oosB={:.3e},|cA|={:.3e},|cB|={:.3e},fA={:?}:{:?},fB={:?}:{:?},alongA={:.3e},acrossA={:.3e},alongB={:.3e},acrossB={:.3e}",
                            span_residual(cross_a),
                            span_residual(cross_b),
                            cross_a.dot(&cross_a).sqrt(),
                            cross_b.dot(&cross_b).sqrt(),
                            key_a.face,
                            key_a.side,
                            key_b.face,
                            key_b.side,
                            polyline_length(&row0_a),
                            transverse_a,
                            polyline_length(&row0_b),
                            transverse_b
                        ));
                    }
                }
            }
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
) -> f64 {
    use super::row_fit::bernstein9;
    let derivative = row_derivative_coefficients(row0);
    let mut worst = 0.0_f64;
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
        worst = worst.max(delta.dot(&plane_normal).abs() / magnitude);
    }
    worst
}
