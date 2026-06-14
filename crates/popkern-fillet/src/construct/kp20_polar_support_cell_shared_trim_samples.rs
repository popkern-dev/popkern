//! Ranked shared-trim sample diagnostics for support-cell quality.

use popkern_core::geometry::Vec3;

use super::Kp20SupportCellWorstTrim;

pub(in crate::construct) const SUPPORT_CELL_TOP_G1_SAMPLE_COUNT: usize = 4;
pub(in crate::construct) const SUPPORT_CELL_TOP_G2_SAMPLE_COUNT: usize = 4;

#[derive(Debug, Clone, Copy)]
pub(in crate::construct) struct Kp20SupportCellG1Sample {
    pub(in crate::construct) trim: Kp20SupportCellWorstTrim,
    pub(in crate::construct) normal_disagreement: f64,
    pub(in crate::construct) g0_gap: f64,
    pub(in crate::construct) uv_a: Vec3,
    pub(in crate::construct) uv_b: Vec3,
}

#[derive(Debug, Clone, Copy)]
pub(in crate::construct) struct Kp20SupportCellG2Sample {
    pub(in crate::construct) trim: Kp20SupportCellWorstTrim,
    pub(in crate::construct) g2_residual: f64,
    pub(in crate::construct) normal_aligned_g2_residual: f64,
    pub(in crate::construct) normal_alignment: f64,
    pub(in crate::construct) curvature_a: f64,
    pub(in crate::construct) curvature_b: f64,
    pub(in crate::construct) normal_disagreement: f64,
    pub(in crate::construct) g0_gap: f64,
    pub(in crate::construct) uv_a: Vec3,
    pub(in crate::construct) uv_b: Vec3,
}

pub(super) fn row_label(row: Option<u8>) -> String {
    row.map(|row| row.to_string())
        .unwrap_or_else(|| "none".to_string())
}

pub(super) fn g1_sample_detail(
    samples: &[Option<Kp20SupportCellG1Sample>; SUPPORT_CELL_TOP_G1_SAMPLE_COUNT],
) -> String {
    let entries = samples
        .iter()
        .flatten()
        .map(|sample| {
            format!(
                "{:?}:{}<->{:?}:{}@f={:.6e},g1={:.6e},g0={:.6e},uva=({:.4e},{:.4e}),uvb=({:.4e},{:.4e})",
                sample.trim.face_a,
                row_label(sample.trim.row_a),
                sample.trim.face_b,
                row_label(sample.trim.row_b),
                sample.trim.sample_fraction,
                sample.normal_disagreement,
                sample.g0_gap,
                sample.uv_a.x,
                sample.uv_a.y,
                sample.uv_b.x,
                sample.uv_b.y
            )
        })
        .collect::<Vec<_>>();
    if entries.is_empty() {
        "none".to_string()
    } else {
        entries.join("|")
    }
}

pub(super) fn g2_sample_detail(
    samples: &[Option<Kp20SupportCellG2Sample>; SUPPORT_CELL_TOP_G2_SAMPLE_COUNT],
) -> String {
    let entries = samples
        .iter()
        .flatten()
        .map(|sample| {
            format!(
                "{:?}:{}<->{:?}:{}@f={:.6e},g2={:.6e},g2_aligned={:.6e},normal_alignment={:.1},ka={:.6e},kb={:.6e},g1={:.6e},g0={:.6e},uva=({:.4e},{:.4e}),uvb=({:.4e},{:.4e})",
                sample.trim.face_a,
                row_label(sample.trim.row_a),
                sample.trim.face_b,
                row_label(sample.trim.row_b),
                sample.trim.sample_fraction,
                sample.g2_residual,
                sample.normal_aligned_g2_residual,
                sample.normal_alignment,
                sample.curvature_a,
                sample.curvature_b,
                sample.normal_disagreement,
                sample.g0_gap,
                sample.uv_a.x,
                sample.uv_a.y,
                sample.uv_b.x,
                sample.uv_b.y
            )
        })
        .collect::<Vec<_>>();
    if entries.is_empty() {
        "none".to_string()
    } else {
        entries.join("|")
    }
}
