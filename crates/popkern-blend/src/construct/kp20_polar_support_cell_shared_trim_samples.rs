//! Ranked shared-trim sample diagnostics for support-cell quality.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::NormalCurvatureComponents;

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
    pub(in crate::construct) decomposition_a: Kp20SupportCellG2Decomposition,
    pub(in crate::construct) decomposition_b: Kp20SupportCellG2Decomposition,
    pub(in crate::construct) normal_disagreement: f64,
    pub(in crate::construct) g0_gap: f64,
    pub(in crate::construct) uv_a: Vec3,
    pub(in crate::construct) uv_b: Vec3,
}

#[derive(Debug, Clone, Copy)]
pub(in crate::construct) struct Kp20SupportCellG2Decomposition {
    pub(in crate::construct) projected_u: f64,
    pub(in crate::construct) projected_v: f64,
    pub(in crate::construct) normal: Vec3,
    pub(in crate::construct) directional_second: Vec3,
    pub(in crate::construct) second_form_l: f64,
    pub(in crate::construct) second_form_m: f64,
    pub(in crate::construct) second_form_n: f64,
}

impl Kp20SupportCellG2Decomposition {
    pub(in crate::construct) fn oriented(
        components: NormalCurvatureComponents,
        orientation: bool,
    ) -> Self {
        let sign = if orientation { 1.0 } else { -1.0 };
        Self {
            projected_u: components.projected_u,
            projected_v: components.projected_v,
            normal: components.normal * sign,
            directional_second: components.directional_second,
            second_form_l: components.second_form_l * sign,
            second_form_m: components.second_form_m * sign,
            second_form_n: components.second_form_n * sign,
        }
    }

    fn uu_term(self) -> f64 {
        self.second_form_l * self.projected_u * self.projected_u
    }

    fn uv_term(self) -> f64 {
        2.0 * self.second_form_m * self.projected_u * self.projected_v
    }

    fn vv_term(self) -> f64 {
        self.second_form_n * self.projected_v * self.projected_v
    }

    fn directional_second_magnitude(self) -> f64 {
        self.directional_second.dot(&self.directional_second).sqrt()
    }
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
            ) + &g2_decomposition_detail("a", sample.decomposition_a)
                + &g2_decomposition_detail("b", sample.decomposition_b)
                + &g2_row_domain_detail(sample)
                + &g2_mechanism_detail(
                    sample.decomposition_a,
                    sample.decomposition_b,
                    sample.normal_alignment,
                )
        })
        .collect::<Vec<_>>();
    if entries.is_empty() {
        "none".to_string()
    } else {
        entries.join("|")
    }
}

fn g2_row_domain_detail(sample: &Kp20SupportCellG2Sample) -> String {
    format!(
        ",row_domain=(a_side={},a_row_t={:.6e},a_solver_f={:.6e},b_side={},b_row_t={:.6e},b_solver_f={:.6e})",
        row_side_label(sample.trim.row_a),
        row_parameter(sample.trim.row_a, sample.uv_a),
        sample.trim.sample_fraction,
        row_side_label(sample.trim.row_b),
        row_parameter(sample.trim.row_b, sample.uv_b),
        1.0 - sample.trim.sample_fraction
    )
}

fn row_side_label(row: Option<u8>) -> &'static str {
    match row {
        Some(0) => "UMin",
        Some(1) => "UMax",
        Some(2) => "VMin",
        Some(3) => "VMax",
        _ => "none",
    }
}

fn row_parameter(row: Option<u8>, uv: Vec3) -> f64 {
    match row {
        Some(0 | 1) => uv.y,
        Some(2 | 3) => uv.x,
        _ => f64::NAN,
    }
}

fn g2_decomposition_detail(prefix: &str, decomposition: Kp20SupportCellG2Decomposition) -> String {
    format!(
        ",{prefix}_proj=(u={:.6e},v={:.6e}),{prefix}_q={:.6e},{prefix}_ii=(l={:.6e},m={:.6e},n={:.6e}),{prefix}_terms=(uu={:.6e},uv={:.6e},vv={:.6e})",
        decomposition.projected_u,
        decomposition.projected_v,
        decomposition.directional_second_magnitude(),
        decomposition.second_form_l,
        decomposition.second_form_m,
        decomposition.second_form_n,
        decomposition.uu_term(),
        decomposition.uv_term(),
        decomposition.vv_term()
    )
}

fn g2_mechanism_detail(
    first: Kp20SupportCellG2Decomposition,
    second: Kp20SupportCellG2Decomposition,
    normal_alignment: f64,
) -> String {
    format!(
        ",q_gap={:.6e},normal_only=(a={:.6e},b={:.6e})",
        directional_second_gap(first, second),
        normal_only_residual(first, second, normal_alignment, first.directional_second),
        normal_only_residual(first, second, normal_alignment, second.directional_second)
    )
}

fn directional_second_gap(
    first: Kp20SupportCellG2Decomposition,
    second: Kp20SupportCellG2Decomposition,
) -> f64 {
    let delta = first.directional_second - second.directional_second;
    let scale = first
        .directional_second_magnitude()
        .max(second.directional_second_magnitude());
    delta.dot(&delta).sqrt() / (1.0 + scale)
}

fn normal_only_residual(
    first: Kp20SupportCellG2Decomposition,
    second: Kp20SupportCellG2Decomposition,
    normal_alignment: f64,
    directional_second: Vec3,
) -> f64 {
    let first_projection = first.normal.dot(&directional_second);
    let second_projection = (second.normal * normal_alignment).dot(&directional_second);
    let scale = first_projection.abs().max(second_projection.abs());
    (first_projection - second_projection).abs() / (1.0 + scale)
}
