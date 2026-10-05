//! KP20 macro-cap validation metrics.

use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp20_bi6_assemble::Kp20Bi6MacroCapOutput;
use super::kp20_bi6_cap::{
    Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId, initialise_kp20_macro_patches_from_tensor_border,
    kp20_apply_g2_strip,
};
use super::kp20_bi6_math::kp20_b_coeffs_near_o;
use popkern_math::PortableMath;

/// Maximum residuals for the KP20 identities enforced by assembly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20Bi6ValidationReport {
    /// Largest position mismatch in equation (12)'s internal C2 rows.
    pub max_internal_c2_position_residual: f64,
    /// Largest first-cross-row mismatch in equation (12)'s internal C2 rows.
    pub max_internal_c2_first_residual: f64,
    /// Largest origin-strip row mismatch against equation (13).
    pub max_origin_g2_strip_residual: f64,
    /// Largest m-side tensor-border row mismatch against equation (14) split data.
    pub max_tensor_border_residual: f64,
    /// Largest absolute z coordinate across all bi-6 control points.
    pub max_abs_z: f64,
}

impl Kp20Bi6ValidationReport {
    /// Returns true when every measured residual is within `tolerance`.
    pub fn passes(self, tolerance: f64) -> bool {
        self.max_internal_c2_position_residual <= tolerance
            && self.max_internal_c2_first_residual <= tolerance
            && self.max_origin_g2_strip_residual <= tolerance
            && self.max_tensor_border_residual <= tolerance
    }
}

/// Error from KP20 validation metric construction.
#[derive(Debug, Clone, PartialEq)]
pub enum Kp20Bi6ValidationError {
    /// Output valence and tensor-border valence differ.
    ValenceMismatch {
        /// Number of assembled macro sectors.
        output: usize,
        /// Number of tensor-border sectors.
        tensor_border: usize,
    },
    /// Tensor-border split or validation failed.
    TensorBorder(String),
}

impl std::fmt::Display for Kp20Bi6ValidationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ValenceMismatch {
                output,
                tensor_border,
            } => write!(
                f,
                "KP20 validation valence mismatch: output has {output}, tensor-border has {tensor_border}"
            ),
            Self::TensorBorder(e) => write!(f, "KP20 validation tensor-border: {e}"),
        }
    }
}

impl std::error::Error for Kp20Bi6ValidationError {}

/// Measure the enforced continuity, boundary, and planarity residuals
/// of an assembled KP20 macro cap.
pub fn measure_kp20_bi6_macro_cap(
    output: &Kp20Bi6MacroCapOutput,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20Bi6ValidationReport, Kp20Bi6ValidationError> {
    if output.macro_patches.len() != tensor_border.sectors.len() {
        return Err(Kp20Bi6ValidationError::ValenceMismatch {
            output: output.macro_patches.len(),
            tensor_border: tensor_border.sectors.len(),
        });
    }

    let boundary = initialise_kp20_macro_patches_from_tensor_border(tensor_border)
        .map_err(Kp20Bi6ValidationError::TensorBorder)?;
    let mut report = Kp20Bi6ValidationReport {
        max_internal_c2_position_residual: 0.0,
        max_internal_c2_first_residual: 0.0,
        max_origin_g2_strip_residual: 0.0,
        max_tensor_border_residual: 0.0,
        max_abs_z: 0.0,
    };

    for (sector, expected_boundary) in output.macro_patches.iter().zip(boundary.iter()) {
        measure_sector_c2(sector, &mut report);
        measure_sector_g2(sector, output.n, &mut report);
        measure_sector_boundary(sector, expected_boundary, &mut report);
        measure_sector_planarity(sector, &mut report);
    }

    Ok(report)
}

fn measure_sector_c2(sector: &Kp20MacroPatch, report: &mut Kp20Bi6ValidationReport) {
    measure_c2_pair(
        sector.piece(Kp20PieceId::NearMidpointLeft),
        sector.piece(Kp20PieceId::NearOriginLeft),
        report,
    );
    measure_c2_pair(
        sector.piece(Kp20PieceId::NearMidpointRight),
        sector.piece(Kp20PieceId::NearOriginRight),
        report,
    );
}

fn measure_c2_pair(
    lower: &Kp20Bi6Patch,
    upper: &Kp20Bi6Patch,
    report: &mut Kp20Bi6ValidationReport,
) {
    // The C² propagation rule (paper eq 12) is overridden in the
    // apex region (i + j ≤ 3) by paper §4.4 step (i)'s Lemma 4
    // cubic-jet overlay. We measure C² only on the rows / columns
    // that the overlay does NOT touch.
    //
    // - upper.cps[0][j] is restored for j ≤ 3, propagated for j ≥ 4.
    // - upper.cps[1][j] is restored for j ≤ 2, propagated for j ≥ 3.
    // - upper.cps[5][j] is propagated for every j (the overlay
    //   never reaches row 5 since 5 + 0 > 3).
    for j in 0..7 {
        let half_delta = (lower.cps[4][j] - lower.cps[2][j]) * 0.25;
        if j >= 4 {
            push_max(
                &mut report.max_internal_c2_position_residual,
                upper.cps[0][j].distance_to(&lower.cps[6][j]),
            );
        }
        if j >= 3 {
            let expected1 = lower.cps[6][j] - half_delta;
            push_max(
                &mut report.max_internal_c2_first_residual,
                upper.cps[1][j].distance_to(&expected1),
            );
        }
        if j != 1 && j != 2 {
            let expected5 = lower.cps[6][j] + half_delta;
            push_max(
                &mut report.max_internal_c2_first_residual,
                upper.cps[5][j].distance_to(&expected5),
            );
        }
    }
}

fn measure_sector_g2(sector: &Kp20MacroPatch, n: usize, report: &mut Kp20Bi6ValidationReport) {
    let mut left = *sector.piece(Kp20PieceId::NearOriginLeft);
    let mut right = *sector.piece(Kp20PieceId::NearOriginRight);
    let c = (core::f64::consts::TAU / n as f64).pcos();
    let b0 = kp20_b_coeffs_near_o(c)[0];
    kp20_apply_g2_strip(&mut left, &mut right, b0);
    let actual_left = sector.piece(Kp20PieceId::NearOriginLeft);
    let actual_right = sector.piece(Kp20PieceId::NearOriginRight);
    for i in 4..=6 {
        push_max(
            &mut report.max_origin_g2_strip_residual,
            actual_left.cps[i][1].distance_to(&left.cps[i][1]),
        );
        push_max(
            &mut report.max_origin_g2_strip_residual,
            actual_right.cps[i][1].distance_to(&right.cps[i][1]),
        );
    }
}

fn measure_sector_boundary(
    sector: &Kp20MacroPatch,
    expected: &Kp20MacroPatch,
    report: &mut Kp20Bi6ValidationReport,
) {
    for id in [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ] {
        let actual_piece = sector.piece(id);
        let expected_piece = expected.piece(id);
        for row in 0..3 {
            for i in 0..7 {
                push_max(
                    &mut report.max_tensor_border_residual,
                    actual_piece.cps[row][i].distance_to(&expected_piece.cps[row][i]),
                );
            }
        }
    }
}

fn measure_sector_planarity(sector: &Kp20MacroPatch, report: &mut Kp20Bi6ValidationReport) {
    for id in Kp20PieceId::all() {
        for row in sector.piece(id).cps {
            for cp in row {
                push_max(&mut report.max_abs_z, cp.z.abs());
            }
        }
    }
}

fn push_max(slot: &mut f64, value: f64) {
    *slot = slot.max(value);
}

#[cfg(test)]
#[path = "kp20_bi6_validate_tests.rs"]
mod tests;
