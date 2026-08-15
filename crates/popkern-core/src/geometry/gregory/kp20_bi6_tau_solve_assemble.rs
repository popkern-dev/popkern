//! Assembly of the KP20 planar τ constraint set.
//!
//! The low-level F3 solver accepts independent row families. This module
//! combines the paper row families for one sector while keeping the
//! Catmull-Clark characteristic tensor-border as an explicit input.

use super::bernstein_bi5::bernstein_bi5_all;
use super::bernstein_bi6::bernstein_bi6_all;
use super::direct_linalg::solve_dense_finite_guarded;
use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use super::kp20_bi6_piece_maps::{
    Kp20PlanarTauValidityReport, kp20_planar_bi6_tau_validity_report,
};
use super::kp20_bi6_tau_f3_solve::{
    Kp20PlanarTauAugmentedF3ScalarSolveInput, Kp20PlanarTauAugmentedScalarCoeff,
    Kp20PlanarTauAugmentedScalarLinearConstraint, Kp20PlanarTauCoeff, Kp20PlanarTauScalarCoeff,
    solve_kp20_planar_tau_augmented_f3_with_scalar_constraints,
};
use super::kp20_bi6_tau_f3_solve::{
    Kp20PlanarTauF3ScalarSolveInput, Kp20PlanarTauF3SolveError, Kp20PlanarTauF3SolveOutput,
    Kp20PlanarTauLinearConstraint, Kp20PlanarTauScalarLinearConstraint,
    kp20_planar_tau_fig8_indicated_averaging_constraints, kp20_planar_tau_internal_c2_constraints,
    kp20_planar_tau_sector_g1_g2_constraints, kp20_planar_tau_tensor_border_constraints,
    solve_kp20_planar_tau_f3_with_scalar_constraints,
    solve_kp20_planar_tau_f3_with_scalar_constraints_diagnosed,
};
use crate::geometry::Point3;

/// Secondary source-movement penalty for the §4.2 equation (14) projection.
///
/// This measures the F3-priority limit of the projected source row space. The
/// current limit is row-clean and cap-fold-free, but sampled determinant tests
/// show it collapses to a near-zero Jacobian. Do not route this as an accepted
/// nonsingular KP20 reparameterization without an additional paper-supported
/// source normalization/rank condition.
const SECTION42_SOURCE_PROJECTION_SECONDARY_WEIGHT: f64 = 1.0e-12;
const PAPER_TAU_VALIDITY_SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
const APPENDIX_A_FIG6C_SQUARE_DOMAIN_TOLERANCE: f64 = 1.0e-6;
const APPENDIX_A_FIG6C_GUIDE_DOMAIN_TOLERANCE: f64 = 1.0e-6;

/// Fully assembled paper row set for one planar KP20 τ sector.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauPaperRowSet {
    /// Vector rows applied independently to planar x/y coordinates.
    pub vector_constraints: Vec<Kp20PlanarTauLinearConstraint>,
    /// Scalar rows that may couple planar x/y coordinates.
    pub scalar_constraints: Vec<Kp20PlanarTauScalarLinearConstraint>,
}

/// KP20 paper row family used by the assembled planar τ solve.
#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20PlanarTauPaperRowFamily {
    /// Internal C² joins across the Fig. 8 four-piece macro-patch layout.
    ///
    /// Not emitted by the live §4.2 residual surface (its RC1-contaminated
    /// emitters were retired with the KP20 recovery plan); kept because
    /// frontier gates match on the full paper row-family vocabulary.
    #[allow(dead_code)]
    InternalC2,
    /// KP20 Fig. 8b indicated coefficient average.
    Fig8IndicatedAverage,
    /// KP20 equation (14) tensor-border rows.
    TensorBorderEq14,
    /// KP20 equations (8)/(9) along the near-origin sector edge.
    ///
    /// Vocabulary-only, like [`Kp20PlanarTauPaperRowFamily::InternalC2`].
    #[allow(dead_code)]
    SectorG1G2NearOrigin,
    /// KP20 equations (8)/(9) along the near-midpoint sector edge.
    ///
    /// Vocabulary-only, like [`Kp20PlanarTauPaperRowFamily::InternalC2`].
    #[allow(dead_code)]
    SectorG1G2NearMidpoint,
    /// KP20 §4.2 C⁰ rows: `\vec p_i0 = p_i0` and underlined analogue.
    Section42C0,
    /// KP20 equation (11), the degree-5 sector-separating curve condition.
    Section42Eq11Degree5,
    /// KP20 equation (12), the pre-solved internal C² strip join.
    Section42Eq12InternalC2,
    /// KP20 equation (13), the pre-solved lower strip compatibility rows.
    Section42Eq13Strip,
    /// Diagonal symmetry rows for the sector.
    DiagonalSymmetry,
}

/// Maximum residual for one KP20 paper row family.
#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20PlanarTauPaperFamilyResidual {
    /// Paper row family that was evaluated.
    pub family: Kp20PlanarTauPaperRowFamily,
    /// Maximum absolute scalar residual in this family.
    pub max_residual: f64,
    /// Constraint row that produced `max_residual` within this family.
    pub worst_row_index: usize,
    /// Coordinate axis for vector rows: `Some(0)` for x, `Some(1)` for y.
    /// Scalar row families use `None`.
    pub worst_axis: Option<usize>,
}

/// Diagnostic summary for the KP20 §4.2 equation (14) source projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20PlanarTauSection42SourceProjectionReport {
    /// Maximum residual in the augmented projection KKT solve.
    pub max_row_residual: f64,
    /// Largest Euclidean movement of one projected tensor-border target.
    pub max_source_delta: f64,
    /// Number of scalar tensor-border target coordinates adjusted.
    pub adjusted_scalar_targets: usize,
}

/// Input for the KP20 Appendix-A/Fig. 6c active-set τ projection.
#[derive(Debug, Clone)]
pub struct Kp20PlanarTauAppendixAFig6cSolveInput {
    /// Valence `n` of the extraordinary vertex.
    pub valence: usize,
    /// Sector index for the global planar spoke basis.
    pub sector: usize,
    /// Global planar τ seed. The paper-faithful seed is the Appendix-A refined
    /// central map sampled on the KP20 Fig. 8 paper quadrant layout.
    pub initial: Kp20MacroPatch,
    /// Small Tikhonov weight for row-space nullspace determinism.
    pub regularization_weight: f64,
}

/// Result of the KP20 Appendix-A/Fig. 6c active-set τ projection.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauAppendixAFig6cSolveOutput {
    /// Projected global planar τ macro-patch.
    pub tau_macro_patch: Kp20MacroPatch,
    /// Number of active-set linearization passes.
    pub iterations: usize,
    /// Remaining sampled square-domain constraints at termination.
    pub active_domain_constraints: usize,
    /// Remaining sampled orientation constraints at termination.
    pub active_orientation_constraints: usize,
    /// Largest coefficient movement in one active-set step.
    pub max_step_delta: f64,
    /// Maximum residual in the final active-set least-squares step.
    pub linear_residual: f64,
    /// Maximum residual of the enforced Fig. 6c §4.2 row set.
    pub max_row_residual: f64,
    /// Maximum movement of one projected source-row coefficient.
    pub max_source_delta: f64,
    /// Sampled square-domain boundedness in sector-local coordinates.
    pub square_domain: Kp20PlanarTauSquareDomainReport,
    /// Sampled orientation/nonsingularity in sector-local coordinates.
    pub orientation: super::kp20_bi6_piece_maps::Kp20PlanarTauOrientationReport,
}

/// Square-domain report for the KP20 Fig. 8 tensor-product subpatch layout.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauSquareDomainReport {
    /// Number of sampled image points.
    pub sample_count: usize,
    /// Number of samples outside `0 <= u <= 1`, `0 <= v <= 1`.
    pub outside_square_samples: usize,
    /// Minimum sampled local `u`.
    pub min_u: f64,
    /// Minimum sampled local `v`.
    pub min_v: f64,
    /// Maximum sampled local `u`.
    pub max_u: f64,
    /// Maximum sampled local `v`.
    pub max_v: f64,
    /// Largest positive square-bound violation.
    pub max_square_excess: f64,
    /// Stable label for the sample with the largest square violation.
    pub worst_sample: String,
}

impl Kp20PlanarTauSquareDomainReport {
    /// Whether every sampled point stays inside the Fig. 8 square domain.
    pub fn stays_inside_square(&self) -> bool {
        self.outside_square_samples == 0
    }
}

/// Input for assembling and solving one paper-style planar τ sector.
#[derive(Debug, Clone)]
pub struct Kp20PlanarTauPaperSolveInput<'a> {
    /// Valence `n` of the extraordinary vertex.
    pub valence: usize,
    /// Sector index whose equation (14) tensor-border rows should be used.
    pub sector: usize,
    /// Planar characteristic tensor-border. In the final paper path this is
    /// the Catmull-Clark characteristic-ring tensor-border.
    pub characteristic_tensor_border: &'a TensorBorderBi3,
    /// Initial τ candidate used only as the regularization target.
    pub initial: Kp20MacroPatch,
    /// Small Tikhonov weight for affine-nullspace determinism.
    pub regularization_weight: f64,
}

/// Error from assembling or solving the paper-style planar τ sector.
#[derive(Debug, Clone, PartialEq)]
pub enum Kp20PlanarTauPaperSolveError {
    /// Valence was invalid for the paper row formulas.
    InvalidValence(usize),
    /// The supplied tensor-border could not be split or indexed.
    TensorBorder(String),
    /// The underlying F3 KKT solve failed.
    Solve(Kp20PlanarTauF3SolveError),
    /// The solved τ satisfies the linear rows but is not a valid sampled
    /// reparameterization.
    DegenerateReparameterization(Box<Kp20PlanarTauValidityReport>),
}

impl std::fmt::Display for Kp20PlanarTauPaperSolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidValence(n) => write!(f, "invalid KP20 planar tau valence {n}"),
            Self::TensorBorder(message) => write!(f, "{message}"),
            Self::Solve(err) => write!(f, "{err}"),
            Self::DegenerateReparameterization(validity) => write!(
                f,
                "KP20 planar tau solve produced an invalid reparameterization: outside={}/{}, negative_det={}/{}, near_zero_det={}/{}, det=[{}, {}]",
                validity.domain.outside_triangle_samples,
                validity.domain.sample_count,
                validity.orientation.negative_determinants,
                validity.orientation.sample_count,
                validity.orientation.near_zero_determinants,
                validity.orientation.sample_count,
                validity.orientation.min_determinant,
                validity.orientation.max_determinant
            ),
        }
    }
}

impl std::error::Error for Kp20PlanarTauPaperSolveError {}

impl From<Kp20PlanarTauF3SolveError> for Kp20PlanarTauPaperSolveError {
    fn from(value: Kp20PlanarTauF3SolveError) -> Self {
        Self::Solve(value)
    }
}

/// Build the KP20 §4.3 planar τ row set for one sector.
///
/// This assembles the row families already verified independently:
/// internal C², equation (14) tensor-border rows, and equation (8)/(9) G¹/G²
/// rows using KP20 equation (10)'s two `b(u)` choices.
///
/// KP20 Fig. 8 imposes a small indicated coefficient averaging relation before
/// minimization. That relation is included here as a single vector row. It is
/// not the same as the stronger full-grid diagonal symmetry diagnostic rows
/// returned by [`kp20_planar_tau_diagonal_symmetry_constraints`], so those rows
/// are not part of this accepted paper row set.
pub fn kp20_planar_tau_paper_row_set(
    valence: usize,
    sector: usize,
    characteristic_tensor_border: &TensorBorderBi3,
) -> Result<Kp20PlanarTauPaperRowSet, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    let c = (core::f64::consts::TAU / valence as f64).cos();
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    vector_constraints.extend(
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?,
    );
    vector_constraints.extend(kp20_planar_tau_sector_g1_g2_constraints(
        Kp20PieceId::NearOriginLeft,
        Kp20PieceId::NearOriginRight,
        kp20_b_coeffs_near_o(c),
    ));
    vector_constraints.extend(kp20_planar_tau_sector_g1_g2_constraints(
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
        kp20_b_coeffs_near_m(c),
    ));
    Ok(Kp20PlanarTauPaperRowSet {
        vector_constraints,
        scalar_constraints: Vec::new(),
    })
}

/// Build the KP20 §4.2/§4.3 row set using the paper's pre-solved strip basis.
///
/// The raw equations (8)/(9) give the continuity equations before elimination.
/// KP20 §4.2 immediately replaces them with equations (11), (12), and (13).
/// This accepted rotational path therefore keeps the Fig. 8 average and
/// equation (14) tensor-border rows, then adds the pre-solved §4.2 scalar rows
/// across the previous sector's right boundary and current sector's left
/// boundary. The Fig. 8 symmetry used by the live row set is the paper's
/// indicated four-quadrant midpoint averaging; full-grid diagonal coefficient
/// equality is retained only as a diagnostic because it is not stated by KP20.
pub fn kp20_planar_tau_rotational_paper_row_set(
    valence: usize,
    sector: usize,
    characteristic_tensor_border: &TensorBorderBi3,
) -> Result<Kp20PlanarTauPaperRowSet, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let scalar_constraints = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let (projected_tensor_constraints, _projection_report) =
        kp20_project_rotational_paper_tensor_band(
            tensor_constraints,
            &vector_constraints,
            &scalar_constraints,
        )?;
    vector_constraints.extend(projected_tensor_constraints);
    Ok(Kp20PlanarTauPaperRowSet {
        vector_constraints,
        scalar_constraints,
    })
}

/// Diagnostic variant that asks whether the sector-local diagonal rows are
/// compatible with §4.2 if every equation (14) tensor target may move.
#[cfg(test)]
pub fn kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_tensor_targets(
    valence: usize,
    sector: usize,
    characteristic_tensor_border: &TensorBorderBi3,
) -> Result<Kp20PlanarTauPaperRowSet, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let free_coeffs = tensor_constraints
        .iter()
        .map(|row| row.terms[0].0)
        .collect::<Vec<_>>();
    kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs(
        valence,
        sector,
        characteristic_tensor_border,
        &free_coeffs,
    )
}

/// Diagnostic variant that projects selected equation (14) tensor targets while
/// enforcing sector-local diagonal symmetry.
#[cfg(test)]
pub fn kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs(
    valence: usize,
    sector: usize,
    characteristic_tensor_border: &TensorBorderBi3,
    free_coeffs: &[Kp20PlanarTauCoeff],
) -> Result<Kp20PlanarTauPaperRowSet, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let mut scalar_constraints = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    scalar_constraints.extend(kp20_planar_tau_sector_local_diagonal_symmetry_constraints(
        valence, sector,
    )?);
    let (projected_tensor_constraints, _projection_report) =
        kp20_project_section42_source_tensor_rows_with_free_coeffs(
            tensor_constraints,
            &vector_constraints,
            &scalar_constraints,
            free_coeffs,
        )?;
    vector_constraints.extend(projected_tensor_constraints);
    Ok(Kp20PlanarTauPaperRowSet {
        vector_constraints,
        scalar_constraints,
    })
}

/// Measure the source-target movement used by the rotational §4.2 projection.
pub fn kp20_planar_tau_rotational_section42_source_projection_report(
    valence: usize,
    sector: usize,
    characteristic_tensor_border: &TensorBorderBi3,
) -> Result<Kp20PlanarTauSection42SourceProjectionReport, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_planar_tau_fig8_indicated_averaging_constraints());
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let scalar_constraints = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let (_projected_tensor_constraints, projection_report) =
        kp20_project_rotational_paper_tensor_band(
            tensor_constraints,
            &vector_constraints,
            &scalar_constraints,
        )?;
    Ok(projection_report)
}

/// Diagonal symmetry rows expressed in the sector-local spoke basis.
///
/// KP20 Fig. 8's diagonal symmetry swaps the two sector parameters. The τ
/// coefficients are stored in global planar coordinates, so the symmetry rows
/// must first project each coefficient onto the sector's two spoke coordinates.
pub fn kp20_planar_tau_sector_local_diagonal_symmetry_constraints(
    valence: usize,
    sector: usize,
) -> Result<Vec<Kp20PlanarTauScalarLinearConstraint>, Kp20PlanarTauPaperSolveError> {
    let (local0, local1) = sector_local_spoke_components(valence, sector)?;
    let mut rows =
        kp20_planar_tau_sector_local_self_diagonal_symmetry_constraints(valence, sector)?;
    add_sector_local_pair_diagonal_symmetry(
        &mut rows,
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
        local0,
        local1,
    );
    Ok(rows)
}

/// Sector-local self-diagonal rows for the two Fig. 8 pieces that lie on the
/// sector diagonal.
///
/// This helper is kept separate so diagnostics can measure the self-diagonal
/// contribution independently from the off-diagonal Fig. 8 quadrant reflection.
pub fn kp20_planar_tau_sector_local_self_diagonal_symmetry_constraints(
    valence: usize,
    sector: usize,
) -> Result<Vec<Kp20PlanarTauScalarLinearConstraint>, Kp20PlanarTauPaperSolveError> {
    let (local0, local1) = sector_local_spoke_components(valence, sector)?;
    let mut rows = Vec::new();
    add_sector_local_self_diagonal_symmetry(&mut rows, Kp20PieceId::NearOriginLeft, local0, local1);
    add_sector_local_self_diagonal_symmetry(
        &mut rows,
        Kp20PieceId::NearMidpointRight,
        local0,
        local1,
    );
    Ok(rows)
}

fn sector_local_spoke_components(
    valence: usize,
    sector: usize,
) -> Result<([f64; 2], [f64; 2]), Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
    let local0 = [spoke0[1] / det, -spoke0[0] / det];
    let local1 = [-spoke1[1] / det, spoke1[0] / det];
    Ok((local0, local1))
}

fn add_sector_local_self_diagonal_symmetry(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    piece: Kp20PieceId,
    local0: [f64; 2],
    local1: [f64; 2],
) {
    for i in 0..=6 {
        rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
            sector_local_component_difference_terms(
                Kp20PlanarTauCoeff::new(piece, i, i),
                local0,
                Kp20PlanarTauCoeff::new(piece, i, i),
                local1,
            ),
        ));
        for j in (i + 1)..=6 {
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                sector_local_component_difference_terms(
                    Kp20PlanarTauCoeff::new(piece, i, j),
                    local0,
                    Kp20PlanarTauCoeff::new(piece, j, i),
                    local1,
                ),
            ));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                sector_local_component_difference_terms(
                    Kp20PlanarTauCoeff::new(piece, i, j),
                    local1,
                    Kp20PlanarTauCoeff::new(piece, j, i),
                    local0,
                ),
            ));
        }
    }
}

fn add_sector_local_pair_diagonal_symmetry(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    lower_right: Kp20PieceId,
    upper_left: Kp20PieceId,
    local0: [f64; 2],
    local1: [f64; 2],
) {
    for i in 0..=6 {
        for j in 0..=6 {
            let upper_i = 6 - j;
            let upper_j = 6 - i;
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                sector_local_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local0,
                    Kp20PlanarTauCoeff::new(upper_left, upper_i, upper_j),
                    local1,
                ),
            ));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(
                sector_local_component_difference_terms(
                    Kp20PlanarTauCoeff::new(lower_right, i, j),
                    local1,
                    Kp20PlanarTauCoeff::new(upper_left, upper_i, upper_j),
                    local0,
                ),
            ));
        }
    }
}

fn sector_local_component_difference_terms(
    lhs: Kp20PlanarTauCoeff,
    lhs_component: [f64; 2],
    rhs: Kp20PlanarTauCoeff,
    rhs_component: [f64; 2],
) -> Vec<(Kp20PlanarTauScalarCoeff, f64)> {
    vec![
        (Kp20PlanarTauScalarCoeff::new(lhs, 0), lhs_component[0]),
        (Kp20PlanarTauScalarCoeff::new(lhs, 1), lhs_component[1]),
        (Kp20PlanarTauScalarCoeff::new(rhs, 0), -rhs_component[0]),
        (Kp20PlanarTauScalarCoeff::new(rhs, 1), -rhs_component[1]),
    ]
}

fn kp20_project_rotational_paper_tensor_band(
    tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    fixed_vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
) -> Result<
    (
        Vec<Kp20PlanarTauLinearConstraint>,
        Kp20PlanarTauSection42SourceProjectionReport,
    ),
    Kp20PlanarTauPaperSolveError,
> {
    let free_coeffs = tensor_rows
        .iter()
        .map(|row| row.terms[0].0)
        .collect::<Vec<_>>();
    kp20_project_section42_source_tensor_rows_with_free_coeffs(
        tensor_rows,
        fixed_vector_constraints,
        scalar_constraints,
        &free_coeffs,
    )
}

fn kp20_project_section42_source_tensor_rows_with_free_coeffs(
    tensor_rows: Vec<Kp20PlanarTauLinearConstraint>,
    fixed_vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
    free_coeffs: &[Kp20PlanarTauCoeff],
) -> Result<
    (
        Vec<Kp20PlanarTauLinearConstraint>,
        Kp20PlanarTauSection42SourceProjectionReport,
    ),
    Kp20PlanarTauPaperSolveError,
> {
    let mut vector_constraints = fixed_vector_constraints.to_vec();
    let mut augmented_scalar_constraints = Vec::new();
    let mut free_scalar_targets = Vec::new();
    let mut extra_count = 0_usize;

    for (row_index, row) in tensor_rows.iter().enumerate() {
        let coeff = row.terms[0].0;
        if free_coeffs.contains(&coeff) {
            for axis in 0..=1 {
                augmented_scalar_constraints.push(Kp20PlanarTauAugmentedScalarLinearConstraint {
                    terms: vec![
                        (
                            Kp20PlanarTauAugmentedScalarCoeff::Tau(Kp20PlanarTauScalarCoeff::new(
                                coeff, axis,
                            )),
                            1.0,
                        ),
                        (Kp20PlanarTauAugmentedScalarCoeff::Extra(extra_count), -1.0),
                    ],
                    target: row.target[axis],
                });
                free_scalar_targets.push((row_index, axis));
                extra_count += 1;
            }
        } else {
            vector_constraints.push(row.clone());
        }
    }

    let solved = solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
        &Kp20PlanarTauAugmentedF3ScalarSolveInput {
            initial: Kp20MacroPatch::zero(),
            vector_constraints,
            scalar_constraints: scalar_constraints.to_vec(),
            augmented_scalar_constraints,
            extra_variable_count: extra_count,
            regularization_weight: 1.0e-10,
            extra_regularization_weight: SECTION42_SOURCE_PROJECTION_SECONDARY_WEIGHT,
        },
    );
    let solved = solved.map_err(Kp20PlanarTauPaperSolveError::Solve)?;

    let mut projected_tensor_rows = tensor_rows;
    for ((row_index, axis), delta) in free_scalar_targets
        .into_iter()
        .zip(solved.extra_values.iter().copied())
    {
        projected_tensor_rows[row_index].target[axis] += delta;
    }

    Ok((
        projected_tensor_rows,
        Kp20PlanarTauSection42SourceProjectionReport {
            max_row_residual: solved.max_constraint_residual,
            max_source_delta: max_point_delta_from_scalar_deltas(&solved.extra_values),
            adjusted_scalar_targets: solved.extra_values.len(),
        },
    ))
}

/// Solve one KP20 §4.3 paper-style planar τ sector against a supplied
/// characteristic tensor-border.
///
/// Skips the KKT condition diagnostic (`kkt_condition_number` is
/// `None`); use [`solve_kp20_planar_tau_paper_sector_f3_diagnosed`]
/// to measure it.
pub fn solve_kp20_planar_tau_paper_sector_f3(
    input: &Kp20PlanarTauPaperSolveInput<'_>,
) -> Result<Kp20PlanarTauF3SolveOutput, Kp20PlanarTauPaperSolveError> {
    let row_set = kp20_planar_tau_paper_row_set(
        input.valence,
        input.sector,
        input.characteristic_tensor_border,
    )?;
    Ok(solve_kp20_planar_tau_f3_with_scalar_constraints(
        &Kp20PlanarTauF3ScalarSolveInput {
            initial: input.initial,
            vector_constraints: row_set.vector_constraints,
            scalar_constraints: row_set.scalar_constraints,
            regularization_weight: input.regularization_weight,
        },
    )?)
}

/// As [`solve_kp20_planar_tau_paper_sector_f3`], additionally
/// measuring the exact-SVD condition number of the scalar KKT
/// matrix — an O(n³) diagnostic, for diagnostic runs only.
pub fn solve_kp20_planar_tau_paper_sector_f3_diagnosed(
    input: &Kp20PlanarTauPaperSolveInput<'_>,
) -> Result<Kp20PlanarTauF3SolveOutput, Kp20PlanarTauPaperSolveError> {
    let row_set = kp20_planar_tau_paper_row_set(
        input.valence,
        input.sector,
        input.characteristic_tensor_border,
    )?;
    Ok(solve_kp20_planar_tau_f3_with_scalar_constraints_diagnosed(
        &Kp20PlanarTauF3ScalarSolveInput {
            initial: input.initial,
            vector_constraints: row_set.vector_constraints,
            scalar_constraints: row_set.scalar_constraints,
            regularization_weight: input.regularization_weight,
        },
    )?)
}

/// Solve one KP20 §4.3 planar τ sector using adjacent-sector equation (8)/(9)
/// rows.
///
/// Skips the KKT condition diagnostic (`kkt_condition_number` is
/// `None`); use
/// [`solve_kp20_planar_tau_rotational_paper_sector_f3_diagnosed`] to
/// measure it.
pub fn solve_kp20_planar_tau_rotational_paper_sector_f3(
    input: &Kp20PlanarTauPaperSolveInput<'_>,
) -> Result<Kp20PlanarTauF3SolveOutput, Kp20PlanarTauPaperSolveError> {
    let row_set = kp20_planar_tau_rotational_paper_row_set(
        input.valence,
        input.sector,
        input.characteristic_tensor_border,
    )?;
    Ok(solve_kp20_planar_tau_f3_with_scalar_constraints(
        &Kp20PlanarTauF3ScalarSolveInput {
            initial: input.initial,
            vector_constraints: row_set.vector_constraints,
            scalar_constraints: row_set.scalar_constraints,
            regularization_weight: input.regularization_weight,
        },
    )?)
}

/// As [`solve_kp20_planar_tau_rotational_paper_sector_f3`],
/// additionally measuring the exact-SVD condition number of the
/// scalar KKT matrix — an O(n³) diagnostic, for diagnostic runs only.
pub fn solve_kp20_planar_tau_rotational_paper_sector_f3_diagnosed(
    input: &Kp20PlanarTauPaperSolveInput<'_>,
) -> Result<Kp20PlanarTauF3SolveOutput, Kp20PlanarTauPaperSolveError> {
    let row_set = kp20_planar_tau_rotational_paper_row_set(
        input.valence,
        input.sector,
        input.characteristic_tensor_border,
    )?;
    Ok(solve_kp20_planar_tau_f3_with_scalar_constraints_diagnosed(
        &Kp20PlanarTauF3ScalarSolveInput {
            initial: input.initial,
            vector_constraints: row_set.vector_constraints,
            scalar_constraints: row_set.scalar_constraints,
            regularization_weight: input.regularization_weight,
        },
    )?)
}

/// Solve one KP20 §4.3 planar τ sector and require the sampled
/// reparameterization gate to pass.
///
/// The unchecked solver remains available for source-candidate diagnostics. This
/// wrapper is the production-facing boundary: a candidate that satisfies the
/// linear paper rows by collapsing the Jacobian is declined with a structured
/// validity report instead of being treated as an accepted τ.
pub fn solve_kp20_planar_tau_rotational_paper_sector_f3_checked(
    input: &Kp20PlanarTauPaperSolveInput<'_>,
) -> Result<Kp20PlanarTauF3SolveOutput, Kp20PlanarTauPaperSolveError> {
    let solved = solve_kp20_planar_tau_rotational_paper_sector_f3(input)?;
    let local_tau =
        global_planar_patch_to_local(input.valence, input.sector, solved.tau_macro_patch)?;
    let validity = kp20_planar_bi6_tau_validity_report(
        &local_tau,
        &PAPER_TAU_VALIDITY_SAMPLES,
        APPENDIX_A_FIG6C_GUIDE_DOMAIN_TOLERANCE,
        1.0e-10,
    );
    if validity.is_strict_sector_bounded_reparameterization() {
        Ok(solved)
    } else {
        Err(Kp20PlanarTauPaperSolveError::DegenerateReparameterization(
            Box::new(validity),
        ))
    }
}

/// Project a KP20 Appendix-A refined τ seed onto the Fig. 6c §4.2 row space
/// while enforcing sampled square-domain and orientation validity.
///
/// This is the live construction for the paper-derived reparameterization
/// substrate. It deliberately does not include equation (14) tensor-border
/// rows: §4.4 uses those gray tensor-border coefficients as cap-geometry input
/// after τ is valid, not as source constraints that are allowed to collapse τ.
pub fn solve_kp20_planar_tau_appendix_a_fig6c_active_set(
    input: &Kp20PlanarTauAppendixAFig6cSolveInput,
) -> Result<Kp20PlanarTauAppendixAFig6cSolveOutput, Kp20PlanarTauPaperSolveError> {
    if input.valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(input.valence));
    }
    if input.sector != 0 {
        let local_initial =
            global_planar_patch_to_local(input.valence, input.sector, input.initial)?;
        let canonical_initial = local_planar_patch_to_global(input.valence, 0, local_initial)?;
        let canonical = solve_kp20_planar_tau_appendix_a_fig6c_active_set(
            &Kp20PlanarTauAppendixAFig6cSolveInput {
                valence: input.valence,
                sector: 0,
                initial: canonical_initial,
                regularization_weight: input.regularization_weight,
            },
        )?;
        let canonical_local =
            global_planar_patch_to_local(input.valence, 0, canonical.tau_macro_patch)?;
        let tau_macro_patch =
            local_planar_patch_to_global(input.valence, input.sector, canonical_local)?;
        return Ok(Kp20PlanarTauAppendixAFig6cSolveOutput {
            tau_macro_patch,
            ..canonical
        });
    }
    let map = kp20_appendix_a_refined_fig6c_strip_map(input.valence);
    let (row_seed, max_source_delta, _) = kp20_fig6c_boundary_row_basis_projected_source_patch(
        &input.initial,
        map,
        SECTION42_SOURCE_PROJECTION_SECONDARY_WEIGHT,
    );
    let rows = kp20_fig6c_full_paper_direction_rows(map);
    let report = kp20_fig6c_active_set_projection(
        &row_seed,
        &rows,
        input.valence,
        input.sector,
        &PAPER_TAU_VALIDITY_SAMPLES,
        input.regularization_weight,
    )?;
    let local_tau =
        global_planar_patch_to_local(input.valence, input.sector, report.tau_macro_patch)?;
    let square_domain = kp20_planar_tau_square_domain_report(
        &local_tau,
        &PAPER_TAU_VALIDITY_SAMPLES,
        APPENDIX_A_FIG6C_SQUARE_DOMAIN_TOLERANCE,
    );
    let validity = kp20_planar_bi6_tau_validity_report(
        &local_tau,
        &PAPER_TAU_VALIDITY_SAMPLES,
        APPENDIX_A_FIG6C_GUIDE_DOMAIN_TOLERANCE,
        1.0e-10,
    );
    if !validity.domain.stays_inside_triangle() || !validity.orientation.is_orientation_preserving()
    {
        return Err(Kp20PlanarTauPaperSolveError::DegenerateReparameterization(
            Box::new(validity),
        ));
    }
    Ok(Kp20PlanarTauAppendixAFig6cSolveOutput {
        tau_macro_patch: report.tau_macro_patch,
        iterations: report.iterations,
        active_domain_constraints: report.active_domain_constraints,
        active_orientation_constraints: report.active_orientation_constraints,
        max_step_delta: report.max_step_delta,
        linear_residual: report.linear_residual,
        max_row_residual: scalar_row_set_residual(&report.tau_macro_patch, &rows),
        max_source_delta,
        square_domain,
        orientation: validity.orientation,
    })
}

fn global_planar_patch_to_local(
    valence: usize,
    sector: usize,
    mut patch: Kp20MacroPatch,
) -> Result<Kp20MacroPatch, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    for piece_id in Kp20PieceId::all() {
        let piece = patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                let local = global_to_local_planar(valence, sector, piece.cps[i][j]);
                piece.cps[i][j] = Point3::new(local[0], local[1], 0.0);
            }
        }
    }
    Ok(patch)
}

fn local_planar_patch_to_global(
    valence: usize,
    sector: usize,
    mut patch: Kp20MacroPatch,
) -> Result<Kp20MacroPatch, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    for piece_id in Kp20PieceId::all() {
        let piece = patch.piece_mut(piece_id);
        for i in 0..=6 {
            for j in 0..=6 {
                piece.cps[i][j] =
                    local_to_global_planar(valence, sector, [piece.cps[i][j].x, piece.cps[i][j].y]);
            }
        }
    }
    Ok(patch)
}

fn local_to_global_planar(valence: usize, sector: usize, local: [f64; 2]) -> Point3 {
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    Point3::new(
        local[0] * spoke1[0] + local[1] * spoke0[0],
        local[0] * spoke1[1] + local[1] * spoke0[1],
        0.0,
    )
}

fn global_to_local_planar(valence: usize, sector: usize, point: Point3) -> [f64; 2] {
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = [(sector as f64 * alpha).cos(), (sector as f64 * alpha).sin()];
    let spoke1 = [
        ((sector + 1) as f64 * alpha).cos(),
        ((sector + 1) as f64 * alpha).sin(),
    ];
    let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
    [
        (point.x * spoke0[1] - point.y * spoke0[0]) / det,
        (spoke1[0] * point.y - spoke1[1] * point.x) / det,
    ]
}

/// Evaluate the rotational KP20 §4.2 row set against each accepted paper family.
#[cfg(test)]
pub fn kp20_planar_tau_rotational_paper_family_residuals(
    valence: usize,
    sector: usize,
    characteristic_tensor_border: &TensorBorderBi3,
    tau_macro_patch: &Kp20MacroPatch,
) -> Result<Vec<Kp20PlanarTauPaperFamilyResidual>, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let map = kp20_section42_oriented_endpoint_aware_map();
    let mut c0_rows = Vec::new();
    add_kp20_section42_oriented_c0_rows(&mut c0_rows, map);
    let mut eq11_rows = Vec::new();
    add_kp20_section42_oriented_eq11_degree5_rows(&mut eq11_rows, map);
    let mut eq12_rows = Vec::new();
    add_kp20_section42_oriented_endpoint_aware_eq12_rows(&mut eq12_rows, map);
    let mut eq13_rows = Vec::new();
    add_kp20_section42_oriented_eq13_rows(&mut eq13_rows, map, 2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let scalar_constraints = kp20_section42_oriented_endpoint_aware_strip_basis_rows(2.0 * c);
    let diagonal_rows =
        kp20_planar_tau_sector_local_diagonal_symmetry_constraints(valence, sector)?;
    let tensor_border_rows =
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let (tensor_border, _projection_report) = kp20_project_rotational_paper_tensor_band(
        tensor_border_rows,
        &fig8_rows,
        &scalar_constraints,
    )?;
    Ok(vec![
        vector_family_residual(
            Kp20PlanarTauPaperRowFamily::Fig8IndicatedAverage,
            &fig8_rows,
            tau_macro_patch,
        ),
        vector_family_residual(
            Kp20PlanarTauPaperRowFamily::TensorBorderEq14,
            &tensor_border,
            tau_macro_patch,
        ),
        scalar_family_residual(
            Kp20PlanarTauPaperRowFamily::Section42C0,
            &c0_rows,
            tau_macro_patch,
        ),
        scalar_family_residual(
            Kp20PlanarTauPaperRowFamily::Section42Eq11Degree5,
            &eq11_rows,
            tau_macro_patch,
        ),
        scalar_family_residual(
            Kp20PlanarTauPaperRowFamily::Section42Eq12InternalC2,
            &eq12_rows,
            tau_macro_patch,
        ),
        scalar_family_residual(
            Kp20PlanarTauPaperRowFamily::Section42Eq13Strip,
            &eq13_rows,
            tau_macro_patch,
        ),
        scalar_family_residual(
            Kp20PlanarTauPaperRowFamily::DiagonalSymmetry,
            &diagonal_rows,
            tau_macro_patch,
        ),
    ])
}

#[derive(Debug, Clone, Copy)]
#[cfg(test)]
struct Kp20Section42PaperSymbolMap {
    p: Kp20Section42MappedPiece,
    arrow_p: Kp20Section42MappedPiece,
    underlined_p: Kp20Section42MappedPiece,
    arrow_underlined_p: Kp20Section42MappedPiece,
    previous_sector_rotation_radians: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Section42MappedPiece {
    piece: Kp20PieceId,
    is_previous_sector: bool,
}

impl Kp20Section42MappedPiece {
    const fn previous(piece: Kp20PieceId) -> Self {
        Self {
            piece,
            is_previous_sector: true,
        }
    }

    const fn current(piece: Kp20PieceId) -> Self {
        Self {
            piece,
            is_previous_sector: false,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Kp20Section42PaperSymbol {
    P,
    ArrowP,
    UnderlinedP,
    ArrowUnderlinedP,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Section42PaperTerm {
    symbol: Kp20Section42PaperSymbol,
    i: usize,
    j: usize,
    weight: f64,
}

#[derive(Debug, Clone, Copy)]
enum Kp20Section42StripEdge {
    I0,
    J0,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Section42Fig6cStripMap {
    mapped: Kp20Section42MappedPiece,
    edge: Kp20Section42StripEdge,
    sector_rotation_steps: i32,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Section42Fig6cSymbolMap {
    p: Kp20Section42Fig6cStripMap,
    arrow_p: Kp20Section42Fig6cStripMap,
    underlined_p: Kp20Section42Fig6cStripMap,
    arrow_underlined_p: Kp20Section42Fig6cStripMap,
    previous_sector_rotation_radians: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Fig6cSourceRowVariable {
    strip: Kp20Section42Fig6cStripMap,
    paper_i: usize,
    paper_j: usize,
    axis: usize,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Fig6cDomainActiveConstraint {
    piece_id: Kp20PieceId,
    u: f64,
    v: f64,
    local_axis_weights: [f64; 2],
    target: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Fig6cOrientationActiveConstraint {
    piece_id: Kp20PieceId,
    u: f64,
    v: f64,
    target_det: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Fig6cActiveSetProjectionReport {
    tau_macro_patch: Kp20MacroPatch,
    iterations: usize,
    active_domain_constraints: usize,
    active_orientation_constraints: usize,
    max_step_delta: f64,
    linear_residual: f64,
}

#[cfg(test)]
fn kp20_section42_standard_map(angle: f64) -> Kp20Section42PaperSymbolMap {
    Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    }
}

fn kp20_appendix_a_refined_fig6c_strip_map(valence: usize) -> Kp20Section42Fig6cSymbolMap {
    let angle = core::f64::consts::TAU / valence as f64;
    Kp20Section42Fig6cSymbolMap {
        p: Kp20Section42Fig6cStripMap {
            mapped: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginLeft),
            edge: Kp20Section42StripEdge::I0,
            sector_rotation_steps: 1,
        },
        arrow_p: Kp20Section42Fig6cStripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
            edge: Kp20Section42StripEdge::J0,
            sector_rotation_steps: 0,
        },
        underlined_p: Kp20Section42Fig6cStripMap {
            mapped: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginRight),
            edge: Kp20Section42StripEdge::J0,
            sector_rotation_steps: 0,
        },
        arrow_underlined_p: Kp20Section42Fig6cStripMap {
            mapped: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
            edge: Kp20Section42StripEdge::I0,
            sector_rotation_steps: 1,
        },
        previous_sector_rotation_radians: angle,
    }
}

#[derive(Debug, Clone, Copy)]
struct Kp20Section42OrientedStripMap {
    piece: Kp20PieceId,
    reverse_i: bool,
    reverse_j: bool,
    transpose: bool,
}

#[derive(Debug, Clone, Copy)]
struct Kp20Section42OrientedPaperSymbolMap {
    p: Kp20Section42OrientedStripMap,
    arrow_p: Kp20Section42OrientedStripMap,
    underlined_p: Kp20Section42OrientedStripMap,
    arrow_underlined_p: Kp20Section42OrientedStripMap,
}

fn kp20_section42_oriented_endpoint_aware_map() -> Kp20Section42OrientedPaperSymbolMap {
    Kp20Section42OrientedPaperSymbolMap {
        p: Kp20Section42OrientedStripMap {
            piece: Kp20PieceId::NearOriginLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_p: Kp20Section42OrientedStripMap {
            piece: Kp20PieceId::NearOriginRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
        underlined_p: Kp20Section42OrientedStripMap {
            piece: Kp20PieceId::NearMidpointLeft,
            reverse_i: false,
            reverse_j: true,
            transpose: true,
        },
        arrow_underlined_p: Kp20Section42OrientedStripMap {
            piece: Kp20PieceId::NearMidpointRight,
            reverse_i: false,
            reverse_j: false,
            transpose: true,
        },
    }
}

fn kp20_section42_oriented_endpoint_aware_strip_basis_rows(
    b0_near_origin: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let map = kp20_section42_oriented_endpoint_aware_map();
    let mut rows = Vec::new();
    add_kp20_section42_oriented_c0_rows(&mut rows, map);
    add_kp20_section42_oriented_eq11_degree5_rows(&mut rows, map);
    add_kp20_section42_oriented_endpoint_aware_eq12_rows(&mut rows, map);
    add_kp20_section42_oriented_eq13_rows(&mut rows, map, b0_near_origin);
    rows
}

#[cfg(test)]
fn kp20_section42_presolved_strip_basis_rows(
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_c0_rows(&mut rows, map);
    add_kp20_section42_eq11_degree5_rows(&mut rows, map);
    add_kp20_section42_eq12_internal_c2_rows(&mut rows, map);
    add_kp20_section42_eq13_rows(&mut rows, map, b0_near_origin);
    rows
}

fn kp20_fig6c_full_paper_direction_rows(
    map: Kp20Section42Fig6cSymbolMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let b0_near_origin = 2.0 * map.previous_sector_rotation_radians.cos();
    let mut rows = Vec::new();
    add_kp20_fig6c_c0_rows(&mut rows, map);
    add_kp20_fig6c_eq11_degree5_rows(&mut rows, map);
    add_kp20_fig6c_paper_direction_eq12_rows(&mut rows, map);
    add_kp20_fig6c_eq13_rows(&mut rows, map, b0_near_origin);
    rows
}

fn add_kp20_fig6c_c0_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42Fig6cSymbolMap,
) {
    for i in 0..=6 {
        add_kp20_fig6c_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
            ],
        );
        add_kp20_fig6c_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
            ],
        );
    }
}

fn add_kp20_fig6c_eq11_degree5_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42Fig6cSymbolMap,
) {
    const DEGREE6_ALT_BINOMIAL: [f64; 7] = [1.0, -6.0, 15.0, -20.0, 15.0, -6.0, 1.0];
    for symbol in [
        Kp20Section42PaperSymbol::P,
        Kp20Section42PaperSymbol::UnderlinedP,
    ] {
        let terms = DEGREE6_ALT_BINOMIAL
            .iter()
            .enumerate()
            .map(|(i, &weight)| section42_term(symbol, i, 0, weight))
            .collect::<Vec<_>>();
        add_kp20_fig6c_vector_rows(rows, map, &terms);
    }
}

fn add_kp20_fig6c_paper_direction_eq12_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42Fig6cSymbolMap,
) {
    for (source, receiving) in [
        (
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ] {
        for j in 0..=2 {
            add_kp20_fig6c_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 0, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            add_kp20_fig6c_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 5, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, -0.25),
                    section42_term(source, 2, j, 0.25),
                ],
            );
            add_kp20_fig6c_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 1, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, 0.25),
                    section42_term(source, 2, j, -0.25),
                ],
            );
        }
    }
}

fn add_kp20_fig6c_eq13_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42Fig6cSymbolMap,
    b0_near_origin: f64,
) {
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = (0.5 * b0_near_origin) / 30.0;

    for i in [5, 6] {
        add_kp20_fig6c_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, -0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, 0.25),
            ],
        );
        add_kp20_fig6c_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, 0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, -0.25),
            ],
        );
    }

    add_kp20_fig6c_vector_rows(
        rows,
        map,
        &[
            section42_term(Kp20Section42PaperSymbol::P, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, -w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, w1),
        ],
    );
    add_kp20_fig6c_vector_rows(
        rows,
        map,
        &[
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, -w1),
        ],
    );
}

fn max_point_delta_from_scalar_deltas(values: &[f64]) -> f64 {
    values
        .chunks(2)
        .map(|chunk| {
            let x = chunk[0];
            let y = chunk.get(1).copied().unwrap_or(0.0);
            (x * x + y * y).sqrt()
        })
        .fold(0.0_f64, f64::max)
}

#[cfg(test)]
fn kp20_section42_only_c0_rows(
    map: Kp20Section42PaperSymbolMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_c0_rows(&mut rows, map);
    rows
}

#[cfg(test)]
fn kp20_section42_only_eq11_rows(
    map: Kp20Section42PaperSymbolMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_eq11_degree5_rows(&mut rows, map);
    rows
}

#[cfg(test)]
fn kp20_section42_only_eq12_rows(
    map: Kp20Section42PaperSymbolMap,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_eq12_internal_c2_rows(&mut rows, map);
    rows
}

#[cfg(test)]
fn kp20_section42_only_eq13_rows(
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_kp20_section42_eq13_rows(&mut rows, map, b0_near_origin);
    rows
}

#[cfg(test)]
fn add_kp20_section42_c0_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
) {
    for i in 0..=6 {
        add_kp20_section42_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
            ],
        );
        add_kp20_section42_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
            ],
        );
    }
}

fn add_kp20_section42_oriented_c0_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42OrientedPaperSymbolMap,
) {
    for i in 0..=6 {
        add_kp20_section42_oriented_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::P, i, 0, -1.0),
            ],
        );
        add_kp20_section42_oriented_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 0, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
            ],
        );
    }
}

#[cfg(test)]
fn add_kp20_section42_eq11_degree5_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
) {
    const DEGREE6_ALT_BINOMIAL: [f64; 7] = [1.0, -6.0, 15.0, -20.0, 15.0, -6.0, 1.0];
    for symbol in [
        Kp20Section42PaperSymbol::P,
        Kp20Section42PaperSymbol::UnderlinedP,
    ] {
        let terms = DEGREE6_ALT_BINOMIAL
            .iter()
            .enumerate()
            .map(|(i, &weight)| section42_term(symbol, i, 0, weight))
            .collect::<Vec<_>>();
        add_kp20_section42_vector_rows(rows, map, &terms);
    }
}

fn add_kp20_section42_oriented_eq11_degree5_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42OrientedPaperSymbolMap,
) {
    const DEGREE6_ALT_BINOMIAL: [f64; 7] = [1.0, -6.0, 15.0, -20.0, 15.0, -6.0, 1.0];
    for symbol in [
        Kp20Section42PaperSymbol::P,
        Kp20Section42PaperSymbol::UnderlinedP,
    ] {
        let terms = DEGREE6_ALT_BINOMIAL
            .iter()
            .enumerate()
            .map(|(i, &weight)| section42_term(symbol, i, 0, weight))
            .collect::<Vec<_>>();
        add_kp20_section42_oriented_vector_rows(rows, map, &terms);
    }
}

#[cfg(test)]
fn add_kp20_section42_eq12_internal_c2_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
) {
    // KP20 eq. (12) writes rows 0, 1, and 5 of the upper/origin piece
    // from rows 2, 4, and 6 of the lower piece closer to m.
    for (upper, lower) in [
        (
            Kp20Section42PaperSymbol::UnderlinedP,
            Kp20Section42PaperSymbol::P,
        ),
        (
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
            Kp20Section42PaperSymbol::ArrowP,
        ),
    ] {
        for j in 0..=2 {
            add_kp20_section42_vector_rows(
                rows,
                map,
                &[
                    section42_term(lower, 0, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            add_kp20_section42_vector_rows(
                rows,
                map,
                &[
                    section42_term(lower, 5, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, -0.25),
                    section42_term(upper, 2, j, 0.25),
                ],
            );
            add_kp20_section42_vector_rows(
                rows,
                map,
                &[
                    section42_term(lower, 1, j, 1.0),
                    section42_term(upper, 6, j, -1.0),
                    section42_term(upper, 4, j, 0.25),
                    section42_term(upper, 2, j, -0.25),
                ],
            );
        }
    }
}

fn add_kp20_section42_oriented_endpoint_aware_eq12_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42OrientedPaperSymbolMap,
) {
    for (source, receiving) in [
        (
            Kp20Section42PaperSymbol::P,
            Kp20Section42PaperSymbol::UnderlinedP,
        ),
        (
            Kp20Section42PaperSymbol::ArrowP,
            Kp20Section42PaperSymbol::ArrowUnderlinedP,
        ),
    ] {
        for j in 0..=2 {
            add_kp20_section42_oriented_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 0, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                ],
            );
        }
        for j in 1..=2 {
            add_kp20_section42_oriented_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 1, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, -0.25),
                    section42_term(source, 2, j, 0.25),
                ],
            );
            add_kp20_section42_oriented_vector_rows(
                rows,
                map,
                &[
                    section42_term(receiving, 5, j, 1.0),
                    section42_term(source, 6, j, -1.0),
                    section42_term(source, 4, j, 0.25),
                    section42_term(source, 2, j, -0.25),
                ],
            );
        }
    }
}

#[cfg(test)]
fn add_kp20_section42_eq13_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
) {
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    // KP20 eq. (13): w3 uses the underlined near-m coefficient b0 = c,
    // while w1/w2 use the near-o coefficient b0 = 2c.
    let b0_near_midpoint = 0.5 * b0_near_origin;
    let w3 = b0_near_midpoint / 30.0;

    for i in [5, 6] {
        add_kp20_section42_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, -0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, 0.25),
            ],
        );
        add_kp20_section42_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, 0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, -0.25),
            ],
        );
    }

    add_kp20_section42_vector_rows(
        rows,
        map,
        &[
            section42_term(Kp20Section42PaperSymbol::P, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, -w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, w1),
        ],
    );
    add_kp20_section42_vector_rows(
        rows,
        map,
        &[
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, -w1),
        ],
    );
}

fn add_kp20_section42_oriented_eq13_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42OrientedPaperSymbolMap,
    b0_near_origin: f64,
) {
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    // KP20 eq. (13): w3 uses the underlined near-m coefficient b0 = c,
    // while w1/w2 use the near-o coefficient b0 = 2c.
    let b0_near_midpoint = 0.5 * b0_near_origin;
    let w3 = b0_near_midpoint / 30.0;

    for i in [5, 6] {
        add_kp20_section42_oriented_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, -0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, 0.25),
            ],
        );
        add_kp20_section42_oriented_vector_rows(
            rows,
            map,
            &[
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 1, 1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 0, -1.0),
                section42_term(Kp20Section42PaperSymbol::UnderlinedP, i, 2, 0.25),
                section42_term(Kp20Section42PaperSymbol::ArrowUnderlinedP, i, 2, -0.25),
            ],
        );
    }

    add_kp20_section42_oriented_vector_rows(
        rows,
        map,
        &[
            section42_term(Kp20Section42PaperSymbol::P, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, -w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, w1),
        ],
    );
    add_kp20_section42_oriented_vector_rows(
        rows,
        map,
        &[
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 1, 1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 0, -1.0),
            section42_term(Kp20Section42PaperSymbol::P, 4, 2, 0.25),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 4, 2, -0.25),
            section42_term(Kp20Section42PaperSymbol::P, 5, 0, w3),
            section42_term(Kp20Section42PaperSymbol::P, 6, 0, -w3),
            section42_term(Kp20Section42PaperSymbol::P, 5, 2, w2),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 5, 2, -w2),
            section42_term(Kp20Section42PaperSymbol::P, 6, 2, w1),
            section42_term(Kp20Section42PaperSymbol::ArrowP, 6, 2, -w1),
        ],
    );
}

const fn section42_term(
    symbol: Kp20Section42PaperSymbol,
    i: usize,
    j: usize,
    weight: f64,
) -> Kp20Section42PaperTerm {
    Kp20Section42PaperTerm {
        symbol,
        i,
        j,
        weight,
    }
}

#[cfg(test)]
fn add_kp20_section42_vector_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    terms: &[Kp20Section42PaperTerm],
) {
    let (sin_theta, cos_theta) = map.previous_sector_rotation_radians.sin_cos();
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for &term in terms {
        let mapped = kp20_section42_mapped_piece(map, term.symbol);
        let coeff = Kp20PlanarTauCoeff::new(mapped.piece, term.i, term.j);
        if mapped.is_previous_sector {
            x_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 0),
                term.weight * cos_theta,
            ));
            x_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 1),
                -term.weight * sin_theta,
            ));
            y_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 0),
                term.weight * sin_theta,
            ));
            y_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 1),
                term.weight * cos_theta,
            ));
        } else {
            x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), term.weight));
            y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), term.weight));
        }
    }
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(x_terms));
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(y_terms));
}

fn add_kp20_fig6c_vector_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42Fig6cSymbolMap,
    terms: &[Kp20Section42PaperTerm],
) {
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for &term in terms {
        let strip = kp20_fig6c_symbol_map(map, term.symbol);
        let (i, j) = strip.stored_index(term.i, term.j);
        let coeff = Kp20PlanarTauCoeff::new(strip.mapped.piece, i, j);
        if strip.sector_rotation_steps != 0 {
            let angle = map.previous_sector_rotation_radians * strip.sector_rotation_steps as f64;
            let (sin_theta, cos_theta) = angle.sin_cos();
            x_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 0),
                term.weight * cos_theta,
            ));
            x_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 1),
                -term.weight * sin_theta,
            ));
            y_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 0),
                term.weight * sin_theta,
            ));
            y_terms.push((
                Kp20PlanarTauScalarCoeff::new(coeff, 1),
                term.weight * cos_theta,
            ));
        } else {
            x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), term.weight));
            y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), term.weight));
        }
    }
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(x_terms));
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(y_terms));
}

fn kp20_fig6c_boundary_row_basis_projected_source_patch(
    patch: &Kp20MacroPatch,
    map: Kp20Section42Fig6cSymbolMap,
    regularization_weight: f64,
) -> (Kp20MacroPatch, f64, f64) {
    let variables = kp20_fig6c_boundary_row_basis_variables(map);
    kp20_fig6c_project_source_patch_with_variables(patch, map, &variables, regularization_weight)
}

fn kp20_fig6c_project_source_patch_with_variables(
    patch: &Kp20MacroPatch,
    map: Kp20Section42Fig6cSymbolMap,
    variables: &[Kp20Fig6cSourceRowVariable],
    regularization_weight: f64,
) -> (Kp20MacroPatch, f64, f64) {
    let mut base = *patch;
    apply_kp20_fig6c_eq12_strip_construction(&mut base, map.p, map.underlined_p, map);
    apply_kp20_fig6c_eq12_strip_construction(&mut base, map.arrow_p, map.arrow_underlined_p, map);
    apply_kp20_fig6c_eq13_strip_construction(&mut base, map);
    apply_kp20_fig6c_eq12_strip_construction(&mut base, map.p, map.underlined_p, map);
    apply_kp20_fig6c_eq12_strip_construction(&mut base, map.arrow_p, map.arrow_underlined_p, map);
    apply_kp20_fig6c_eq13_strip_construction(&mut base, map);

    let rows = kp20_fig6c_full_paper_direction_rows(map);
    let base_residuals = rows
        .iter()
        .map(|row| scalar_residual(&base, row))
        .collect::<Vec<_>>();
    let mut design = vec![vec![0.0_f64; variables.len()]; rows.len()];
    for (variable_index, variable) in variables.iter().copied().enumerate() {
        let mut perturbed = base;
        apply_kp20_fig6c_source_row_variable_delta(&mut perturbed, variable, 1.0, map);
        apply_kp20_fig6c_eq12_strip_construction(&mut perturbed, map.p, map.underlined_p, map);
        apply_kp20_fig6c_eq12_strip_construction(
            &mut perturbed,
            map.arrow_p,
            map.arrow_underlined_p,
            map,
        );
        apply_kp20_fig6c_eq13_strip_construction(&mut perturbed, map);
        for (row_index, row) in rows.iter().enumerate() {
            design[row_index][variable_index] =
                scalar_residual(&perturbed, row) - base_residuals[row_index];
        }
    }
    let rhs = base_residuals
        .iter()
        .map(|residual| -residual)
        .collect::<Vec<_>>();
    let solution = solve_regularized_least_squares(&design, &rhs, regularization_weight);
    let row_lsq_residual = max_dense_row_residual(&design, &solution, &rhs);

    let mut projected = base;
    for (variable, delta) in variables.iter().copied().zip(solution.iter().copied()) {
        apply_kp20_fig6c_source_row_variable_delta(&mut projected, variable, delta, map);
    }
    apply_kp20_fig6c_eq12_strip_construction(&mut projected, map.p, map.underlined_p, map);
    apply_kp20_fig6c_eq12_strip_construction(
        &mut projected,
        map.arrow_p,
        map.arrow_underlined_p,
        map,
    );
    apply_kp20_fig6c_eq13_strip_construction(&mut projected, map);

    let max_source_delta = solution
        .chunks(2)
        .map(|chunk| {
            let x = chunk[0];
            let y = chunk.get(1).copied().unwrap_or(0.0);
            (x * x + y * y).sqrt()
        })
        .fold(0.0_f64, f64::max);
    (projected, max_source_delta, row_lsq_residual)
}

fn kp20_fig6c_boundary_row_basis_variables(
    map: Kp20Section42Fig6cSymbolMap,
) -> Vec<Kp20Fig6cSourceRowVariable> {
    let mut variables = Vec::new();
    for strip in [map.p, map.arrow_p] {
        for paper_i in [2_usize, 4, 6] {
            for axis in [0_usize, 1] {
                variables.push(Kp20Fig6cSourceRowVariable {
                    strip,
                    paper_i,
                    paper_j: 1,
                    axis,
                });
            }
        }
    }
    for strip in [map.p, map.arrow_p, map.underlined_p, map.arrow_underlined_p] {
        for paper_i in 0..=6 {
            for axis in [0_usize, 1] {
                variables.push(Kp20Fig6cSourceRowVariable {
                    strip,
                    paper_i,
                    paper_j: 0,
                    axis,
                });
            }
        }
    }
    variables
}

fn apply_kp20_fig6c_source_row_variable_delta(
    patch: &mut Kp20MacroPatch,
    variable: Kp20Fig6cSourceRowVariable,
    delta: f64,
    map: Kp20Section42Fig6cSymbolMap,
) {
    let mut point = transformed_fig6c_strip_point(
        patch,
        variable.strip,
        variable.paper_i,
        variable.paper_j,
        map,
    );
    if variable.axis == 0 {
        point.x += delta;
    } else {
        point.y += delta;
    }
    set_transformed_fig6c_strip_point(
        patch,
        variable.strip,
        variable.paper_i,
        variable.paper_j,
        point,
        map,
    );
}

fn add_kp20_section42_oriented_vector_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42OrientedPaperSymbolMap,
    terms: &[Kp20Section42PaperTerm],
) {
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for &term in terms {
        let mapped = kp20_section42_oriented_mapped_piece(map, term.symbol);
        let i = if mapped.reverse_i { 6 - term.i } else { term.i };
        let j = if mapped.reverse_j { 6 - term.j } else { term.j };
        let coeff = if mapped.transpose {
            Kp20PlanarTauCoeff::new(mapped.piece, j, i)
        } else {
            Kp20PlanarTauCoeff::new(mapped.piece, i, j)
        };
        x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), term.weight));
        y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), term.weight));
    }
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(x_terms));
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(y_terms));
}

#[cfg(test)]
fn kp20_section42_mapped_piece(
    map: Kp20Section42PaperSymbolMap,
    symbol: Kp20Section42PaperSymbol,
) -> Kp20Section42MappedPiece {
    match symbol {
        Kp20Section42PaperSymbol::P => map.p,
        Kp20Section42PaperSymbol::ArrowP => map.arrow_p,
        Kp20Section42PaperSymbol::UnderlinedP => map.underlined_p,
        Kp20Section42PaperSymbol::ArrowUnderlinedP => map.arrow_underlined_p,
    }
}

fn kp20_section42_oriented_mapped_piece(
    map: Kp20Section42OrientedPaperSymbolMap,
    symbol: Kp20Section42PaperSymbol,
) -> Kp20Section42OrientedStripMap {
    match symbol {
        Kp20Section42PaperSymbol::P => map.p,
        Kp20Section42PaperSymbol::ArrowP => map.arrow_p,
        Kp20Section42PaperSymbol::UnderlinedP => map.underlined_p,
        Kp20Section42PaperSymbol::ArrowUnderlinedP => map.arrow_underlined_p,
    }
}

fn kp20_fig6c_symbol_map(
    map: Kp20Section42Fig6cSymbolMap,
    symbol: Kp20Section42PaperSymbol,
) -> Kp20Section42Fig6cStripMap {
    match symbol {
        Kp20Section42PaperSymbol::P => map.p,
        Kp20Section42PaperSymbol::ArrowP => map.arrow_p,
        Kp20Section42PaperSymbol::UnderlinedP => map.underlined_p,
        Kp20Section42PaperSymbol::ArrowUnderlinedP => map.arrow_underlined_p,
    }
}

impl Kp20Section42Fig6cStripMap {
    fn stored_index(self, paper_i: usize, paper_j: usize) -> (usize, usize) {
        match self.edge {
            Kp20Section42StripEdge::I0 => (paper_j, paper_i),
            Kp20Section42StripEdge::J0 => (paper_i, paper_j),
        }
    }
}

fn transformed_fig6c_strip_point(
    patch: &Kp20MacroPatch,
    strip: Kp20Section42Fig6cStripMap,
    paper_i: usize,
    paper_j: usize,
    map: Kp20Section42Fig6cSymbolMap,
) -> Point3 {
    let (i, j) = strip.stored_index(paper_i, paper_j);
    let cp = patch.piece(strip.mapped.piece).cps[i][j];
    let sector_rotation_steps = if strip.mapped.is_previous_sector {
        strip.sector_rotation_steps
    } else {
        0
    };
    if sector_rotation_steps != 0 {
        rotate_xy(
            cp,
            map.previous_sector_rotation_radians * sector_rotation_steps as f64,
        )
    } else {
        cp
    }
}

fn set_transformed_fig6c_strip_point(
    patch: &mut Kp20MacroPatch,
    strip: Kp20Section42Fig6cStripMap,
    paper_i: usize,
    paper_j: usize,
    transformed_point: Point3,
    map: Kp20Section42Fig6cSymbolMap,
) {
    let sector_rotation_steps = if strip.mapped.is_previous_sector {
        strip.sector_rotation_steps
    } else {
        0
    };
    let stored_point = if sector_rotation_steps != 0 {
        rotate_xy(
            transformed_point,
            -map.previous_sector_rotation_radians * sector_rotation_steps as f64,
        )
    } else {
        transformed_point
    };
    let (i, j) = strip.stored_index(paper_i, paper_j);
    patch.piece_mut(strip.mapped.piece).cps[i][j] = stored_point;
}

fn apply_kp20_fig6c_eq12_strip_construction(
    patch: &mut Kp20MacroPatch,
    source: Kp20Section42Fig6cStripMap,
    receiving: Kp20Section42Fig6cStripMap,
    map: Kp20Section42Fig6cSymbolMap,
) {
    for j in 0..=2 {
        let row6 = transformed_fig6c_strip_point(patch, source, 6, j, map);
        set_transformed_fig6c_strip_point(patch, receiving, 0, j, row6, map);
        if j >= 1 {
            let row4 = transformed_fig6c_strip_point(patch, source, 4, j, map);
            let row2 = transformed_fig6c_strip_point(patch, source, 2, j, map);
            let quarter_delta = scale_point(sub_point(row4, row2), 0.25);
            set_transformed_fig6c_strip_point(
                patch,
                receiving,
                5,
                j,
                add_point(row6, quarter_delta),
                map,
            );
            set_transformed_fig6c_strip_point(
                patch,
                receiving,
                1,
                j,
                sub_point(row6, quarter_delta),
                map,
            );
        }
    }
}

fn apply_kp20_fig6c_eq13_strip_construction(
    patch: &mut Kp20MacroPatch,
    map: Kp20Section42Fig6cSymbolMap,
) {
    let b0_near_origin = 2.0 * map.previous_sector_rotation_radians.cos();
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = (0.5 * b0_near_origin) / 30.0;

    for i in [5_usize, 6] {
        let underlined_i0 = transformed_fig6c_strip_point(patch, map.underlined_p, i, 0, map);
        let underlined_i2 = transformed_fig6c_strip_point(patch, map.underlined_p, i, 2, map);
        let arrow_underlined_i2 =
            transformed_fig6c_strip_point(patch, map.arrow_underlined_p, i, 2, map);
        let quarter_delta = scale_point(sub_point(underlined_i2, arrow_underlined_i2), 0.25);
        set_transformed_fig6c_strip_point(
            patch,
            map.underlined_p,
            i,
            1,
            add_point(underlined_i0, quarter_delta),
            map,
        );
        set_transformed_fig6c_strip_point(
            patch,
            map.arrow_underlined_p,
            i,
            1,
            sub_point(underlined_i0, quarter_delta),
            map,
        );
    }

    let p40 = transformed_fig6c_strip_point(patch, map.p, 4, 0, map);
    let p42 = transformed_fig6c_strip_point(patch, map.p, 4, 2, map);
    let arrow_p42 = transformed_fig6c_strip_point(patch, map.arrow_p, 4, 2, map);
    let p50 = transformed_fig6c_strip_point(patch, map.p, 5, 0, map);
    let p60 = transformed_fig6c_strip_point(patch, map.p, 6, 0, map);
    let p52 = transformed_fig6c_strip_point(patch, map.p, 5, 2, map);
    let arrow_p52 = transformed_fig6c_strip_point(patch, map.arrow_p, 5, 2, map);
    let p62 = transformed_fig6c_strip_point(patch, map.p, 6, 2, map);
    let arrow_p62 = transformed_fig6c_strip_point(patch, map.arrow_p, 6, 2, map);

    let quarter = scale_point(sub_point(p42, arrow_p42), 0.25);
    let w3_term = scale_point(sub_point(p50, p60), w3);
    let w2_term = scale_point(sub_point(p52, arrow_p52), w2);
    let w1_term = scale_point(sub_point(p62, arrow_p62), w1);
    let correction = add_point(sub_point(add_point(quarter, w2_term), w3_term), w1_term);

    set_transformed_fig6c_strip_point(patch, map.p, 4, 1, add_point(p40, correction), map);
    set_transformed_fig6c_strip_point(patch, map.arrow_p, 4, 1, sub_point(p40, correction), map);
}

fn kp20_fig6c_active_set_projection(
    patch: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauScalarLinearConstraint],
    valence: usize,
    sector: usize,
    samples: &[f64],
    regularization_weight: f64,
) -> Result<Kp20Fig6cActiveSetProjectionReport, Kp20PlanarTauPaperSolveError> {
    let mut projected = *patch;
    let mut total_max_delta = 0.0_f64;
    let mut last_active_domain_constraints = 0_usize;
    let mut last_active_orientation_constraints = 0_usize;
    let mut last_linear_residual = 0.0_f64;
    const MAX_ACTIVE_SET_ITERATIONS: usize = 1024;
    for iteration in 0..MAX_ACTIVE_SET_ITERATIONS {
        let active = kp20_guide_domain_active_constraints(&projected, valence, sector, samples)?;
        let orientation_active =
            kp20_orientation_active_constraints(&projected, valence, sector, samples)?;
        last_active_domain_constraints = active.len();
        last_active_orientation_constraints = orientation_active.len();
        if active.is_empty()
            && orientation_active.is_empty()
            && scalar_row_set_residual(&projected, rows) < 1.0e-12
        {
            return Ok(Kp20Fig6cActiveSetProjectionReport {
                tau_macro_patch: projected,
                iterations: iteration,
                active_domain_constraints: 0,
                active_orientation_constraints: 0,
                max_step_delta: total_max_delta,
                linear_residual: last_linear_residual,
            });
        }
        let (next, max_delta, linear_residual) = solve_kp20_fig6c_active_set_step(
            &projected,
            valence,
            sector,
            &active,
            &orientation_active,
            rows,
            regularization_weight,
        )?;
        projected = next;
        total_max_delta = total_max_delta.max(max_delta);
        last_linear_residual = linear_residual;
    }
    Ok(Kp20Fig6cActiveSetProjectionReport {
        tau_macro_patch: projected,
        iterations: MAX_ACTIVE_SET_ITERATIONS,
        active_domain_constraints: last_active_domain_constraints,
        active_orientation_constraints: last_active_orientation_constraints,
        max_step_delta: total_max_delta,
        linear_residual: last_linear_residual,
    })
}

// Index form: `global_axis` addresses the same component of both sector-local
// axis vectors and feeds the flattened variable-index map; keep indices explicit.
#[allow(clippy::needless_range_loop)]
fn solve_kp20_fig6c_active_set_step(
    patch: &Kp20MacroPatch,
    valence: usize,
    sector: usize,
    active: &[Kp20Fig6cDomainActiveConstraint],
    orientation_active: &[Kp20Fig6cOrientationActiveConstraint],
    rows: &[Kp20PlanarTauScalarLinearConstraint],
    regularization_weight: f64,
) -> Result<(Kp20MacroPatch, f64, f64), Kp20PlanarTauPaperSolveError> {
    const ROW_WEIGHT: f64 = 100.0;
    const MAX_CONTROL_STEP: f64 = 0.25;
    let variables = all_macro_patch_axis_variables();
    let row_count = rows.len() + active.len() + orientation_active.len();
    let mut design = vec![vec![0.0_f64; variables.len()]; row_count];
    let mut rhs = vec![0.0_f64; row_count];
    for (row_index, row) in rows.iter().enumerate() {
        rhs[row_index] = -ROW_WEIGHT * scalar_residual(patch, row);
        for &(scalar_coeff, weight) in &row.terms {
            let variable_index =
                macro_patch_axis_variable_index(scalar_coeff.coeff, scalar_coeff.axis);
            design[row_index][variable_index] += ROW_WEIGHT * weight;
        }
    }
    let local_axes = sector_local_axis_components(valence, sector)?;
    let local_patch = global_planar_patch_to_local(valence, sector, *patch)?;
    for (active_index, active_constraint) in active.iter().copied().enumerate() {
        let row_index = rows.len() + active_index;
        let current_value = active_constraint.local_axis_weights[0]
            * sample_local_patch_axis(
                &local_patch,
                active_constraint.piece_id,
                active_constraint.u,
                active_constraint.v,
                0,
            )
            + active_constraint.local_axis_weights[1]
                * sample_local_patch_axis(
                    &local_patch,
                    active_constraint.piece_id,
                    active_constraint.u,
                    active_constraint.v,
                    1,
                );
        rhs[row_index] = active_constraint.target - current_value;
        let bu = bernstein_bi6_all(active_constraint.u);
        let bv = bernstein_bi6_all(active_constraint.v);
        for (i, &bu_i) in bu.iter().enumerate() {
            for (j, &bv_j) in bv.iter().enumerate() {
                let coeff = Kp20PlanarTauCoeff::new(active_constraint.piece_id, i, j);
                let bernstein_weight = bu_i * bv_j;
                for global_axis in 0..=1 {
                    let variable_index = macro_patch_axis_variable_index(coeff, global_axis);
                    design[row_index][variable_index] += bernstein_weight
                        * (active_constraint.local_axis_weights[0] * local_axes[0][global_axis]
                            + active_constraint.local_axis_weights[1] * local_axes[1][global_axis]);
                }
            }
        }
    }
    let orientation_start = rows.len() + active.len();
    for (active_index, active_constraint) in orientation_active.iter().copied().enumerate() {
        let row_index = orientation_start + active_index;
        let local_piece = local_patch.piece(active_constraint.piece_id);
        let (du, dv) =
            sample_local_patch_derivatives(local_piece, active_constraint.u, active_constraint.v);
        let det = du.x * dv.y - du.y * dv.x;
        rhs[row_index] = active_constraint.target_det - det;
        add_orientation_active_design_row(
            &mut design[row_index],
            active_constraint,
            du,
            dv,
            &local_axes,
        );
    }
    let solution = solve_regularized_least_squares(&design, &rhs, regularization_weight);
    let linear_residual = max_dense_row_residual(&design, &solution, &rhs);
    let raw_max_delta = solution
        .chunks(2)
        .map(|chunk| {
            let dx = chunk[0];
            let dy = chunk.get(1).copied().unwrap_or(0.0);
            (dx * dx + dy * dy).sqrt()
        })
        .fold(0.0_f64, f64::max);
    let step_scale = if raw_max_delta.is_finite() && raw_max_delta > MAX_CONTROL_STEP {
        MAX_CONTROL_STEP / raw_max_delta
    } else {
        1.0
    };
    let mut projected = *patch;
    for (&(coeff, axis), delta) in variables.iter().zip(solution.iter().copied()) {
        let cp = &mut projected.piece_mut(coeff.piece).cps[coeff.i][coeff.j];
        let delta = step_scale * delta;
        if axis == 0 {
            cp.x += delta;
        } else {
            cp.y += delta;
        }
    }
    let max_delta = solution
        .chunks(2)
        .map(|chunk| {
            let dx = chunk[0];
            let dy = chunk.get(1).copied().unwrap_or(0.0);
            step_scale * (dx * dx + dy * dy).sqrt()
        })
        .fold(0.0_f64, f64::max);
    Ok((projected, max_delta, linear_residual))
}

fn kp20_guide_domain_active_constraints(
    patch: &Kp20MacroPatch,
    valence: usize,
    sector: usize,
    samples: &[f64],
) -> Result<Vec<Kp20Fig6cDomainActiveConstraint>, Kp20PlanarTauPaperSolveError> {
    let local_patch = global_planar_patch_to_local(valence, sector, *patch)?;
    let mut constraints = Vec::new();
    for piece_id in Kp20PieceId::all() {
        for &u in samples {
            for &v in samples {
                let x = sample_local_patch_axis(&local_patch, piece_id, u, v, 0);
                let y = sample_local_patch_axis(&local_patch, piece_id, u, v, 1);
                if x < -APPENDIX_A_FIG6C_GUIDE_DOMAIN_TOLERANCE {
                    constraints.push(Kp20Fig6cDomainActiveConstraint {
                        piece_id,
                        u,
                        v,
                        local_axis_weights: [1.0, 0.0],
                        target: 0.0,
                    });
                }
                if y < -APPENDIX_A_FIG6C_GUIDE_DOMAIN_TOLERANCE {
                    constraints.push(Kp20Fig6cDomainActiveConstraint {
                        piece_id,
                        u,
                        v,
                        local_axis_weights: [0.0, 1.0],
                        target: 0.0,
                    });
                }
                if x + y > 1.0 + APPENDIX_A_FIG6C_GUIDE_DOMAIN_TOLERANCE {
                    constraints.push(Kp20Fig6cDomainActiveConstraint {
                        piece_id,
                        u,
                        v,
                        local_axis_weights: [1.0, 1.0],
                        target: 1.0,
                    });
                }
            }
        }
    }
    Ok(constraints)
}

fn kp20_planar_tau_square_domain_report(
    tau: &Kp20MacroPatch,
    samples: &[f64],
    tolerance: f64,
) -> Kp20PlanarTauSquareDomainReport {
    let mut report = Kp20PlanarTauSquareDomainReport {
        sample_count: 0,
        outside_square_samples: 0,
        min_u: f64::INFINITY,
        min_v: f64::INFINITY,
        max_u: f64::NEG_INFINITY,
        max_v: f64::NEG_INFINITY,
        max_square_excess: 0.0,
        worst_sample: String::new(),
    };
    for piece_id in Kp20PieceId::all() {
        for &u in samples {
            let bu = bernstein_bi6_all(u);
            for &v in samples {
                let bv = bernstein_bi6_all(v);
                let mut x = 0.0_f64;
                let mut y = 0.0_f64;
                for (i, &bu_i) in bu.iter().enumerate() {
                    for (j, &bv_j) in bv.iter().enumerate() {
                        let cp = tau.piece(piece_id).cps[i][j];
                        let weight = bu_i * bv_j;
                        x += weight * cp.x;
                        y += weight * cp.y;
                    }
                }
                let lower_excess = (-x).max(0.0).max((-y).max(0.0));
                let upper_excess = (x - 1.0).max(0.0).max((y - 1.0).max(0.0));
                let excess = lower_excess.max(upper_excess);
                report.sample_count += 1;
                report.outside_square_samples += usize::from(excess > tolerance);
                report.min_u = report.min_u.min(x);
                report.min_v = report.min_v.min(y);
                report.max_u = report.max_u.max(x);
                report.max_v = report.max_v.max(y);
                if excess > report.max_square_excess {
                    report.max_square_excess = excess;
                    report.worst_sample = format!("{piece_id:?}@({u:.2},{v:.2})=({x:.6},{y:.6})");
                }
            }
        }
    }
    report
}

fn kp20_orientation_active_constraints(
    patch: &Kp20MacroPatch,
    valence: usize,
    sector: usize,
    samples: &[f64],
) -> Result<Vec<Kp20Fig6cOrientationActiveConstraint>, Kp20PlanarTauPaperSolveError> {
    let local_patch = global_planar_patch_to_local(valence, sector, *patch)?;
    let mut constraints = Vec::new();
    for piece_id in Kp20PieceId::all() {
        let piece = local_patch.piece(piece_id);
        for &u in samples {
            for &v in samples {
                let (du, dv) = sample_local_patch_derivatives(piece, u, v);
                let det = du.x * dv.y - du.y * dv.x;
                if det < -1.0e-10 {
                    constraints.push(Kp20Fig6cOrientationActiveConstraint {
                        piece_id,
                        u,
                        v,
                        target_det: 1.0e-3,
                    });
                }
            }
        }
    }
    Ok(constraints)
}

// Index form: `global_axis` addresses the same component of both sector-local
// axis vectors and feeds the flattened variable-index map; keep indices explicit.
#[allow(clippy::needless_range_loop)]
fn add_orientation_active_design_row(
    design_row: &mut [f64],
    active_constraint: Kp20Fig6cOrientationActiveConstraint,
    du: Point3,
    dv: Point3,
    local_axes: &[[f64; 2]; 2],
) {
    let bu5 = bernstein_bi5_all(active_constraint.u);
    let bu6 = bernstein_bi6_all(active_constraint.u);
    let bv5 = bernstein_bi5_all(active_constraint.v);
    let bv6 = bernstein_bi6_all(active_constraint.v);
    for i in 0..=6 {
        for j in 0..=6 {
            let mut local_sensitivity = [0.0_f64; 2];
            if i > 0 {
                let du_weight = 6.0 * bu5[i - 1] * bv6[j];
                local_sensitivity[0] += du_weight * dv.y;
                local_sensitivity[1] -= du_weight * dv.x;
            }
            if i < 6 {
                let du_weight = -6.0 * bu5[i] * bv6[j];
                local_sensitivity[0] += du_weight * dv.y;
                local_sensitivity[1] -= du_weight * dv.x;
            }
            if j > 0 {
                let dv_weight = 6.0 * bu6[i] * bv5[j - 1];
                local_sensitivity[0] -= dv_weight * du.y;
                local_sensitivity[1] += dv_weight * du.x;
            }
            if j < 6 {
                let dv_weight = -6.0 * bu6[i] * bv5[j];
                local_sensitivity[0] -= dv_weight * du.y;
                local_sensitivity[1] += dv_weight * du.x;
            }
            let coeff = Kp20PlanarTauCoeff::new(active_constraint.piece_id, i, j);
            for global_axis in 0..=1 {
                let variable_index = macro_patch_axis_variable_index(coeff, global_axis);
                design_row[variable_index] += local_sensitivity[0] * local_axes[0][global_axis]
                    + local_sensitivity[1] * local_axes[1][global_axis];
            }
        }
    }
}

fn sample_local_patch_axis(
    patch: &Kp20MacroPatch,
    piece_id: Kp20PieceId,
    u: f64,
    v: f64,
    local_axis: usize,
) -> f64 {
    let piece = patch.piece(piece_id);
    let bu = bernstein_bi6_all(u);
    let bv = bernstein_bi6_all(v);
    let mut value = 0.0_f64;
    for (i, &bu_i) in bu.iter().enumerate() {
        for (j, &bv_j) in bv.iter().enumerate() {
            let cp = piece.cps[i][j];
            value += bu_i * bv_j * if local_axis == 0 { cp.x } else { cp.y };
        }
    }
    value
}

fn sample_local_patch_derivatives(
    patch: &super::kp20_bi6_cap::Kp20Bi6Patch,
    u: f64,
    v: f64,
) -> (Point3, Point3) {
    let bu5 = bernstein_bi5_all(u);
    let bu6 = bernstein_bi6_all(u);
    let bv5 = bernstein_bi5_all(v);
    let bv6 = bernstein_bi6_all(v);
    let mut du = Point3::new(0.0, 0.0, 0.0);
    let mut dv = Point3::new(0.0, 0.0, 0.0);
    for (i, &bu5_i) in bu5.iter().enumerate() {
        for (j, &bv_j) in bv6.iter().enumerate() {
            let delta = sub_point(patch.cps[i + 1][j], patch.cps[i][j]);
            du = add_point(du, scale_point(delta, 6.0 * bu5_i * bv_j));
        }
    }
    for (i, &bu_i) in bu6.iter().enumerate() {
        for (j, &bv5_j) in bv5.iter().enumerate() {
            let delta = sub_point(patch.cps[i][j + 1], patch.cps[i][j]);
            dv = add_point(dv, scale_point(delta, 6.0 * bu_i * bv5_j));
        }
    }
    (du, dv)
}

fn scalar_row_set_residual(
    patch: &Kp20MacroPatch,
    rows: &[Kp20PlanarTauScalarLinearConstraint],
) -> f64 {
    rows.iter()
        .map(|row| scalar_residual(patch, row).abs())
        .fold(0.0_f64, f64::max)
}

fn scalar_residual(patch: &Kp20MacroPatch, row: &Kp20PlanarTauScalarLinearConstraint) -> f64 {
    let value = row
        .terms
        .iter()
        .map(|&(scalar_coeff, weight)| {
            let cp = patch.piece(scalar_coeff.coeff.piece).cps[scalar_coeff.coeff.i]
                [scalar_coeff.coeff.j];
            weight * if scalar_coeff.axis == 0 { cp.x } else { cp.y }
        })
        .sum::<f64>();
    value - row.target
}

fn all_macro_patch_axis_variables() -> Vec<(Kp20PlanarTauCoeff, usize)> {
    let mut variables = Vec::with_capacity(4 * 7 * 7 * 2);
    for piece_id in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let coeff = Kp20PlanarTauCoeff::new(piece_id, i, j);
                variables.push((coeff, 0));
                variables.push((coeff, 1));
            }
        }
    }
    variables
}

fn macro_patch_axis_variable_index(coeff: Kp20PlanarTauCoeff, axis: usize) -> usize {
    (((coeff.piece.index() * 7 + coeff.i) * 7 + coeff.j) * 2) + axis
}

fn sector_local_axis_components(
    valence: usize,
    sector: usize,
) -> Result<[[f64; 2]; 2], Kp20PlanarTauPaperSolveError> {
    let (local0, local1) = sector_local_spoke_components(valence, sector)?;
    Ok([local0, local1])
}

fn solve_regularized_least_squares(
    design: &[Vec<f64>],
    rhs: &[f64],
    regularization: f64,
) -> Vec<f64> {
    let variable_count = design.first().map(|row| row.len()).unwrap_or(0);
    let mut normal = vec![vec![0.0_f64; variable_count]; variable_count];
    let mut target = vec![0.0_f64; variable_count];
    for (row, &rhs_value) in design.iter().zip(rhs) {
        if !rhs_value.is_finite() || row.iter().any(|value| !value.is_finite()) {
            return vec![0.0; variable_count];
        }
        // Constraint rows are sparse over the ~400 scalar variables:
        // accumulate only the nonzero support (PK-289). Skipped
        // terms are exact ±0.0 products, which cannot change an
        // accumulator's value — only, in principle, the sign of an
        // exactly-zero entry, which no downstream consumer
        // distinguishes.
        let support: Vec<usize> = (0..variable_count).filter(|&i| row[i] != 0.0).collect();
        for &i in &support {
            target[i] += row[i] * rhs_value;
            for &j in &support {
                normal[i][j] += row[i] * row[j];
            }
        }
    }
    for (i, row) in normal.iter_mut().enumerate() {
        row[i] += regularization;
    }
    solve_dense_finite_guarded(normal, target).unwrap_or_else(|| vec![0.0; variable_count])
}

fn max_dense_row_residual(design: &[Vec<f64>], solution: &[f64], rhs: &[f64]) -> f64 {
    design
        .iter()
        .zip(rhs)
        .map(|(row, &target)| {
            let value = row
                .iter()
                .zip(solution)
                .map(|(&coeff, &variable)| coeff * variable)
                .sum::<f64>();
            (value - target).abs()
        })
        .fold(0.0_f64, f64::max)
}

fn rotate_xy(point: Point3, angle: f64) -> Point3 {
    let (sin_theta, cos_theta) = angle.sin_cos();
    Point3::new(
        cos_theta * point.x - sin_theta * point.y,
        sin_theta * point.x + cos_theta * point.y,
        point.z,
    )
}

fn add_point(left: Point3, right: Point3) -> Point3 {
    Point3::new(left.x + right.x, left.y + right.y, left.z + right.z)
}

fn sub_point(left: Point3, right: Point3) -> Point3 {
    Point3::new(left.x - right.x, left.y - right.y, left.z - right.z)
}

fn scale_point(point: Point3, scale: f64) -> Point3 {
    Point3::new(point.x * scale, point.y * scale, point.z * scale)
}

#[cfg(test)]
fn vector_family_residual(
    family: Kp20PlanarTauPaperRowFamily,
    constraints: &[Kp20PlanarTauLinearConstraint],
    tau_macro_patch: &Kp20MacroPatch,
) -> Kp20PlanarTauPaperFamilyResidual {
    let mut out = Kp20PlanarTauPaperFamilyResidual {
        family,
        max_residual: 0.0,
        worst_row_index: 0,
        worst_axis: Some(0),
    };
    for (row_index, row) in constraints.iter().enumerate() {
        let mut value = [0.0_f64; 2];
        for &(coeff, weight) in &row.terms {
            let cp = tau_macro_patch.piece(coeff.piece).cps[coeff.i][coeff.j];
            value[0] += weight * cp.x;
            value[1] += weight * cp.y;
        }
        for (axis, axis_value) in value.iter().enumerate() {
            let residual = (axis_value - row.target[axis]).abs();
            if residual > out.max_residual {
                out.max_residual = residual;
                out.worst_row_index = row_index;
                out.worst_axis = Some(axis);
            }
        }
    }
    out
}

#[cfg(test)]
fn scalar_family_residual(
    family: Kp20PlanarTauPaperRowFamily,
    constraints: &[Kp20PlanarTauScalarLinearConstraint],
    tau_macro_patch: &Kp20MacroPatch,
) -> Kp20PlanarTauPaperFamilyResidual {
    let mut out = Kp20PlanarTauPaperFamilyResidual {
        family,
        max_residual: 0.0,
        worst_row_index: 0,
        worst_axis: None,
    };
    for (row_index, row) in constraints.iter().enumerate() {
        let mut value = 0.0_f64;
        for &(coeff, weight) in &row.terms {
            let cp = tau_macro_patch.piece(coeff.coeff.piece).cps[coeff.coeff.i][coeff.coeff.j];
            let axis_value = if coeff.axis == 0 { cp.x } else { cp.y };
            value += weight * axis_value;
        }
        let residual = (value - row.target).abs();
        if residual > out.max_residual {
            out.max_residual = residual;
            out.worst_row_index = row_index;
        }
    }
    out
}

#[cfg(test)]
#[path = "kp20_bi6_tau_solve_assemble_tests.rs"]
mod tests;
