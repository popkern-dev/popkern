//! Assembly of the KP20 planar τ constraint set.
//!
//! The low-level F3 solver accepts independent row families. This module
//! combines the paper row families for one sector while keeping the
//! Catmull-Clark characteristic tensor-border as an explicit input.

use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
#[cfg(test)]
use super::kp20_bi6_tau_f3_solve::kp20_planar_tau_diagonal_symmetry_constraints;
#[cfg(test)]
use super::kp20_bi6_tau_f3_solve::kp20_planar_tau_tensor_border_constraints_with_d1_3;
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
};

/// Secondary source-movement penalty for the §4.2 equation (14) projection.
///
/// This measures the F3-priority limit of the projected source row space. The
/// current limit is row-clean and cap-fold-free, but sampled determinant tests
/// show it collapses to a near-zero Jacobian. Do not route this as an accepted
/// nonsingular KP20 reparameterization without an additional paper-supported
/// source normalization/rank condition.
const SECTION42_SOURCE_PROJECTION_SECONDARY_WEIGHT: f64 = 1.0e-12;

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
    InternalC2,
    /// KP20 Fig. 8b indicated coefficient average.
    Fig8IndicatedAverage,
    /// KP20 equation (14) tensor-border rows.
    TensorBorderEq14,
    /// KP20 equations (8)/(9) along the near-origin sector edge.
    SectorG1G2NearOrigin,
    /// KP20 equations (8)/(9) along the near-midpoint sector edge.
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
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Kp20PlanarTauPaperSolveError {
    /// Valence was invalid for the paper row formulas.
    InvalidValence(usize),
    /// The supplied tensor-border could not be split or indexed.
    TensorBorder(String),
    /// The underlying F3 KKT solve failed.
    Solve(Kp20PlanarTauF3SolveError),
}

impl std::fmt::Display for Kp20PlanarTauPaperSolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidValence(n) => write!(f, "invalid KP20 planar tau valence {n}"),
            Self::TensorBorder(message) => write!(f, "{message}"),
            Self::Solve(err) => write!(f, "{err}"),
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
/// boundary. The diagonal rows use the KP20 Fig. 8 four-quadrant layout:
/// quadrants on the sector diagonal are self-symmetric, and the off-diagonal
/// pair is reflected across the 2x2 macro-square diagonal with both local
/// directions reversed.
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
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(kp20_planar_tau_sector_local_diagonal_symmetry_constraints(
        valence, sector,
    )?);
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
    let map = kp20_section42_standard_map(angle);
    let mut vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let mut scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
    scalar_constraints.extend(kp20_planar_tau_sector_local_diagonal_symmetry_constraints(
        valence, sector,
    )?);
    let (projected_tensor_constraints, _projection_report) =
        kp20_project_section42_source_tensor_rows_with_free_coeffs(
            tensor_constraints,
            &vector_constraints,
            &scalar_constraints,
            &free_coeffs,
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
    let map = kp20_section42_standard_map(angle);
    let vector_constraints = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let tensor_constraints =
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
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

/// Solve one KP20 §4.3 planar τ sector using adjacent-sector equation (8)/(9)
/// rows.
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

/// Evaluate the solved macro-patch against each KP20 paper row family.
///
/// This is a diagnostic localization surface: it does not alter the KKT solve,
/// but it keeps failed candidate comparisons tied to source-level constraints.
#[cfg(test)]
pub fn kp20_planar_tau_paper_family_residuals(
    valence: usize,
    sector: usize,
    characteristic_tensor_border: &TensorBorderBi3,
    tau_macro_patch: &Kp20MacroPatch,
) -> Result<Vec<Kp20PlanarTauPaperFamilyResidual>, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    let c = (core::f64::consts::TAU / valence as f64).cos();
    let tensor_border =
        kp20_planar_tau_tensor_border_constraints(characteristic_tensor_border, sector)
            .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    Ok(vec![
        vector_family_residual(
            Kp20PlanarTauPaperRowFamily::InternalC2,
            &kp20_planar_tau_internal_c2_constraints(),
            tau_macro_patch,
        ),
        vector_family_residual(
            Kp20PlanarTauPaperRowFamily::Fig8IndicatedAverage,
            &kp20_planar_tau_fig8_indicated_averaging_constraints(),
            tau_macro_patch,
        ),
        vector_family_residual(
            Kp20PlanarTauPaperRowFamily::TensorBorderEq14,
            &tensor_border,
            tau_macro_patch,
        ),
        vector_family_residual(
            Kp20PlanarTauPaperRowFamily::SectorG1G2NearOrigin,
            &kp20_planar_tau_sector_g1_g2_constraints(
                Kp20PieceId::NearOriginLeft,
                Kp20PieceId::NearOriginRight,
                kp20_b_coeffs_near_o(c),
            ),
            tau_macro_patch,
        ),
        vector_family_residual(
            Kp20PlanarTauPaperRowFamily::SectorG1G2NearMidpoint,
            &kp20_planar_tau_sector_g1_g2_constraints(
                Kp20PieceId::NearMidpointLeft,
                Kp20PieceId::NearMidpointRight,
                kp20_b_coeffs_near_m(c),
            ),
            tau_macro_patch,
        ),
        scalar_family_residual(
            Kp20PlanarTauPaperRowFamily::DiagonalSymmetry,
            &kp20_planar_tau_diagonal_symmetry_constraints(),
            tau_macro_patch,
        ),
    ])
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
    let map = kp20_section42_standard_map(angle);
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let eq12_rows = kp20_section42_only_eq12_rows(map);
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);
    let fig8_rows = kp20_planar_tau_fig8_indicated_averaging_constraints();
    let scalar_constraints = kp20_section42_presolved_strip_basis_rows(map, 2.0 * c);
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

/// Evaluate the adjacent-sector row set with an explicit KP20 eq. (14)
/// `d¹_3` value against each paper family.
#[cfg(test)]
pub fn kp20_planar_tau_rotational_paper_family_residuals_with_d1_3(
    valence: usize,
    sector: usize,
    characteristic_tensor_border: &TensorBorderBi3,
    d1_3: f64,
    tau_macro_patch: &Kp20MacroPatch,
) -> Result<Vec<Kp20PlanarTauPaperFamilyResidual>, Kp20PlanarTauPaperSolveError> {
    if valence < 3 {
        return Err(Kp20PlanarTauPaperSolveError::InvalidValence(valence));
    }
    let angle = core::f64::consts::TAU / valence as f64;
    let c = angle.cos();
    let tensor_border = kp20_planar_tau_tensor_border_constraints_with_d1_3(
        characteristic_tensor_border,
        sector,
        d1_3,
    )
    .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let map = kp20_section42_standard_map(angle);
    let c0_rows = kp20_section42_only_c0_rows(map);
    let eq11_rows = kp20_section42_only_eq11_rows(map);
    let eq12_rows = kp20_section42_only_eq12_rows(map);
    let eq13_rows = kp20_section42_only_eq13_rows(map, 2.0 * c);
    Ok(vec![
        vector_family_residual(
            Kp20PlanarTauPaperRowFamily::Fig8IndicatedAverage,
            &kp20_planar_tau_fig8_indicated_averaging_constraints(),
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
            &kp20_planar_tau_diagonal_symmetry_constraints(),
            tau_macro_patch,
        ),
    ])
}

#[derive(Debug, Clone, Copy)]
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

fn kp20_section42_standard_map(angle: f64) -> Kp20Section42PaperSymbolMap {
    Kp20Section42PaperSymbolMap {
        p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearOriginRight),
        arrow_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearOriginLeft),
        underlined_p: Kp20Section42MappedPiece::previous(Kp20PieceId::NearMidpointLeft),
        arrow_underlined_p: Kp20Section42MappedPiece::current(Kp20PieceId::NearMidpointRight),
        previous_sector_rotation_radians: angle,
    }
}

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

fn add_kp20_section42_eq13_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    map: Kp20Section42PaperSymbolMap,
    b0_near_origin: f64,
) {
    let w1 = -b0_near_origin / 100.0;
    let w2 = -w1;
    let w3 = b0_near_origin / 30.0;

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
