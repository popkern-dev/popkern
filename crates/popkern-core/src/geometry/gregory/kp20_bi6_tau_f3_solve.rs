//! Constrained `F3` minimization for planar KP20 bi-6 τ candidates.
//!
//! KP20 §4.3 chooses the remaining coefficients of the planar
//! reparameterization τ by minimizing
//! `F3 = ∫∫ Σ_{a+b=3} 3!/(a!b!) ||∂u^a ∂v^b τ||²`.
//! This module owns the dense KKT substrate for that minimization.

use super::direct_linalg::{condition_number_dense, solve_dense};
use super::kp_polar_boundary_input::TensorBorderBi3;
#[cfg(test)]
use super::kp20_bi6_cap::split_kp20_tensor_border_with_d1_3;
use super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId, split_kp20_tensor_border};
#[cfg(test)]
use super::kp20_bi6_local_constraint::{
    Kp20Eq9DiagnosticTerms, Kp20Eq9TermSigns, kp20_eq9_jacobian_dense_degree8,
    kp20_eq9_jacobian_dense_with_diagnostic_terms, kp20_eq9_jacobian_dense_with_term_signs,
};
use super::kp20_bi6_local_constraint::{kp20_eq8_jacobian_dense, kp20_eq9_jacobian_dense};
use crate::geometry::Point3;

/// One scalar/vector coefficient in a planar τ macro-patch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Kp20PlanarTauCoeff {
    /// Semantic KP20 2x2 piece.
    pub piece: Kp20PieceId,
    /// Tensor-product BB `u` index.
    pub i: usize,
    /// Tensor-product BB `v` index.
    pub j: usize,
}

impl Kp20PlanarTauCoeff {
    /// Construct a coefficient handle.
    pub const fn new(piece: Kp20PieceId, i: usize, j: usize) -> Self {
        Self { piece, i, j }
    }
}

/// One coordinate axis of a planar τ coefficient.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Kp20PlanarTauScalarCoeff {
    /// Tensor-product coefficient handle.
    pub coeff: Kp20PlanarTauCoeff,
    /// Coordinate axis: `0` for x, `1` for y.
    pub axis: usize,
}

impl Kp20PlanarTauScalarCoeff {
    /// Construct an axis-specific coefficient handle.
    pub const fn new(coeff: Kp20PlanarTauCoeff, axis: usize) -> Self {
        Self { coeff, axis }
    }
}

/// Linear vector-valued constraint on planar τ coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauLinearConstraint {
    /// Sparse scalar terms applied independently to x and y.
    pub terms: Vec<(Kp20PlanarTauCoeff, f64)>,
    /// Target `(x, y)` value for the row.
    pub target: [f64; 2],
}

impl Kp20PlanarTauLinearConstraint {
    /// Pin one coefficient to a target value.
    pub fn pin(coeff: Kp20PlanarTauCoeff, target: [f64; 2]) -> Self {
        Self {
            terms: vec![(coeff, 1.0)],
            target,
        }
    }

    /// Construct a homogeneous row.
    pub fn homogeneous(terms: Vec<(Kp20PlanarTauCoeff, f64)>) -> Self {
        Self {
            terms,
            target: [0.0, 0.0],
        }
    }
}

/// Linear scalar constraint on axis-specific planar τ coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauScalarLinearConstraint {
    /// Sparse scalar terms over x/y coefficient axes.
    pub terms: Vec<(Kp20PlanarTauScalarCoeff, f64)>,
    /// Target scalar value for the row.
    pub target: f64,
}

impl Kp20PlanarTauScalarLinearConstraint {
    /// Construct a homogeneous scalar row.
    pub fn homogeneous(terms: Vec<(Kp20PlanarTauScalarCoeff, f64)>) -> Self {
        Self { terms, target: 0.0 }
    }
}

/// Input to the constrained planar τ `F3` solve.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauF3SolveInput {
    /// Initial τ candidate. Used as a tiny regularization target and as the
    /// source for unconstrained `z` values, which are otherwise ignored.
    pub initial: Kp20MacroPatch,
    /// Linear constraints to enforce exactly.
    pub constraints: Vec<Kp20PlanarTauLinearConstraint>,
    /// Tikhonov weight toward the initial candidate. Use a small positive
    /// value to remove the affine nullspace without materially changing the
    /// `F3` optimum.
    pub regularization_weight: f64,
}

/// Input to the constrained planar τ `F3` solve with axis-coupled rows.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauF3ScalarSolveInput {
    /// Initial τ candidate. Used as a tiny regularization target.
    pub initial: Kp20MacroPatch,
    /// Vector-valued rows applied independently to x and y.
    pub vector_constraints: Vec<Kp20PlanarTauLinearConstraint>,
    /// Scalar rows that may couple x and y coefficient axes.
    pub scalar_constraints: Vec<Kp20PlanarTauScalarLinearConstraint>,
    /// Tikhonov weight toward the initial candidate.
    pub regularization_weight: f64,
}

/// Result of a constrained planar τ `F3` solve.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauF3SolveOutput {
    /// Solved planar τ macro-patch.
    pub tau_macro_patch: Kp20MacroPatch,
    /// Maximum absolute constraint residual across x/y rows.
    pub max_constraint_residual: f64,
    /// Condition number of the scalar KKT matrix.
    pub kkt_condition_number: f64,
}

/// One scalar coefficient in an augmented planar τ solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20PlanarTauAugmentedScalarCoeff {
    /// Ordinary planar τ coefficient.
    Tau(Kp20PlanarTauScalarCoeff),
    /// Additional scalar diagnostic unknown.
    Extra(usize),
}

/// Scalar row for an augmented planar τ solve.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauAugmentedScalarLinearConstraint {
    /// Sparse scalar terms over τ coefficients and extra diagnostic unknowns.
    pub terms: Vec<(Kp20PlanarTauAugmentedScalarCoeff, f64)>,
    /// Target scalar value.
    pub target: f64,
}

/// Input to the augmented constrained planar τ `F3` solve.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauAugmentedF3ScalarSolveInput {
    /// Initial τ candidate. Used as a tiny regularization target.
    pub initial: Kp20MacroPatch,
    /// Vector-valued rows applied independently to x and y.
    pub vector_constraints: Vec<Kp20PlanarTauLinearConstraint>,
    /// Scalar rows that may couple x and y coefficient axes.
    pub scalar_constraints: Vec<Kp20PlanarTauScalarLinearConstraint>,
    /// Scalar rows that may also reference extra diagnostic unknowns.
    pub augmented_scalar_constraints: Vec<Kp20PlanarTauAugmentedScalarLinearConstraint>,
    /// Number of extra scalar unknowns.
    pub extra_variable_count: usize,
    /// Tikhonov weight toward the initial τ candidate.
    pub regularization_weight: f64,
    /// Tikhonov weight toward zero for extra diagnostic unknowns.
    pub extra_regularization_weight: f64,
}

/// Output from the augmented constrained planar τ `F3` solve.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PlanarTauAugmentedF3SolveOutput {
    /// Solved planar τ macro-patch.
    pub tau_macro_patch: Kp20MacroPatch,
    /// Solved extra scalar unknowns.
    pub extra_values: Vec<f64>,
    /// Maximum absolute constraint residual across all scalar rows.
    pub max_constraint_residual: f64,
    /// Condition number of the augmented scalar KKT matrix.
    pub kkt_condition_number: f64,
}

/// Test-only rank/nullity report for a planar τ constraint system.
#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Kp20PlanarTauConstraintSpaceReport {
    /// Number of scalar variables in the linear system.
    pub scalar_variable_count: usize,
    /// Number of scalar constraint rows before rank filtering.
    pub dense_row_count: usize,
    /// Number of linearly independent scalar constraint rows.
    pub independent_row_count: usize,
    /// Algebraic nullity of the scalar constraint matrix.
    pub nullity: usize,
}

/// Error from the constrained planar τ `F3` solve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Kp20PlanarTauF3SolveError {
    /// A coefficient handle referenced an invalid tensor index.
    InvalidCoeff(Kp20PlanarTauCoeff),
    /// A scalar coefficient handle referenced an invalid coordinate axis.
    InvalidAxis(usize),
    /// Regularization weight was negative or non-finite.
    InvalidRegularization,
    /// KKT system was singular.
    SingularSystem,
}

impl std::fmt::Display for Kp20PlanarTauF3SolveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidCoeff(coeff) => write!(f, "invalid planar tau coefficient {coeff:?}"),
            Self::InvalidAxis(axis) => write!(f, "invalid planar tau coordinate axis {axis}"),
            Self::InvalidRegularization => write!(f, "invalid planar tau regularization weight"),
            Self::SingularSystem => write!(f, "planar tau F3 KKT system singular"),
        }
    }
}

impl std::error::Error for Kp20PlanarTauF3SolveError {}

/// Minimize KP20 §4.3 `F3` subject to caller-supplied linear constraints.
pub fn solve_kp20_planar_tau_f3(
    input: &Kp20PlanarTauF3SolveInput,
) -> Result<Kp20PlanarTauF3SolveOutput, Kp20PlanarTauF3SolveError> {
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: input.initial,
        vector_constraints: input.constraints.clone(),
        scalar_constraints: Vec::new(),
        regularization_weight: input.regularization_weight,
    })
}

/// Minimize KP20 §4.3 `F3` with rows that may couple x and y axes.
pub fn solve_kp20_planar_tau_f3_with_scalar_constraints(
    input: &Kp20PlanarTauF3ScalarSolveInput,
) -> Result<Kp20PlanarTauF3SolveOutput, Kp20PlanarTauF3SolveError> {
    if !input.regularization_weight.is_finite() || input.regularization_weight < 0.0 {
        return Err(Kp20PlanarTauF3SolveError::InvalidRegularization);
    }
    let constraints =
        dense_scalar_constraints(&input.vector_constraints, &input.scalar_constraints)?;
    let targets = scalar_constraint_targets(&input.vector_constraints, &input.scalar_constraints);
    let active_constraint_indices = independent_constraint_indices(&constraints);
    let active_constraints = active_constraint_indices
        .iter()
        .map(|&idx| constraints[idx].clone())
        .collect::<Vec<_>>();
    let active_targets = active_constraint_indices
        .iter()
        .map(|&idx| targets[idx])
        .collect::<Vec<_>>();
    let objective = scalar_f3_objective_matrix();
    let initial = flatten_scalar_macro_patch(&input.initial);
    let solution = solve_scalar_system(
        &objective,
        &active_constraints,
        &active_targets,
        &initial,
        input.regularization_weight,
    )?;

    let mut tau_macro_patch = input.initial;
    for piece in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let idx = var_index(piece, i, j);
                tau_macro_patch.piece_mut(piece).cps[i][j] =
                    Point3::new(solution.values[idx], solution.values[VARS + idx], 0.0);
            }
        }
    }
    let max_constraint_residual =
        max_scalar_constraint_residual(&constraints, &targets, &solution.values);
    Ok(Kp20PlanarTauF3SolveOutput {
        tau_macro_patch,
        max_constraint_residual,
        kkt_condition_number: solution.kkt_condition_number,
    })
}

/// Minimize KP20 §4.3 `F3` with extra scalar unknowns.
///
/// The live KP20 §4.2 source projection uses this as a lexicographic solve:
/// close the accepted paper rows first, minimize F3, and use source movement
/// only as a secondary selector in the remaining null space.
pub fn solve_kp20_planar_tau_augmented_f3_with_scalar_constraints(
    input: &Kp20PlanarTauAugmentedF3ScalarSolveInput,
) -> Result<Kp20PlanarTauAugmentedF3SolveOutput, Kp20PlanarTauF3SolveError> {
    if !input.regularization_weight.is_finite()
        || input.regularization_weight < 0.0
        || !input.extra_regularization_weight.is_finite()
        || input.extra_regularization_weight < 0.0
    {
        return Err(Kp20PlanarTauF3SolveError::InvalidRegularization);
    }
    let scalar_variable_count = SCALAR_VARS + input.extra_variable_count;
    let mut constraints =
        dense_scalar_constraints(&input.vector_constraints, &input.scalar_constraints)?
            .into_iter()
            .map(|mut row| {
                row.resize(scalar_variable_count, 0.0);
                row
            })
            .collect::<Vec<_>>();
    constraints.extend(dense_augmented_scalar_constraints(
        &input.augmented_scalar_constraints,
        scalar_variable_count,
    )?);
    let mut targets =
        scalar_constraint_targets(&input.vector_constraints, &input.scalar_constraints);
    targets.extend(
        input
            .augmented_scalar_constraints
            .iter()
            .map(|constraint| constraint.target),
    );
    let active_constraint_indices = independent_constraint_indices(&constraints);
    let active_constraints = active_constraint_indices
        .iter()
        .map(|&idx| constraints[idx].clone())
        .collect::<Vec<_>>();
    let active_targets = active_constraint_indices
        .iter()
        .map(|&idx| targets[idx])
        .collect::<Vec<_>>();
    let objective = augmented_scalar_f3_objective_matrix(input.extra_variable_count);
    let mut initial = flatten_scalar_macro_patch(&input.initial);
    initial.resize(scalar_variable_count, 0.0);
    let mut regularization_weights = vec![input.regularization_weight; SCALAR_VARS];
    regularization_weights.resize(scalar_variable_count, input.extra_regularization_weight);
    let solution = solve_general_scalar_system(
        &objective,
        &active_constraints,
        &active_targets,
        &initial,
        &regularization_weights,
    )?;

    let mut tau_macro_patch = input.initial;
    for piece in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let idx = var_index(piece, i, j);
                tau_macro_patch.piece_mut(piece).cps[i][j] =
                    Point3::new(solution.values[idx], solution.values[VARS + idx], 0.0);
            }
        }
    }
    let max_constraint_residual =
        max_scalar_constraint_residual(&constraints, &targets, &solution.values);
    Ok(Kp20PlanarTauAugmentedF3SolveOutput {
        tau_macro_patch,
        extra_values: solution.values[SCALAR_VARS..].to_vec(),
        max_constraint_residual,
        kkt_condition_number: solution.kkt_condition_number,
    })
}

/// Measure the scalar linear free space of a planar τ row system.
#[cfg(test)]
pub fn kp20_planar_tau_constraint_space_report(
    vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
) -> Result<Kp20PlanarTauConstraintSpaceReport, Kp20PlanarTauF3SolveError> {
    let constraints = dense_scalar_constraints(vector_constraints, scalar_constraints)?;
    Ok(constraint_space_report_for_dense_rows(
        SCALAR_VARS,
        &constraints,
    ))
}

/// Measure the scalar linear free space of an augmented planar τ row system.
#[cfg(test)]
#[allow(dead_code)]
pub fn kp20_planar_tau_augmented_constraint_space_report(
    vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
    augmented_scalar_constraints: &[Kp20PlanarTauAugmentedScalarLinearConstraint],
    extra_variable_count: usize,
) -> Result<Kp20PlanarTauConstraintSpaceReport, Kp20PlanarTauF3SolveError> {
    let scalar_variable_count = SCALAR_VARS + extra_variable_count;
    let mut constraints = dense_scalar_constraints(vector_constraints, scalar_constraints)?
        .into_iter()
        .map(|mut row| {
            row.resize(scalar_variable_count, 0.0);
            row
        })
        .collect::<Vec<_>>();
    constraints.extend(dense_augmented_scalar_constraints(
        augmented_scalar_constraints,
        scalar_variable_count,
    )?);
    Ok(constraint_space_report_for_dense_rows(
        scalar_variable_count,
        &constraints,
    ))
}

/// Internal C² constraints for the KP20 Fig. 8 natural 2x2 quadrant layout.
///
/// Piece order is the paper/quadrant order used by [`Kp20PieceId`]:
/// lower-left, lower-right, upper-right, upper-left. The constraints join
/// lower-left to lower-right, upper-left to upper-right, lower-left to
/// upper-left, and lower-right to upper-right through value, first, and
/// second BB difference equality along the shared edge.
pub fn kp20_planar_tau_internal_c2_constraints() -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut rows = Vec::new();
    add_u_c2_join(
        &mut rows,
        Kp20PieceId::NearOriginLeft,
        Kp20PieceId::NearOriginRight,
    );
    add_u_c2_join(
        &mut rows,
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    );
    add_v_c2_join(
        &mut rows,
        Kp20PieceId::NearOriginLeft,
        Kp20PieceId::NearMidpointLeft,
    );
    add_v_c2_join(
        &mut rows,
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointRight,
    );
    rows
}

/// Diagonal symmetry rows for one KP20 Fig. 8 τ sector.
///
/// The sector diagonal swaps domain coordinates and range coordinates. In the
/// Fig. 8 piece layout, pieces 1 and 3 are self-symmetric while pieces 2 and 4
/// mirror each other.
pub fn kp20_planar_tau_diagonal_symmetry_constraints() -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let mut rows = Vec::new();
    add_self_diagonal_symmetry(&mut rows, Kp20PieceId::NearOriginLeft);
    add_self_diagonal_symmetry(&mut rows, Kp20PieceId::NearMidpointRight);
    add_pair_diagonal_symmetry(
        &mut rows,
        Kp20PieceId::NearOriginRight,
        Kp20PieceId::NearMidpointLeft,
    );
    rows
}

/// KP20 Fig. 8b indicated coefficient averaging relations.
///
/// The figure marks the first edge coefficient of each quadrant as the midpoint
/// of its two neighboring edge coefficients before the remaining coefficients
/// minimize `F3`. In the paper quadrant layout this is the bottom-edge row
/// `tau_10 = (tau_00 + tau_20) / 2` of each oriented bi-6 quadrant.
pub fn kp20_planar_tau_fig8_indicated_averaging_constraints() -> Vec<Kp20PlanarTauLinearConstraint>
{
    Kp20PieceId::all()
        .into_iter()
        .map(|piece| {
            Kp20PlanarTauLinearConstraint::homogeneous(vec![
                (Kp20PlanarTauCoeff::new(piece, 1, 0), 1.0),
                (Kp20PlanarTauCoeff::new(piece, 0, 0), -0.5),
                (Kp20PlanarTauCoeff::new(piece, 2, 0), -0.5),
            ])
        })
        .collect()
}

/// KP20 equation (14) tensor-border rows for a planar τ sector.
///
/// The caller supplies the planar tensor-border to which τ should be G²
/// connected at the outer `c,m` boundary. For the paper path this tensor-border
/// is the Catmull-Clark characteristic-ring tensor-border. This helper only
/// encodes the split/reparameterized row constraints from KP20 Fig. 6(a) and
/// equation (14); it intentionally does not synthesize the characteristic ring.
pub fn kp20_planar_tau_tensor_border_constraints(
    tensor_border: &TensorBorderBi3,
    sector: usize,
) -> Result<Vec<Kp20PlanarTauLinearConstraint>, String> {
    let split = split_kp20_tensor_border(tensor_border)?;
    kp20_planar_tau_tensor_border_constraints_from_split(&split, sector)
}

/// KP20 equation (14) tensor-border rows with an explicit `d¹_3`.
#[cfg(test)]
pub fn kp20_planar_tau_tensor_border_constraints_with_d1_3(
    tensor_border: &TensorBorderBi3,
    sector: usize,
    d1_3: f64,
) -> Result<Vec<Kp20PlanarTauLinearConstraint>, String> {
    let split = split_kp20_tensor_border_with_d1_3(tensor_border, d1_3)?;
    kp20_planar_tau_tensor_border_constraints_from_split(&split, sector)
}

/// KP20 equation (14) tensor-border rows from an already-split (and
/// reparameterised) border. Public entry for natively-split sources
/// such as the characteristic ring's two-piece sector boundary.
pub fn kp20_planar_tau_tensor_border_constraints_from_split(
    split: &[super::kp20_bi6_cap::Kp20SplitTensorBorder],
    sector: usize,
) -> Result<Vec<Kp20PlanarTauLinearConstraint>, String> {
    let halves = split.get(sector).ok_or_else(|| {
        format!(
            "KP20 planar tau tensor-border sector {sector} outside valence {}",
            split.len()
        )
    })?;
    let mut rows = Vec::with_capacity(42);
    let previous = (sector + split.len() - 1) % split.len();
    add_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointLeft,
        split[previous].near_c.rows,
    );
    add_tensor_border_half_rows(
        &mut rows,
        Kp20PieceId::NearMidpointRight,
        halves.near_c.rows,
    );
    Ok(rows)
}

/// KP20 §4.1 equation (8)/(9) sector-boundary rows for planar τ pieces.
///
/// `previous_piece` is the paper's `p` and `current_piece` is `p̃`.
/// The shared sector-separating boundary is represented as `v = 0` of both
/// tensor-product pieces, matching the existing KP20 residual extractors.
/// `b_coeffs` must be the degree-2 BB coefficients from KP20 eq. (10):
/// use `[2c, 2c, c]` near `o` and `[c, 0, 0]` near `m`.
pub fn kp20_planar_tau_sector_g1_g2_constraints(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    b_coeffs: [f64; 3],
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let mut rows = Vec::with_capacity(15);
    for jac_row in kp20_eq8_jacobian_dense(b_coeffs) {
        rows.push(jacobian_row_constraint(
            previous_piece,
            current_piece,
            &jac_row,
        ));
    }
    for jac_row in kp20_eq9_jacobian_dense(b_coeffs) {
        rows.push(jacobian_row_constraint(
            previous_piece,
            current_piece,
            &jac_row,
        ));
    }
    rows
}

#[cfg(test)]
/// KP20 §4.1 equation (8)/(9) rows across a sector-separating curve
/// between rotationally adjacent planar τ sectors.
///
/// `previous_piece` is the right-side piece of sector `s - 1`,
/// represented by the canonical sector's coefficient variables rotated by
/// `previous_sector_rotation_radians`. `current_piece` is the left-side piece
/// of sector `s`. This preserves the paper's adjacent-sector meaning while
/// still solving one rotationally symmetric canonical sector.
pub fn kp20_planar_tau_rotated_sector_g1_g2_constraints(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    b_coeffs: [f64; 3],
    previous_sector_rotation_radians: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let (sin_theta, cos_theta) = previous_sector_rotation_radians.sin_cos();
    let mut rows = Vec::with_capacity(30);
    for jac_row in kp20_eq8_jacobian_dense(b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    for jac_row in kp20_eq9_jacobian_dense(b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    rows
}

#[cfg(test)]
/// Variant of [`kp20_planar_tau_rotated_sector_g1_g2_constraints`] that lets
/// diagnostic code orient the eq. (8) and eq. (9) Bernstein `b(u)` rows
/// independently. This is useful because eq. (9) contains `b'(u)`, so a
/// boundary-parameter reversal may expose a G2-only orientation mistake.
pub fn kp20_planar_tau_rotated_sector_g1_g2_constraints_split_b(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    eq8_b_coeffs: [f64; 3],
    eq9_b_coeffs: [f64; 3],
    previous_sector_rotation_radians: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let (sin_theta, cos_theta) = previous_sector_rotation_radians.sin_cos();
    let mut rows = Vec::with_capacity(30);
    for jac_row in kp20_eq8_jacobian_dense(eq8_b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    for jac_row in kp20_eq9_jacobian_dense(eq9_b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    rows
}

#[cfg(test)]
/// Diagnostic variant that keeps eq. (8) fixed but allows eq. (9)'s
/// polynomial terms to use measured sign/orientation variants.
pub fn kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_term_signs(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    eq8_b_coeffs: [f64; 3],
    eq9_b_coeffs: [f64; 3],
    eq9_signs: Kp20Eq9TermSigns,
    previous_sector_rotation_radians: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let (sin_theta, cos_theta) = previous_sector_rotation_radians.sin_cos();
    let mut rows = Vec::with_capacity(30);
    for jac_row in kp20_eq8_jacobian_dense(eq8_b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    for jac_row in kp20_eq9_jacobian_dense_with_term_signs(eq9_b_coeffs, eq9_signs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    rows
}

#[cfg(test)]
/// Diagnostic variant for testing which side supplies the eq. (9) derivative
/// terms under adjacent-sector boundary reversal.
pub fn kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_diagnostic_terms(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    eq8_b_coeffs: [f64; 3],
    eq9_b_coeffs: [f64; 3],
    eq9_terms: Kp20Eq9DiagnosticTerms,
    previous_sector_rotation_radians: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let (sin_theta, cos_theta) = previous_sector_rotation_radians.sin_cos();
    let mut rows = Vec::with_capacity(30);
    for jac_row in kp20_eq8_jacobian_dense(eq8_b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    for jac_row in kp20_eq9_jacobian_dense_with_diagnostic_terms(eq9_b_coeffs, eq9_terms) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    rows
}

#[cfg(test)]
pub fn kp20_planar_tau_rotated_sector_g1_constraints(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    b_coeffs: [f64; 3],
    previous_sector_rotation_radians: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let (sin_theta, cos_theta) = previous_sector_rotation_radians.sin_cos();
    let mut rows = Vec::with_capacity(14);
    for jac_row in kp20_eq8_jacobian_dense(b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    rows
}

#[cfg(test)]
/// Diagnostic variant that enforces eq. (9) as the full degree-8
/// polynomial before KP20's degree-reduction identity is assumed.
pub fn kp20_planar_tau_rotated_sector_g1_g2_constraints_eq9_degree8(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    eq8_b_coeffs: [f64; 3],
    eq9_b_coeffs: [f64; 3],
    previous_sector_rotation_radians: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let (sin_theta, cos_theta) = previous_sector_rotation_radians.sin_cos();
    let mut rows = Vec::with_capacity(32);
    for jac_row in kp20_eq8_jacobian_dense(eq8_b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    for jac_row in kp20_eq9_jacobian_dense_degree8(eq9_b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    rows
}

#[cfg(test)]
pub fn kp20_planar_tau_rotated_sector_g2_constraints(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    b_coeffs: [f64; 3],
    previous_sector_rotation_radians: f64,
) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let (sin_theta, cos_theta) = previous_sector_rotation_radians.sin_cos();
    let mut rows = Vec::with_capacity(16);
    for jac_row in kp20_eq9_jacobian_dense(b_coeffs) {
        add_rotated_jacobian_rows(
            &mut rows,
            previous_piece,
            current_piece,
            &jac_row,
            sin_theta,
            cos_theta,
        );
    }
    rows
}

const VARS_PER_PIECE: usize = 49;
const VARS: usize = 4 * VARS_PER_PIECE;
const SCALAR_VARS: usize = 2 * VARS;

fn sx(piece: Kp20PieceId, i: usize, j: usize) -> Kp20PlanarTauScalarCoeff {
    Kp20PlanarTauScalarCoeff::new(Kp20PlanarTauCoeff::new(piece, i, j), 0)
}

fn sy(piece: Kp20PieceId, i: usize, j: usize) -> Kp20PlanarTauScalarCoeff {
    Kp20PlanarTauScalarCoeff::new(Kp20PlanarTauCoeff::new(piece, i, j), 1)
}

fn add_self_diagonal_symmetry(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    piece: Kp20PieceId,
) {
    for i in 0..=6 {
        rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
            (sx(piece, i, i), 1.0),
            (sy(piece, i, i), -1.0),
        ]));
        for j in (i + 1)..=6 {
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
                (sx(piece, i, j), 1.0),
                (sy(piece, j, i), -1.0),
            ]));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
                (sy(piece, i, j), 1.0),
                (sx(piece, j, i), -1.0),
            ]));
        }
    }
}

fn add_pair_diagonal_symmetry(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    lower_right: Kp20PieceId,
    upper_left: Kp20PieceId,
) {
    for i in 0..=6 {
        for j in 0..=6 {
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
                (sx(lower_right, i, j), 1.0),
                (sy(upper_left, j, i), -1.0),
            ]));
            rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
                (sy(lower_right, i, j), 1.0),
                (sx(upper_left, j, i), -1.0),
            ]));
        }
    }
}

fn add_tensor_border_half_rows(
    rows: &mut Vec<Kp20PlanarTauLinearConstraint>,
    piece: Kp20PieceId,
    half_rows: [[Point3; 7]; 3],
) {
    for (depth, row) in half_rows.iter().enumerate() {
        for (u_index, cp) in row.iter().enumerate() {
            rows.push(Kp20PlanarTauLinearConstraint::pin(
                Kp20PlanarTauCoeff::new(piece, depth, u_index),
                [cp.x, cp.y],
            ));
        }
    }
}

fn jacobian_row_constraint(
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    jacobian_row: &[f64; 98],
) -> Kp20PlanarTauLinearConstraint {
    let mut terms = Vec::new();
    for side in 0..2 {
        let piece = if side == 0 {
            previous_piece
        } else {
            current_piece
        };
        for i in 0..=6 {
            for j in 0..=6 {
                let weight = jacobian_row[side * VARS_PER_PIECE + i * 7 + j];
                if weight != 0.0 {
                    terms.push((Kp20PlanarTauCoeff::new(piece, i, j), weight));
                }
            }
        }
    }
    Kp20PlanarTauLinearConstraint::homogeneous(terms)
}

#[cfg(test)]
fn add_rotated_jacobian_rows(
    rows: &mut Vec<Kp20PlanarTauScalarLinearConstraint>,
    previous_piece: Kp20PieceId,
    current_piece: Kp20PieceId,
    jacobian_row: &[f64; 98],
    sin_theta: f64,
    cos_theta: f64,
) {
    let mut x_terms = Vec::new();
    let mut y_terms = Vec::new();
    for side in 0..2 {
        let piece = if side == 0 {
            previous_piece
        } else {
            current_piece
        };
        for i in 0..=6 {
            for j in 0..=6 {
                let weight = jacobian_row[side * VARS_PER_PIECE + i * 7 + j];
                if weight == 0.0 {
                    continue;
                }
                let coeff = Kp20PlanarTauCoeff::new(piece, i, j);
                if side == 0 {
                    x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), weight * cos_theta));
                    x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), -weight * sin_theta));
                    y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), weight * sin_theta));
                    y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), weight * cos_theta));
                } else {
                    x_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 0), weight));
                    y_terms.push((Kp20PlanarTauScalarCoeff::new(coeff, 1), weight));
                }
            }
        }
    }
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(x_terms));
    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(y_terms));
}

fn add_u_c2_join(
    rows: &mut Vec<Kp20PlanarTauLinearConstraint>,
    left: Kp20PieceId,
    right: Kp20PieceId,
) {
    for j in 0..=6 {
        rows.push(Kp20PlanarTauLinearConstraint::homogeneous(vec![
            (Kp20PlanarTauCoeff::new(left, 6, j), 1.0),
            (Kp20PlanarTauCoeff::new(right, 0, j), -1.0),
        ]));
        rows.push(Kp20PlanarTauLinearConstraint::homogeneous(vec![
            (Kp20PlanarTauCoeff::new(left, 6, j), 1.0),
            (Kp20PlanarTauCoeff::new(left, 5, j), -1.0),
            (Kp20PlanarTauCoeff::new(right, 1, j), -1.0),
            (Kp20PlanarTauCoeff::new(right, 0, j), 1.0),
        ]));
        rows.push(Kp20PlanarTauLinearConstraint::homogeneous(vec![
            (Kp20PlanarTauCoeff::new(left, 6, j), 1.0),
            (Kp20PlanarTauCoeff::new(left, 5, j), -2.0),
            (Kp20PlanarTauCoeff::new(left, 4, j), 1.0),
            (Kp20PlanarTauCoeff::new(right, 2, j), -1.0),
            (Kp20PlanarTauCoeff::new(right, 1, j), 2.0),
            (Kp20PlanarTauCoeff::new(right, 0, j), -1.0),
        ]));
    }
}

fn add_v_c2_join(
    rows: &mut Vec<Kp20PlanarTauLinearConstraint>,
    bottom: Kp20PieceId,
    top: Kp20PieceId,
) {
    for i in 0..=6 {
        rows.push(Kp20PlanarTauLinearConstraint::homogeneous(vec![
            (Kp20PlanarTauCoeff::new(bottom, i, 6), 1.0),
            (Kp20PlanarTauCoeff::new(top, i, 0), -1.0),
        ]));
        rows.push(Kp20PlanarTauLinearConstraint::homogeneous(vec![
            (Kp20PlanarTauCoeff::new(bottom, i, 6), 1.0),
            (Kp20PlanarTauCoeff::new(bottom, i, 5), -1.0),
            (Kp20PlanarTauCoeff::new(top, i, 1), -1.0),
            (Kp20PlanarTauCoeff::new(top, i, 0), 1.0),
        ]));
        rows.push(Kp20PlanarTauLinearConstraint::homogeneous(vec![
            (Kp20PlanarTauCoeff::new(bottom, i, 6), 1.0),
            (Kp20PlanarTauCoeff::new(bottom, i, 5), -2.0),
            (Kp20PlanarTauCoeff::new(bottom, i, 4), 1.0),
            (Kp20PlanarTauCoeff::new(top, i, 2), -1.0),
            (Kp20PlanarTauCoeff::new(top, i, 1), 2.0),
            (Kp20PlanarTauCoeff::new(top, i, 0), -1.0),
        ]));
    }
}

struct ScalarSolution {
    values: Vec<f64>,
    kkt_condition_number: f64,
}

fn solve_scalar_system(
    objective: &[Vec<f64>],
    constraints: &[Vec<f64>],
    targets: &[f64],
    initial: &[f64],
    regularization_weight: f64,
) -> Result<ScalarSolution, Kp20PlanarTauF3SolveError> {
    let regularization_weights = vec![regularization_weight; SCALAR_VARS];
    solve_general_scalar_system(
        objective,
        constraints,
        targets,
        initial,
        &regularization_weights,
    )
}

fn solve_general_scalar_system(
    objective: &[Vec<f64>],
    constraints: &[Vec<f64>],
    targets: &[f64],
    initial: &[f64],
    regularization_weights: &[f64],
) -> Result<ScalarSolution, Kp20PlanarTauF3SolveError> {
    let variable_count = objective.len();
    let size = variable_count + constraints.len();
    let mut a = vec![vec![0.0_f64; size]; size];
    let mut b = vec![0.0_f64; size];
    for r in 0..variable_count {
        for c in 0..variable_count {
            a[r][c] = objective[r][c];
        }
        let weight = regularization_weights[r];
        if weight > 0.0 {
            a[r][r] += weight;
            b[r] += weight * initial[r];
        }
    }
    for (constraint_idx, row) in constraints.iter().enumerate() {
        let kkt_row = variable_count + constraint_idx;
        b[kkt_row] = targets[constraint_idx];
        for (var, &weight) in row.iter().enumerate() {
            if weight != 0.0 {
                a[var][kkt_row] += weight;
                a[kkt_row][var] += weight;
            }
        }
    }
    let kkt_condition_number = condition_number_dense(&a)
        .map(|report| report.condition_number)
        .unwrap_or(f64::INFINITY);
    let values = solve_dense(a, b).ok_or(Kp20PlanarTauF3SolveError::SingularSystem)?
        [..variable_count]
        .to_vec();
    Ok(ScalarSolution {
        values,
        kkt_condition_number,
    })
}

fn dense_augmented_scalar_constraints(
    constraints: &[Kp20PlanarTauAugmentedScalarLinearConstraint],
    scalar_variable_count: usize,
) -> Result<Vec<Vec<f64>>, Kp20PlanarTauF3SolveError> {
    let mut rows = Vec::with_capacity(constraints.len());
    for constraint in constraints {
        let mut row = vec![0.0_f64; scalar_variable_count];
        for &(coeff, weight) in &constraint.terms {
            let index = match coeff {
                Kp20PlanarTauAugmentedScalarCoeff::Tau(coeff) => scalar_var_index_checked(coeff)?,
                Kp20PlanarTauAugmentedScalarCoeff::Extra(extra_index) => {
                    if SCALAR_VARS + extra_index >= scalar_variable_count {
                        return Err(Kp20PlanarTauF3SolveError::InvalidAxis(extra_index));
                    }
                    SCALAR_VARS + extra_index
                }
            };
            row[index] += weight;
        }
        rows.push(row);
    }
    Ok(rows)
}

fn dense_scalar_constraints(
    vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
) -> Result<Vec<Vec<f64>>, Kp20PlanarTauF3SolveError> {
    let mut rows = Vec::with_capacity(2 * vector_constraints.len() + scalar_constraints.len());
    for constraint in vector_constraints {
        for axis in 0..=1 {
            let mut row = vec![0.0_f64; SCALAR_VARS];
            for &(coeff, weight) in &constraint.terms {
                row[scalar_var_index_checked(Kp20PlanarTauScalarCoeff::new(coeff, axis))?] +=
                    weight;
            }
            rows.push(row);
        }
    }
    for constraint in scalar_constraints {
        let mut row = vec![0.0_f64; SCALAR_VARS];
        for &(coeff, weight) in &constraint.terms {
            row[scalar_var_index_checked(coeff)?] += weight;
        }
        rows.push(row);
    }
    Ok(rows)
}

#[cfg(test)]
fn constraint_space_report_for_dense_rows(
    scalar_variable_count: usize,
    constraints: &[Vec<f64>],
) -> Kp20PlanarTauConstraintSpaceReport {
    let independent_row_count = independent_constraint_indices(constraints).len();
    Kp20PlanarTauConstraintSpaceReport {
        scalar_variable_count,
        dense_row_count: constraints.len(),
        independent_row_count,
        nullity: scalar_variable_count.saturating_sub(independent_row_count),
    }
}

fn scalar_constraint_targets(
    vector_constraints: &[Kp20PlanarTauLinearConstraint],
    scalar_constraints: &[Kp20PlanarTauScalarLinearConstraint],
) -> Vec<f64> {
    let mut targets = Vec::with_capacity(2 * vector_constraints.len() + scalar_constraints.len());
    for constraint in vector_constraints {
        targets.push(constraint.target[0]);
        targets.push(constraint.target[1]);
    }
    targets.extend(
        scalar_constraints
            .iter()
            .map(|constraint| constraint.target),
    );
    targets
}

fn independent_constraint_indices(rows: &[Vec<f64>]) -> Vec<usize> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut indices = Vec::new();
    for (row_idx, row) in rows.iter().enumerate() {
        let mut residual = row.clone();
        for basis_row in &basis {
            let projection = dot(&residual, basis_row);
            for (value, basis_value) in residual.iter_mut().zip(basis_row) {
                *value -= projection * basis_value;
            }
        }
        let norm = dot(&residual, &residual).sqrt();
        if norm > 1.0e-10 {
            for value in &mut residual {
                *value /= norm;
            }
            basis.push(residual);
            indices.push(row_idx);
        }
    }
    indices
}

fn scalar_f3_objective_matrix() -> Vec<Vec<f64>> {
    let axis_objective = f3_objective_matrix();
    let mut q = vec![vec![0.0_f64; SCALAR_VARS]; SCALAR_VARS];
    for axis in 0..=1 {
        let offset = axis * VARS;
        for r in 0..VARS {
            for c in 0..VARS {
                q[offset + r][offset + c] = axis_objective[r][c];
            }
        }
    }
    q
}

fn augmented_scalar_f3_objective_matrix(extra_variable_count: usize) -> Vec<Vec<f64>> {
    let mut q = scalar_f3_objective_matrix();
    let variable_count = SCALAR_VARS + extra_variable_count;
    for row in &mut q {
        row.resize(variable_count, 0.0);
    }
    q.resize(variable_count, vec![0.0_f64; variable_count]);
    q
}

fn f3_objective_matrix() -> Vec<Vec<f64>> {
    let mut q = vec![vec![0.0_f64; VARS]; VARS];
    for piece in Kp20PieceId::all() {
        for &(u, wu) in GAUSS7_01 {
            for &(v, wv) in GAUSS7_01 {
                let bu = [
                    bernstein1d_derivative_basis(6, 0, u),
                    bernstein1d_derivative_basis(6, 1, u),
                    bernstein1d_derivative_basis(6, 2, u),
                    bernstein1d_derivative_basis(6, 3, u),
                ];
                let bv = [
                    bernstein1d_derivative_basis(6, 0, v),
                    bernstein1d_derivative_basis(6, 1, v),
                    bernstein1d_derivative_basis(6, 2, v),
                    bernstein1d_derivative_basis(6, 3, v),
                ];
                for (a, bu_a) in bu.iter().enumerate() {
                    let b = 3 - a;
                    let scale = wu * wv * binom(3, a);
                    let mut weights = Vec::with_capacity(VARS_PER_PIECE);
                    for (i, &bu_i) in bu_a.iter().enumerate() {
                        for (j, &bv_j) in bv[b].iter().enumerate() {
                            weights.push((var_index(piece, i, j), bu_i * bv_j));
                        }
                    }
                    for &(row, row_weight) in &weights {
                        for &(col, col_weight) in &weights {
                            q[row][col] += scale * row_weight * col_weight;
                        }
                    }
                }
            }
        }
    }
    q
}

fn flatten_macro_patch(tau: &Kp20MacroPatch) -> Vec<[f64; 2]> {
    let mut out = vec![[0.0_f64; 2]; VARS];
    for piece in Kp20PieceId::all() {
        for i in 0..=6 {
            for j in 0..=6 {
                let cp = tau.piece(piece).cps[i][j];
                out[var_index(piece, i, j)] = [cp.x, cp.y];
            }
        }
    }
    out
}

fn flatten_scalar_macro_patch(tau: &Kp20MacroPatch) -> Vec<f64> {
    let planar = flatten_macro_patch(tau);
    let mut out = vec![0.0_f64; SCALAR_VARS];
    for (idx, value) in planar.iter().enumerate() {
        out[idx] = value[0];
        out[VARS + idx] = value[1];
    }
    out
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn max_scalar_constraint_residual(
    constraints: &[Vec<f64>],
    targets: &[f64],
    solved_values: &[f64],
) -> f64 {
    let mut max_residual = 0.0_f64;
    for (row, target) in constraints.iter().zip(targets) {
        let value = row
            .iter()
            .zip(solved_values)
            .map(|(weight, value)| weight * value)
            .sum::<f64>();
        max_residual = max_residual.max((value - target).abs());
    }
    max_residual
}

fn var_index_checked(coeff: Kp20PlanarTauCoeff) -> Result<usize, Kp20PlanarTauF3SolveError> {
    if coeff.i > 6 || coeff.j > 6 {
        return Err(Kp20PlanarTauF3SolveError::InvalidCoeff(coeff));
    }
    Ok(var_index(coeff.piece, coeff.i, coeff.j))
}

fn scalar_var_index_checked(
    coeff: Kp20PlanarTauScalarCoeff,
) -> Result<usize, Kp20PlanarTauF3SolveError> {
    if coeff.axis > 1 {
        return Err(Kp20PlanarTauF3SolveError::InvalidAxis(coeff.axis));
    }
    Ok(coeff.axis * VARS + var_index_checked(coeff.coeff)?)
}

fn var_index(piece: Kp20PieceId, i: usize, j: usize) -> usize {
    piece.index() * VARS_PER_PIECE + i * 7 + j
}

fn bernstein1d_derivative_basis(degree: usize, order: usize, t: f64) -> [f64; 7] {
    let mut out = [0.0_f64; 7];
    if order > degree {
        return out;
    }
    let scale = falling_factorial(degree, order);
    let lower_degree = degree - order;
    for (i, out_i) in out.iter_mut().enumerate().take(degree + 1) {
        let start = i.saturating_sub(order);
        let end = i.min(lower_degree);
        let mut value = 0.0_f64;
        for r in start..=end {
            let sign = if ((order - (i - r)) & 1) == 0 {
                1.0
            } else {
                -1.0
            };
            value += sign * binom(order, i - r) * bernstein1d(lower_degree, r, t);
        }
        *out_i = scale * value;
    }
    out
}

fn bernstein1d(degree: usize, index: usize, t: f64) -> f64 {
    binom(degree, index) * t.powi(index as i32) * (1.0 - t).powi((degree - index) as i32)
}

fn binom(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

fn falling_factorial(n: usize, order: usize) -> f64 {
    (0..order).fold(1.0, |acc, k| acc * (n - k) as f64)
}

const GAUSS7_01: &[(f64, f64); 7] = &[
    (0.025446043828620737, 0.06474248308443485),
    (0.12923440720030277, 0.13985269574463832),
    (0.2970774243113014, 0.19091502525255946),
    (0.5, 0.2089795918367347),
    (0.7029225756886985, 0.19091502525255946),
    (0.8707655927996972, 0.13985269574463832),
    (0.9745539561713793, 0.06474248308443485),
];

#[cfg(test)]
#[path = "kp20_bi6_tau_f3_solve_tests.rs"]
mod tests;
