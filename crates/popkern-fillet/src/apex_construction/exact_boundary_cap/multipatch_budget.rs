//! Local budgets for multi-patch exact-boundary cap candidates.
#![allow(dead_code)]

use popkern_core::tolerance::Tolerance;

use crate::apex_construction::exact_boundary_cap::validation::ExactBoundaryValidationBudget;

const DEFAULT_MAX_SUPPORT_NORMAL_DEFICIT: f64 = 5.0e-2;
const DEFAULT_INTERFACE_G1_DEFICIT: f64 = 5.0e-2;
const DEFAULT_MIN_NORMAL_LENGTH: f64 = 1.0e-12;
const DEFAULT_MAX_WIDTH_DIRECTION_CONTROLS: usize = 4;
const DEFAULT_MAX_INTERFACE_ITERATIONS: usize = 16;
const DEFAULT_MAX_INNER_PATCH_ITERATIONS: usize = 16;

/// Local tolerance and bounded-iteration budget for dormant multi-patch work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactBoundaryMultiPatchBudget {
    /// Maximum residual for locked outer strip boundaries.
    pub outer_boundary_residual: f64,
    /// Maximum support-normal deficit accepted along the outer strip band.
    pub support_normal_deficit: f64,
    /// Maximum positional gap across strip/inner interfaces.
    pub interface_g0: f64,
    /// Maximum tangent-plane deficit across strip/inner interfaces.
    pub interface_g1_deficit: f64,
    /// Maximum positional gap where generated spokes meet the pole.
    pub pole_g0: f64,
    /// Minimum sampled normal length before foldover decline.
    pub min_normal_length: f64,
    /// Width-direction control cap for the first strip builder.
    pub max_width_direction_controls: usize,
    /// Fixed iteration cap for interface construction/refinement.
    pub max_interface_iterations: usize,
    /// Fixed iteration cap for inner-patch construction/refinement.
    pub max_inner_patch_iterations: usize,
}

impl ExactBoundaryMultiPatchBudget {
    /// Build the multi-patch budget from the caller's local exact-boundary
    /// tolerance contract.
    pub(crate) fn from_tolerance(tolerance: Tolerance) -> Self {
        Self {
            outer_boundary_residual: tolerance.linear(),
            support_normal_deficit: DEFAULT_MAX_SUPPORT_NORMAL_DEFICIT,
            interface_g0: tolerance.linear(),
            interface_g1_deficit: DEFAULT_INTERFACE_G1_DEFICIT,
            pole_g0: tolerance.linear(),
            min_normal_length: DEFAULT_MIN_NORMAL_LENGTH,
            max_width_direction_controls: DEFAULT_MAX_WIDTH_DIRECTION_CONTROLS,
            max_interface_iterations: DEFAULT_MAX_INTERFACE_ITERATIONS,
            max_inner_patch_iterations: DEFAULT_MAX_INNER_PATCH_ITERATIONS,
        }
    }

    /// Convert the shared validation subset into the current exact-boundary
    /// validation budget.
    pub(crate) fn strip_validation_budget(self) -> ExactBoundaryValidationBudget {
        ExactBoundaryValidationBudget {
            max_outer_boundary_residual: self.outer_boundary_residual,
            max_support_boundary_residual: self.outer_boundary_residual,
            max_pair_boundary_residual: self.outer_boundary_residual,
            max_internal_g0_gap: self.pole_g0,
            max_support_normal_deficit: self.support_normal_deficit,
            min_normal_length: self.min_normal_length,
        }
    }
}

impl Default for ExactBoundaryMultiPatchBudget {
    fn default() -> Self {
        Self::from_tolerance(Tolerance::default())
    }
}

#[cfg(test)]
mod tests {
    use popkern_core::tolerance::Tolerance;

    use super::ExactBoundaryMultiPatchBudget;

    #[test]
    fn budget_defaults_follow_exact_boundary_tolerance_contract() {
        let tolerance = Tolerance::default().with_linear(2.5e-6);
        let budget = ExactBoundaryMultiPatchBudget::from_tolerance(tolerance);

        assert_eq!(budget.outer_boundary_residual, tolerance.linear());
        assert_eq!(budget.interface_g0, tolerance.linear());
        assert_eq!(budget.pole_g0, tolerance.linear());
        assert_eq!(budget.support_normal_deficit, 5.0e-2);
        assert_eq!(budget.interface_g1_deficit, 5.0e-2);
        assert_eq!(budget.min_normal_length, 1.0e-12);
        assert_eq!(budget.max_width_direction_controls, 4);
        assert_eq!(budget.max_interface_iterations, 16);
        assert_eq!(budget.max_inner_patch_iterations, 16);
    }

    #[test]
    fn strip_validation_budget_reuses_matching_exact_boundary_gates() {
        let budget = ExactBoundaryMultiPatchBudget::default();
        let validation = budget.strip_validation_budget();

        assert_eq!(
            validation.max_outer_boundary_residual,
            budget.outer_boundary_residual
        );
        assert_eq!(
            validation.max_support_boundary_residual,
            budget.outer_boundary_residual
        );
        assert_eq!(
            validation.max_pair_boundary_residual,
            budget.outer_boundary_residual
        );
        assert_eq!(validation.max_internal_g0_gap, budget.pole_g0);
        assert_eq!(
            validation.max_support_normal_deficit,
            budget.support_normal_deficit
        );
        assert_eq!(validation.min_normal_length, budget.min_normal_length);
    }
}
