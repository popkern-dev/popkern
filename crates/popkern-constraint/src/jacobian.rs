//! Sparse Jacobian assembly for the constraint equation system.
//!
//! Iterates all **driving** constraints in a [`Sketch`], evaluates their
//! residuals and analytic Jacobian entries, and collects everything into
//! a single [`EquationSystem`] ready for the Newton–Raphson solver.

use crate::constraint::ConstraintKind;
use crate::equations::{equation_count, evaluate_constraint};
use crate::Sketch;

/// Assembled equation system: residuals F(x) and sparse Jacobian dF/dx.
pub(crate) struct EquationSystem {
    /// Residual vector F(x); length = `num_equations`.
    pub residuals: Vec<f64>,
    /// Sparse Jacobian entries as `(row, col, value)` triplets.
    pub jacobian_triplets: Vec<(usize, usize, f64)>,
    /// Total number of scalar equations.
    pub num_equations: usize,
    /// Total number of parameters in the sketch.
    pub num_params: usize,
}

/// Assemble the full equation system from all driving constraints.
///
/// Only constraints with [`ConstraintKind::Driving`] are included; driven
/// (display-only) constraints are skipped.
pub(crate) fn assemble(sketch: &Sketch) -> EquationSystem {
    let num_params = sketch.params.len();

    // Pre-count total equations so we can reserve capacity.
    let total_equations: usize = sketch
        .constraints
        .values()
        .filter(|e| e.kind == ConstraintKind::Driving)
        .map(|e| equation_count(&e.constraint))
        .sum();

    let mut residuals = Vec::with_capacity(total_equations);
    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    let mut row_offset: usize = 0;

    for (handle, entry) in sketch.constraints.iter() {
        if entry.kind != ConstraintKind::Driving {
            continue;
        }

        let output = evaluate_constraint(sketch, Some(handle), entry, row_offset);

        residuals.extend_from_slice(&output.residuals);
        triplets.extend_from_slice(&output.jacobian);

        row_offset += equation_count(&entry.constraint);
    }

    EquationSystem {
        residuals,
        jacobian_triplets: triplets,
        num_equations: row_offset,
        num_params,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint::Constraint;

    #[test]
    fn empty_sketch_produces_empty_system() {
        let sketch = Sketch::new();
        let sys = assemble(&sketch);
        assert_eq!(sys.num_equations, 0);
        assert_eq!(sys.num_params, 0);
        assert!(sys.residuals.is_empty());
        assert!(sys.jacobian_triplets.is_empty());
    }

    #[test]
    fn driven_constraints_are_skipped() {
        let mut sketch = Sketch::new();
        let a = sketch.add_point(0.0, 0.0);
        let b = sketch.add_point(1.0, 0.0);
        // Driven constraint — should not appear in the system.
        sketch.add_driven(Constraint::Coincident(a, b));

        let sys = assemble(&sketch);
        assert_eq!(sys.num_equations, 0);
        assert_eq!(sys.num_params, 4);
    }

    #[test]
    fn single_coincident_constraint() {
        let mut sketch = Sketch::new();
        let a = sketch.add_point(1.0, 2.0);
        let b = sketch.add_point(3.0, 4.0);
        sketch.add_constraint(Constraint::Coincident(a, b));

        let sys = assemble(&sketch);
        assert_eq!(sys.num_equations, 2);
        assert_eq!(sys.num_params, 4);
        // Residuals: x1-x2 = 1-3 = -2, y1-y2 = 2-4 = -2
        assert!((sys.residuals[0] - (-2.0)).abs() < 1e-12);
        assert!((sys.residuals[1] - (-2.0)).abs() < 1e-12);
        // Should have 4 Jacobian entries (2 per equation).
        assert_eq!(sys.jacobian_triplets.len(), 4);
    }

    #[test]
    fn multiple_constraints_stacked() {
        let mut sketch = Sketch::new();
        let a = sketch.add_point(0.0, 0.0);
        let b = sketch.add_point(1.0, 1.0);
        let l = sketch.add_line(a, b);
        // Coincident: 2 equations (rows 0,1)
        sketch.add_constraint(Constraint::Coincident(a, b));
        // Horizontal: 1 equation (row 2)
        sketch.add_constraint(Constraint::Horizontal(l));

        let sys = assemble(&sketch);
        assert_eq!(sys.num_equations, 3);
        assert_eq!(sys.residuals.len(), 3);

        // Verify row indices in triplets span 0..3.
        let max_row = sys
            .jacobian_triplets
            .iter()
            .map(|(r, _, _)| *r)
            .max()
            .unwrap_or(0);
        assert_eq!(max_row, 2);
    }
}
