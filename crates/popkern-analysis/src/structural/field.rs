//! The solved displacement field's value and gradient at a point: the
//! data a submodel's cut faces are held to.

use crate::error::AnalysisError;
use crate::thermal::sample::point_eval;

use super::DisplacementField;

/// Displacement (mm) and its gradient ∂u_i/∂x_j at a point.
pub type ValueAndGradient = ([f64; 3], [[f64; 3]; 3]);

impl DisplacementField {
    /// Displacement and its gradient at `p` from an active cell holding
    /// it, in the first solid that has one; `None` outside every solid.
    pub fn value_and_gradient_at(
        &self,
        p: [f64; 3],
    ) -> Result<Option<ValueAndGradient>, AnalysisError> {
        for (moments, extension, values) in &self.solids {
            let mut u = [0.0; 3];
            let mut g = [[0.0; 3]; 3];
            let mut found = true;
            for (i, v) in values.iter().enumerate() {
                match point_eval(&self.grid, moments, extension, v, p)? {
                    Some((x, grad)) => {
                        u[i] = x;
                        g[i] = grad;
                    }
                    None => found = false,
                }
            }
            if found {
                return Ok(Some((u, g)));
            }
        }
        Ok(None)
    }
}
