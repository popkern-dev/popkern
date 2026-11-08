//! Booleans whose result is several bodies (PK-1115): under
//! `NonManifoldPolicy::SplitToManifold` parts that only touch combine
//! into separate manifold parts, one per piece, each valid on its own.

use popkern_bool::{BoolError, BoolOperation, BoolOptions};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::PartError;
use crate::part::{Part, finer};

impl Part {
    /// [`Part::boolean_with`] as a list of parts: one part, except under
    /// `NonManifoldPolicy::SplitToManifold` a result whose manifold
    /// pieces touch, each piece its own part (with both inputs' feature
    /// logs), in the kernel's order.
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn boolean_bodies_with(
        &self,
        other: &Part,
        op: BoolOperation,
        options: &BoolOptions,
    ) -> Result<Vec<Part>, PartError> {
        let name = match op {
            BoolOperation::Union => "union",
            BoolOperation::Subtraction => "subtract",
            BoolOperation::Intersection => "intersect",
        };
        self.combine_into_parts(other, name, |a, b, tolerance| {
            popkern_bool::boolean_bodies(op, a, b, tolerance, None, options)
        })
    }

    /// [`Part::symmetric_difference`] as a list of parts: under
    /// `NonManifoldPolicy::SplitToManifold` the two differences, where
    /// they touch, separate parts.
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn symmetric_difference_bodies(
        &self,
        other: &Part,
        options: &BoolOptions,
    ) -> Result<Vec<Part>, PartError> {
        self.combine_into_parts(other, "symmetric_difference", |a, b, tolerance| {
            popkern_bool::symmetric_difference_bodies(a, b, tolerance, options)
        })
    }

    fn combine_into_parts(
        &self,
        other: &Part,
        name: &'static str,
        run: impl FnOnce(BRep, BRep, Tolerance) -> Result<Vec<BRep>, BoolError>,
    ) -> Result<Vec<Part>, PartError> {
        let tolerance = finer(self.tolerance(), other.tolerance());
        let bodies =
            run(self.brep().clone(), other.brep().clone(), tolerance).map_err(|source| {
                PartError::Boolean {
                    operation: name,
                    source,
                }
            })?;
        bodies
            .into_iter()
            .map(|brep| Part::derived(brep, tolerance, name, &[self, other]))
            .collect()
    }
}
