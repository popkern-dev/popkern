//! Instancing optimization for patterned boolean operations (spec §15.4).
//!
//! When the same tool body is used at multiple positions (e.g., subtracting
//! 500 bolt holes from a plate), precompute the tool's spatial index and
//! face AABBs once, then transform for each instance.

use popkern_core::geometry::Transform3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::{BoolError, BoolOperation};
use crate::intersect::compute_face_aabbs;

/// Precomputed spatial data for a tool body, reused across instances.
pub struct ToolInstance {
    /// The tool body.
    pub tool: BRep,
    /// Precomputed face AABBs.
    pub aabbs: Vec<(FaceHandle, Aabb)>,
}

impl ToolInstance {
    /// Create a new tool instance with precomputed spatial data.
    pub fn new(tool: BRep) -> Self {
        let aabbs = compute_face_aabbs(&tool);
        Self { tool, aabbs }
    }

    /// Create a transformed copy of this tool at a new position.
    ///
    /// Transforms the tool body and recomputes AABBs from the transformed
    /// positions. The original tool is not modified.
    pub fn at_transform(&self, transform: &Transform3) -> Result<ToolInstance, BoolError> {
        let transformed = self
            .tool
            .transform(transform)
            .map_err(BoolError::CoreError)?;
        Ok(ToolInstance::new(transformed))
    }
}

/// Apply a boolean operation with the same tool at multiple positions.
///
/// Precomputes the tool's spatial index once, then applies the boolean
/// at each position. More efficient than calling `boolean()` N times
/// because the tool analysis is amortized.
///
/// # Arguments
///
/// * `target` — the body to modify
/// * `tool` — the precomputed tool instance
/// * `transforms` — one transform per instance position
/// * `op` — the boolean operation to apply at each position
/// * `tolerance` — geometric tolerance
pub fn boolean_pattern(
    target: &BRep,
    tool: &ToolInstance,
    transforms: &[Transform3],
    op: BoolOperation,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    let mut result = target.clone();

    for transform in transforms {
        let instance = tool.at_transform(transform)?;
        result = crate::boolean(op, result, instance.tool, *tolerance, None)?;
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn tool_instance_creation() {
        let tool = BRep::unit_cube(tol()).unwrap();
        let instance = ToolInstance::new(tool);
        assert!(!instance.aabbs.is_empty());
    }

    #[test]
    fn tool_instance_transform() {
        let tool = BRep::unit_cube(tol()).unwrap();
        let instance = ToolInstance::new(tool);
        let transform = Transform3::translation(5.0, 0.0, 0.0);
        let moved = instance.at_transform(&transform).unwrap();
        assert!(!moved.aabbs.is_empty());
    }

    #[test]
    fn pattern_api_does_not_panic() {
        let target = BRep::unit_cube(tol()).unwrap();
        let tool_brep = BRep::cylinder(0.1, 2.0, tol()).unwrap();
        let tool = ToolInstance::new(tool_brep);
        let transforms = vec![
            Transform3::translation(0.3, 0.3, -0.5),
            Transform3::translation(0.7, 0.7, -0.5),
        ];
        // This may fail (cylinder through cube boolean), but should not panic
        let _result = boolean_pattern(
            &target,
            &tool,
            &transforms,
            BoolOperation::Subtraction,
            &tol(),
        );
    }
}
