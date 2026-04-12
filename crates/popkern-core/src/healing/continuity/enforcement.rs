//! Continuity enforcement — adjust geometry to achieve target continuity.
//!
//! Currently supports G0 enforcement (close positional gaps) by moving
//! vertices toward each other at shared edges.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::healing::quality::ContinuityLevel;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle};

/// Result of a continuity enforcement operation.
#[derive(Debug, Clone)]
pub struct EnforcementResult {
    /// Edge that was enforced.
    pub edge: EdgeHandle,
    /// Maximum displacement introduced.
    pub max_displacement: f64,
    /// Continuity level achieved after enforcement.
    pub achieved_level: ContinuityLevel,
}

/// Enforce G0 continuity at a shared edge by adjusting vertex positions.
///
/// Measures the gap between the edge curve endpoints and their topological
/// vertices. If the gap exceeds tolerance, moves vertices toward the curve
/// endpoints. Reports maximum displacement. Returns error if displacement
/// would exceed `max_displacement`.
///
/// # Errors
///
/// Returns [`CoreError::InvalidHandle`] if the edge is missing from the BRep.
/// Returns [`CoreError::Healing`] if the required displacement exceeds
/// `max_displacement`.
pub fn enforce_g0(
    brep: &mut BRep,
    edge: EdgeHandle,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> Result<EnforcementResult, CoreError> {
    let edge_data = brep
        .edges
        .get(edge)
        .ok_or_else(|| CoreError::InvalidHandle(format!("edge {edge:?}")))?;

    // If no curve geometry, nothing to enforce — already trivially G0.
    let curve_on_edge = match &edge_data.curve {
        Some(c) => c.clone(),
        None => {
            return Ok(EnforcementResult {
                edge,
                max_displacement: 0.0,
                achieved_level: ContinuityLevel::G0,
            });
        }
    };

    // Get start vertex via the half-edge.
    let he_handle = edge_data.half_edge;
    let he = brep
        .half_edges
        .get(he_handle)
        .ok_or_else(|| CoreError::InvalidHandle("half-edge".to_string()))?;
    let start_vh = he.origin;
    let twin_handle = he.twin;

    // Get end vertex via the twin half-edge.
    let twin = brep
        .half_edges
        .get(twin_handle)
        .ok_or_else(|| CoreError::InvalidHandle("twin half-edge".to_string()))?;
    let end_vh = twin.origin;

    // Read current vertex positions.
    let start_pos = brep
        .vertices
        .get(start_vh)
        .ok_or_else(|| CoreError::InvalidHandle("start vertex".to_string()))?
        .position;
    let end_pos = brep
        .vertices
        .get(end_vh)
        .ok_or_else(|| CoreError::InvalidHandle("end vertex".to_string()))?
        .position;

    // Evaluate curve endpoints.
    let curve_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
    let curve_end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);

    // Compute gap at each end.
    let gap_start = start_pos.distance_to(&curve_start);
    let gap_end = end_pos.distance_to(&curve_end);

    // If both gaps are within tolerance, already G0.
    if gap_start < tolerance.linear() && gap_end < tolerance.linear() {
        return Ok(EnforcementResult {
            edge,
            max_displacement: 0.0,
            achieved_level: ContinuityLevel::G0,
        });
    }

    // Compute midpoints and displacements for each gap that exceeds tolerance.
    let mut actual_max_displacement = 0.0_f64;

    let start_adjustment = if gap_start >= tolerance.linear() {
        let mid = midpoint(&start_pos, &curve_start);
        let disp = start_pos.distance_to(&mid);
        actual_max_displacement = actual_max_displacement.max(disp);
        Some((start_vh, mid, disp))
    } else {
        None
    };

    let end_adjustment = if gap_end >= tolerance.linear() {
        let mid = midpoint(&end_pos, &curve_end);
        let disp = end_pos.distance_to(&mid);
        actual_max_displacement = actual_max_displacement.max(disp);
        Some((end_vh, mid, disp))
    } else {
        None
    };

    // Check displacement limit before applying any changes.
    if actual_max_displacement > max_displacement {
        return Err(CoreError::Healing(
            crate::error::HealingError::DisplacementExceeded {
                edge,
                displacement: actual_max_displacement,
                limit: max_displacement,
            },
        ));
    }

    // Apply adjustments.
    if let Some((vh, new_pos, _)) = start_adjustment
        && let Some(v) = brep.vertices.get_mut(vh)
    {
        v.position = new_pos;
    }
    if let Some((vh, new_pos, _)) = end_adjustment
        && let Some(v) = brep.vertices.get_mut(vh)
    {
        v.position = new_pos;
    }

    Ok(EnforcementResult {
        edge,
        max_displacement: actual_max_displacement,
        achieved_level: ContinuityLevel::G0,
    })
}

/// Enforce G0 continuity across all edges in a BRep.
///
/// Iterates every edge and attempts [`enforce_g0`]. Edges that fail
/// (e.g. displacement exceeds limit or missing topology) are silently
/// skipped. Returns results only for edges where enforcement succeeded
/// and had a non-zero displacement.
pub fn enforce_g0_all(
    brep: &mut BRep,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> Vec<EnforcementResult> {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    let mut results = Vec::new();
    for eh in edge_handles {
        if let Ok(result) = enforce_g0(brep, eh, tolerance, max_displacement)
            && result.max_displacement > 0.0
        {
            results.push(result);
        }
    }
    results
}

/// Compute the midpoint of two 3D points.
fn midpoint(a: &Point3, b: &Point3) -> Point3 {
    Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn enforce_g0_on_valid_cube_is_noop() {
        let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
        let tol = Tolerance::default();
        let edge = brep.edges.keys().next().unwrap();
        let result = enforce_g0(&mut brep, edge, &tol, 1.0).unwrap();
        assert!(
            result.max_displacement < f64::EPSILON,
            "valid cube should need no displacement, got {}",
            result.max_displacement
        );
        assert_eq!(result.achieved_level, ContinuityLevel::G0);
    }

    #[test]
    fn enforce_g0_all_on_valid_cube_returns_empty() {
        let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
        let tol = Tolerance::default();
        let results = enforce_g0_all(&mut brep, &tol, 1.0);
        assert!(
            results.is_empty(),
            "valid cube should have no edges needing enforcement"
        );
    }

    #[test]
    fn enforce_g0_on_invalid_edge_returns_error() {
        let mut brep = BRep::new(Tolerance::default());
        let bogus = EdgeHandle::default();
        let tol = Tolerance::default();
        let result = enforce_g0(&mut brep, bogus, &tol, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn enforce_g0_no_curve_returns_zero_displacement() {
        let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
        let tol = Tolerance::default();
        let edge = brep.edges.keys().next().unwrap();
        // Remove curve geometry.
        brep.edges[edge].curve = None;
        let result = enforce_g0(&mut brep, edge, &tol, 1.0).unwrap();
        assert!(
            result.max_displacement < f64::EPSILON,
            "edge without curve should report zero displacement"
        );
        assert_eq!(result.achieved_level, ContinuityLevel::G0);
    }

    #[test]
    fn enforce_g0_displacement_limit_exceeded() {
        let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
        let edge = brep.edges.keys().next().unwrap();

        // Deliberately move a vertex to create a gap.
        let he_handle = brep.edges[edge].half_edge;
        let vh = brep.half_edges[he_handle].origin;
        let original = brep.vertices[vh].position;
        brep.vertices[vh].position =
            Point3::new(original.x + 0.5, original.y + 0.5, original.z + 0.5);

        // Use a very small max_displacement limit.
        let tol = Tolerance::default();
        let result = enforce_g0(&mut brep, edge, &tol, 1e-12);
        assert!(
            result.is_err(),
            "should fail when displacement exceeds limit"
        );
    }

    #[test]
    fn enforce_g0_closes_small_gap() {
        let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
        let edge = brep.edges.keys().next().unwrap();

        // Deliberately move a vertex to create a small gap.
        let he_handle = brep.edges[edge].half_edge;
        let vh = brep.half_edges[he_handle].origin;
        let original = brep.vertices[vh].position;
        let offset = 0.01;
        brep.vertices[vh].position = Point3::new(original.x + offset, original.y, original.z);

        let tol = Tolerance::default();
        let result = enforce_g0(&mut brep, edge, &tol, 1.0).unwrap();
        assert!(
            result.max_displacement > 0.0,
            "should have non-zero displacement"
        );
        assert!(
            result.max_displacement <= offset,
            "displacement should not exceed original gap"
        );
        assert_eq!(result.achieved_level, ContinuityLevel::G0);
    }
}
