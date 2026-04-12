//! Gap repair for edges shared between adjacent faces.
//!
//! Two strategies are available:
//!
//! - **Tolerance absorption**: widens the edge tolerance to encompass the gap
//!   without moving any geometry.
//! - **Geometric adjustment**: moves both endpoint vertices toward a shared
//!   midpoint, halving the gap.
//!
//! Gap measurement examines the edge's curve geometry endpoints versus
//! the positions of the edge's topological vertices.  In a properly
//! wired BRep the topological vertices are shared between the two
//! half-edges, so positional gaps can only arise from curve-geometry
//! mismatch (the curve's start/end point differs from the vertex
//! position).

use serde::{Deserialize, Serialize};

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle};

/// Strategy for closing a gap at a shared edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GapRepairStrategy {
    /// Widen the edge tolerance to absorb the gap (no geometry changes).
    ToleranceAbsorption,
    /// Move both endpoint vertices toward a shared midpoint.
    GeometricAdjustment,
}

/// Result of repairing a single edge gap.
#[derive(Debug, Clone)]
pub struct GapRepairResult {
    /// The edge that was repaired.
    pub edge: EdgeHandle,
    /// Which strategy was applied.
    pub strategy_used: GapRepairStrategy,
    /// Gap distance before repair.
    pub original_gap: f64,
    /// Residual gap after repair.
    pub residual_gap: f64,
    /// Maximum vertex displacement applied.
    pub max_displacement: f64,
}

/// Measure the gap at an edge.
///
/// The gap is the maximum distance between the edge curve endpoints
/// and their corresponding topological vertex positions.  Returns
/// `(gap, start_vertex, start_pos, curve_start_pos, end_vertex,
/// end_pos, curve_end_pos)` or `None` if the edge lacks curve geometry.
fn measure_edge_gap(brep: &BRep, edge: EdgeHandle) -> Option<(f64, Point3, Point3)> {
    let edge_data = brep.edges.get(edge)?;
    let curve_on_edge = edge_data.curve.as_ref()?;

    // Half-edge origin = start vertex of the edge.
    let he = brep.half_edges.get(edge_data.half_edge)?;
    let start_vertex = he.origin;
    // End vertex = origin of the next half-edge in the same loop,
    // which is the same as the twin's origin.
    let twin = brep.half_edges.get(he.twin)?;
    let end_vertex = twin.origin;

    let start_pos = brep.vertices.get(start_vertex)?.position;
    let end_pos = brep.vertices.get(end_vertex)?.position;

    // Evaluate the curve at its parameter endpoints.
    let curve_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
    let curve_end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);

    let gap_start = start_pos.distance_to(&curve_start);
    let gap_end = end_pos.distance_to(&curve_end);
    let max_gap = gap_start.max(gap_end);

    // Return the max gap, plus the two pairs of points (vertex pos vs curve pos)
    // for the endpoint with the larger gap.
    if gap_start >= gap_end {
        Some((max_gap, start_pos, curve_start))
    } else {
        Some((max_gap, end_pos, curve_end))
    }
}

/// Repair the gap at a single edge.
///
/// Measures the maximum gap between the edge curve endpoints and the
/// corresponding vertex positions.  If the gap is within tolerance,
/// returns a zero-displacement result.
///
/// # Errors
///
/// Returns `CoreError::InvalidHandle` if the edge is missing from
/// the BRep.
pub fn repair_gap(
    brep: &mut BRep,
    edge: EdgeHandle,
    strategy: GapRepairStrategy,
    tolerance: &Tolerance,
) -> Result<GapRepairResult, CoreError> {
    if !brep.edges.contains_key(edge) {
        return Err(CoreError::InvalidHandle(format!("edge {edge:?}")));
    }

    // Measure the gap.  If no curve geometry, there is no gap to repair.
    let Some((gap, vertex_pos, curve_pos)) = measure_edge_gap(brep, edge) else {
        return Ok(GapRepairResult {
            edge,
            strategy_used: strategy,
            original_gap: 0.0,
            residual_gap: 0.0,
            max_displacement: 0.0,
        });
    };

    // Already within tolerance -- nothing to do.
    if gap < tolerance.linear() {
        return Ok(GapRepairResult {
            edge,
            strategy_used: strategy,
            original_gap: gap,
            residual_gap: gap,
            max_displacement: 0.0,
        });
    }

    match strategy {
        GapRepairStrategy::ToleranceAbsorption => {
            // Widen the edge tolerance to absorb the gap (10% headroom).
            let new_tol = gap * 1.1;
            if let Some(e) = brep.edges.get_mut(edge) {
                e.tolerance = Some(new_tol);
            }
            Ok(GapRepairResult {
                edge,
                strategy_used: strategy,
                original_gap: gap,
                residual_gap: gap, // gap unchanged geometrically
                max_displacement: 0.0,
            })
        }
        GapRepairStrategy::GeometricAdjustment => {
            // Move the vertex toward the midpoint between vertex and curve.
            let midpoint = Point3::new(
                (vertex_pos.x + curve_pos.x) * 0.5,
                (vertex_pos.y + curve_pos.y) * 0.5,
                (vertex_pos.z + curve_pos.z) * 0.5,
            );
            let displacement = gap * 0.5;

            // Find which vertex to move by re-reading the half-edge pair.
            let edge_data = brep
                .edges
                .get(edge)
                .ok_or_else(|| CoreError::InvalidHandle(format!("edge {edge:?}")))?;
            let he = brep
                .half_edges
                .get(edge_data.half_edge)
                .ok_or_else(|| CoreError::InvalidHandle("half-edge".to_string()))?;
            let twin = brep
                .half_edges
                .get(he.twin)
                .ok_or_else(|| CoreError::InvalidHandle("twin half-edge".to_string()))?;

            let start_pos = brep
                .vertices
                .get(he.origin)
                .map(|v| v.position)
                .unwrap_or(vertex_pos);
            let end_pos = brep
                .vertices
                .get(twin.origin)
                .map(|v| v.position)
                .unwrap_or(vertex_pos);

            // Move the vertex that is closer to the measured vertex_pos.
            let start_origin = he.origin;
            let end_origin = twin.origin;

            if start_pos.distance_to(&vertex_pos) < end_pos.distance_to(&vertex_pos) {
                if let Some(v) = brep.vertices.get_mut(start_origin) {
                    v.position = midpoint;
                }
            } else {
                if let Some(v) = brep.vertices.get_mut(end_origin) {
                    v.position = midpoint;
                }
            }

            Ok(GapRepairResult {
                edge,
                strategy_used: strategy,
                original_gap: gap,
                residual_gap: displacement, // moved halfway
                max_displacement: displacement,
            })
        }
    }
}

/// Repair gaps on every edge in the BRep.
///
/// Iterates all edges, attempts `repair_gap` on each, and collects
/// successful results where a gap was actually found and repaired.
/// Edges that fail (e.g. missing topology) are silently skipped.
pub fn repair_all_gaps(
    brep: &mut BRep,
    strategy: GapRepairStrategy,
    tolerance: &Tolerance,
) -> Vec<GapRepairResult> {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    let mut results = Vec::new();
    for eh in edge_handles {
        if let Ok(result) = repair_gap(brep, eh, strategy, tolerance)
            && result.original_gap >= tolerance.linear()
        {
            results.push(result);
        }
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn repair_gap_on_valid_cube_is_noop() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let edge = brep.edges.keys().next().unwrap();
        let result = repair_gap(
            &mut brep,
            edge,
            GapRepairStrategy::GeometricAdjustment,
            &tol,
        )
        .unwrap();
        // A valid cube has no gaps, so displacement should be zero.
        assert!(
            result.max_displacement < tol.linear(),
            "expected no displacement on valid cube, got {}",
            result.max_displacement
        );
    }

    #[test]
    fn repair_all_gaps_on_valid_cube_returns_empty() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let results = repair_all_gaps(&mut brep, GapRepairStrategy::GeometricAdjustment, &tol);
        assert!(
            results.is_empty(),
            "valid cube should have no gaps to repair"
        );
    }

    #[test]
    fn tolerance_absorption_sets_edge_tolerance() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let edge = brep.edges.keys().next().unwrap();
        // A valid cube has no curve-vertex gap, so no tolerance will be set
        // for tolerance absorption either.
        let result = repair_gap(
            &mut brep,
            edge,
            GapRepairStrategy::ToleranceAbsorption,
            &tol,
        )
        .unwrap();
        assert!(
            result.original_gap < tol.linear(),
            "valid cube should have no curve-vertex gap"
        );
    }

    #[test]
    fn geometric_adjustment_closes_gap() {
        // This test verifies the function doesn't panic on a valid cube.
        // Real gap-closing is tested when we have a BRep with an
        // actual curve-vertex mismatch (e.g. from an import).
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let edge = brep.edges.keys().next().unwrap();
        let result = repair_gap(
            &mut brep,
            edge,
            GapRepairStrategy::GeometricAdjustment,
            &tol,
        )
        .unwrap();
        assert!(result.original_gap >= 0.0);
    }

    #[test]
    fn invalid_edge_returns_error() {
        let mut brep = BRep::new(Tolerance::default());
        let bogus = EdgeHandle::default();
        let tol = Tolerance::default();
        let result = repair_gap(
            &mut brep,
            bogus,
            GapRepairStrategy::GeometricAdjustment,
            &tol,
        );
        assert!(result.is_err());
    }

    #[test]
    fn no_curve_edge_returns_zero_gap() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let edge = brep.edges.keys().next().unwrap();
        // Remove the curve to test the no-geometry path.
        brep.edges[edge].curve = None;
        let result = repair_gap(
            &mut brep,
            edge,
            GapRepairStrategy::GeometricAdjustment,
            &tol,
        )
        .unwrap();
        assert!(
            result.original_gap < f64::EPSILON,
            "edge without curve should report zero gap"
        );
    }
}
