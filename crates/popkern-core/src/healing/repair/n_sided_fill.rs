//! N-sided hole filling by decomposition into 4-sided Coons patches.
//!
//! For holes with more than 4 boundary edges, introduces a central
//! point and decomposes the hole into a fan of 4-sided sub-patches.
//! Each sub-patch is bounded by one boundary edge, two radial lines
//! to the center, and a degenerate top edge collapsed at the center.
//!
//! For N = 3 (triangular holes), the same fan decomposition applies:
//! three triangular sub-patches radiating from the centroid.

use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};
use crate::topology::BRep;

use super::coons_patch::{CoonsBoundary, CoonsPatchResult, construct_coons_patch, point_to_vec};
use super::hole_fill::DetectedHole;

/// Construct a degree-1 NURBS line segment between two points.
fn line_nurbs(start: Point3, end: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    }
}

/// Compute the centroid of a set of points.
fn centroid(points: &[Point3]) -> Option<Point3> {
    if points.is_empty() {
        return None;
    }
    let n = points.len() as f64;
    let sum = points
        .iter()
        .fold(Vec3::zero(), |acc, p| acc + point_to_vec(*p));
    Some(Point3::new(sum.x / n, sum.y / n, sum.z / n))
}

/// Collect boundary vertex positions for a detected hole.
///
/// For each edge, extracts the start point of the edge's curve.
/// Returns the positions in boundary-loop order.
fn hole_vertex_positions(brep: &BRep, hole: &DetectedHole) -> Vec<Point3> {
    let mut positions = Vec::with_capacity(hole.boundary_edges.len());
    for &edge_handle in &hole.boundary_edges {
        let Some(edge_data) = brep.edges.get(edge_handle) else {
            continue;
        };
        let Some(ref curve_on_edge) = edge_data.curve else {
            continue;
        };
        let (t_start, _) = (curve_on_edge.t_start, curve_on_edge.t_end);
        positions.push(curve_on_edge.geometry.point_at(t_start));
    }
    positions
}

/// Fill an N-sided hole by decomposing into 4-sided patches.
///
/// For holes with N >= 3 sides, introduces a central point and
/// splits the hole into N 4-sided sub-patches, each bounded by
/// one hole edge (bottom), two radial lines to the center (left
/// and right), and a degenerate line at the center (top).
///
/// Each sub-patch is filled with a standard Coons patch. The fan
/// arrangement ensures complete coverage of the hole interior.
///
/// Returns `None` if the hole has fewer than 3 sides, any edge
/// lacks curve geometry, or the centroid computation fails.
pub fn fill_n_sided_hole(
    brep: &BRep,
    hole: &DetectedHole,
    tolerance: f64,
) -> Option<Vec<CoonsPatchResult>> {
    if hole.num_sides < 3 {
        return None;
    }

    // Collect boundary curves.
    let mut curves: Vec<NurbsCurve> = Vec::with_capacity(hole.num_sides);
    for &edge_handle in &hole.boundary_edges {
        let edge_data = brep.edges.get(edge_handle)?;
        let curve_on_edge = edge_data.curve.as_ref()?;
        let nurbs = curve_on_edge
            .geometry
            .to_nurbs_in_range((curve_on_edge.t_start, curve_on_edge.t_end));
        curves.push(nurbs);
    }

    // Compute centroid from boundary vertex positions.
    let positions = hole_vertex_positions(brep, hole);
    let center = centroid(&positions)?;

    // For each boundary edge, build a fan sub-patch.
    let n = hole.num_sides;
    let mut patches = Vec::with_capacity(n);

    for curve in &curves {
        // Edge start and end points.
        let (t0, t1) = curve.parameter_range();
        let edge_start = curve.point_at(t0);
        let edge_end = curve.point_at(t1);

        // Small offset from center for the degenerate top edge to avoid
        // a true singularity. Use a point slightly toward the edge midpoint.
        let edge_mid = curve.point_at(t0 + 0.5 * (t1 - t0));
        let toward_edge = point_to_vec(edge_mid) - point_to_vec(center);
        let offset = toward_edge * 0.01;
        let center_near = Point3::new(
            center.x + offset.x,
            center.y + offset.y,
            center.z + offset.z,
        );

        let boundary = CoonsBoundary {
            bottom: curve.clone(),
            left: line_nurbs(edge_start, center),
            right: line_nurbs(edge_end, center),
            top: line_nurbs(center_near, center_near),
        };

        // Use a modest grid; the degenerate top edge limits quality anyway.
        if let Some(patch) = construct_coons_patch(&boundary, 6, 6, tolerance) {
            patches.push(patch);
        }
    }

    if patches.is_empty() {
        return None;
    }

    Some(patches)
}

#[cfg(test)]
#[path = "n_sided_fill_tests.rs"]
mod tests;
