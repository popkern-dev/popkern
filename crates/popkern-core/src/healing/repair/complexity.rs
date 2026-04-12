//! NURBS complexity reduction: remove unnecessary internal structure
//! from curves and surfaces while staying within a deviation tolerance.
//!
//! Uses the existing knot removal infrastructure to progressively simplify
//! NURBS representations, reducing control point count and improving
//! evaluation performance.

use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::SurfaceEval;
use crate::topology::{BRep, EdgeHandle};

/// Number of parameter samples used to measure curve deviation.
const CURVE_SAMPLE_COUNT: usize = 32;

/// Number of parameter samples per direction used to measure surface deviation.
const SURFACE_SAMPLE_COUNT: usize = 16;

/// Result of reducing a NURBS curve's complexity.
#[derive(Debug, Clone)]
pub struct CurveReductionResult {
    /// Number of control points before reduction.
    pub points_before: usize,
    /// Number of control points after reduction.
    pub points_after: usize,
    /// Number of knots removed.
    pub knots_removed: usize,
    /// Maximum deviation from the original curve.
    pub max_deviation: f64,
}

/// Result of reducing a NURBS surface's complexity.
#[derive(Debug, Clone)]
pub struct SurfaceReductionResult {
    /// Control points before: (num_u, num_v).
    pub points_before: (usize, usize),
    /// Control points after: (num_u, num_v).
    pub points_after: (usize, usize),
    /// Number of u-knots removed.
    pub u_knots_removed: usize,
    /// Number of v-knots removed.
    pub v_knots_removed: usize,
    /// Maximum deviation from the original surface.
    pub max_deviation: f64,
}

/// Result of reducing complexity across all entities in a BRep.
#[derive(Debug, Clone)]
pub struct ComplexityReductionResult {
    /// Number of edge curves simplified.
    pub curves_reduced: usize,
    /// Number of face surfaces simplified.
    pub surfaces_reduced: usize,
    /// Total knots removed across all entities.
    pub total_knots_removed: usize,
}

// ---------------------------------------------------------------------------
// Interior knot helpers
// ---------------------------------------------------------------------------

/// Return the distinct interior knot indices for a 1D knot vector.
///
/// Interior knots are those strictly inside the clamped region, i.e.
/// indices in `(degree, knots.len() - degree - 1)` exclusive.
fn interior_knot_indices(degree: usize, knots: &[f64]) -> Vec<usize> {
    let first_interior = degree + 1;
    let last_interior = knots.len().saturating_sub(degree + 1);

    if first_interior >= last_interior {
        return Vec::new();
    }

    // Deduplicate: only keep the first occurrence of each distinct knot value.
    let mut indices = Vec::new();
    let mut prev: Option<f64> = None;
    for (idx, val) in knots
        .iter()
        .enumerate()
        .take(last_interior)
        .skip(first_interior)
    {
        if prev.is_none_or(|p| (val - p).abs() > f64::EPSILON) {
            indices.push(idx);
            prev = Some(*val);
        }
    }
    indices
}

/// Measure the maximum point-wise deviation between two curves by sampling.
fn measure_curve_deviation(original: &NurbsCurve, reduced: &NurbsCurve) -> f64 {
    let (t0, t1) = original.parameter_range();
    let mut max_dev = 0.0_f64;
    for i in 0..CURVE_SAMPLE_COUNT {
        let t = t0 + (t1 - t0) * (i as f64) / ((CURVE_SAMPLE_COUNT - 1) as f64);
        let p_orig = original.point_at(t);
        let p_red = reduced.point_at(t);
        max_dev = max_dev.max(p_orig.distance_to(&p_red));
    }
    max_dev
}

/// Measure the maximum point-wise deviation between two surfaces by sampling.
fn measure_surface_deviation(original: &NurbsSurface, reduced: &NurbsSurface) -> f64 {
    let ((u0, u1), (v0, v1)) = original.parameter_range();
    let mut max_dev = 0.0_f64;
    for i in 0..SURFACE_SAMPLE_COUNT {
        let u = u0 + (u1 - u0) * (i as f64) / ((SURFACE_SAMPLE_COUNT - 1) as f64);
        for j in 0..SURFACE_SAMPLE_COUNT {
            let v = v0 + (v1 - v0) * (j as f64) / ((SURFACE_SAMPLE_COUNT - 1) as f64);
            let p_orig = original.point_at(u, v);
            let p_red = reduced.point_at(u, v);
            max_dev = max_dev.max(p_orig.distance_to(&p_red));
        }
    }
    max_dev
}

// ---------------------------------------------------------------------------
// Internal: curve reduction returning both the reduced curve and metrics
// ---------------------------------------------------------------------------

/// Attempt to reduce a NURBS curve, returning the reduced curve and metrics.
///
/// Returns `None` if no interior knots exist or none could be removed.
fn try_reduce_curve(
    curve: &NurbsCurve,
    tolerance: f64,
) -> Option<(NurbsCurve, CurveReductionResult)> {
    let points_before = curve.control_points.len();
    let interior = interior_knot_indices(curve.degree, &curve.knots);

    if interior.is_empty() {
        return None;
    }

    let original = curve.clone();
    let mut current = curve.clone();
    let mut knots_removed = 0usize;

    // Repeatedly scan for removable knots. Each successful removal
    // restarts the scan since the knot vector has changed.
    loop {
        let candidates = interior_knot_indices(current.degree, &current.knots);
        if candidates.is_empty() {
            break;
        }

        let mut removed_any = false;
        // Try from the last candidate backward so index-shifting
        // from removal doesn't affect untried lower indices.
        for &idx in candidates.iter().rev() {
            if let Some(reduced) = current.remove_knot(idx, tolerance) {
                current = reduced;
                knots_removed += 1;
                removed_any = true;
                break;
            }
        }

        if !removed_any {
            break;
        }
    }

    if knots_removed == 0 {
        return None;
    }

    let max_deviation = measure_curve_deviation(&original, &current);
    let points_after = current.control_points.len();

    let result = CurveReductionResult {
        points_before,
        points_after,
        knots_removed,
        max_deviation,
    };

    Some((current, result))
}

// ---------------------------------------------------------------------------
// Curve complexity reduction
// ---------------------------------------------------------------------------

/// Reduce the complexity of a NURBS curve by removing unnecessary knots.
///
/// Iteratively attempts to remove each interior knot. A knot is removed
/// if the resulting curve stays within `tolerance` of the original.
/// Knots are tried from the end of the knot vector toward the start to
/// avoid index-shifting issues.
///
/// Returns `None` if no knots could be removed (curve is already minimal).
pub fn reduce_curve_complexity(curve: &NurbsCurve, tolerance: f64) -> Option<CurveReductionResult> {
    try_reduce_curve(curve, tolerance).map(|(_curve, result)| result)
}

// ---------------------------------------------------------------------------
// Surface complexity reduction
// ---------------------------------------------------------------------------

/// Reduce the complexity of a NURBS surface by removing unnecessary knots.
///
/// Attempts to remove interior knots in both u and v directions.
/// Each removal must stay within `tolerance` of the original surface.
///
/// Returns `None` if no knots could be removed.
pub fn reduce_surface_complexity(
    surface: &NurbsSurface,
    tolerance: f64,
) -> Option<SurfaceReductionResult> {
    let points_before = (surface.num_u, surface.num_v);
    let original = surface.clone();
    let mut current = surface.clone();
    let mut u_knots_removed = 0usize;
    let mut v_knots_removed = 0usize;

    // --- Remove u-direction interior knots ---
    loop {
        let candidates = interior_knot_indices(current.degree_u, &current.knots_u);
        if candidates.is_empty() {
            break;
        }

        let mut removed_any = false;
        for &idx in candidates.iter().rev() {
            if let Some(reduced) = current.remove_knot_u(idx, tolerance) {
                current = reduced;
                u_knots_removed += 1;
                removed_any = true;
                break;
            }
        }

        if !removed_any {
            break;
        }
    }

    // --- Remove v-direction interior knots ---
    loop {
        let candidates = interior_knot_indices(current.degree_v, &current.knots_v);
        if candidates.is_empty() {
            break;
        }

        let mut removed_any = false;
        for &idx in candidates.iter().rev() {
            if let Some(reduced) = current.remove_knot_v(idx, tolerance) {
                current = reduced;
                v_knots_removed += 1;
                removed_any = true;
                break;
            }
        }

        if !removed_any {
            break;
        }
    }

    let total = u_knots_removed + v_knots_removed;
    if total == 0 {
        return None;
    }

    let max_deviation = measure_surface_deviation(&original, &current);
    let points_after = (current.num_u, current.num_v);

    Some(SurfaceReductionResult {
        points_before,
        points_after,
        u_knots_removed,
        v_knots_removed,
        max_deviation,
    })
}

// ---------------------------------------------------------------------------
// BRep-wide curve reduction
// ---------------------------------------------------------------------------

/// Reduce complexity of all NURBS edge curves in a BRep.
///
/// For each edge with a NURBS curve, attempts knot removal within
/// the given tolerance. Returns aggregate results.
pub fn reduce_all_curves(brep: &mut BRep, tolerance: f64) -> ComplexityReductionResult {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    let mut curves_reduced = 0usize;
    let mut total_knots_removed = 0usize;

    for eh in edge_handles {
        let nurbs = match brep.edges.get(eh).and_then(|e| e.curve.as_ref()) {
            Some(coe) => match &coe.geometry {
                CurveGeometry::Nurbs(n) => n.clone(),
                _ => continue,
            },
            None => continue,
        };

        if let Some((reduced_curve, result)) = try_reduce_curve(&nurbs, tolerance) {
            if let Some(edge) = brep.edges.get_mut(eh)
                && let Some(coe) = edge.curve.as_mut()
            {
                coe.geometry = CurveGeometry::Nurbs(reduced_curve);
            }
            curves_reduced += 1;
            total_knots_removed += result.knots_removed;
        }
    }

    ComplexityReductionResult {
        curves_reduced,
        surfaces_reduced: 0,
        total_knots_removed,
    }
}

#[cfg(test)]
#[path = "complexity_tests.rs"]
mod tests;
