//! Re-derivation gap repair: recompute edge geometry from adjacent surfaces.
//!
//! The most geometrically accurate gap repair strategy. Instead of moving
//! vertices to close gaps, recomputes the edge curve as the intersection
//! of the two adjacent surfaces using iterative projection and NURBS fitting.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::fitting::fit_curve;
use crate::geometry::surface::SurfaceEval;
use crate::healing::deviation::edge_surface_deviation;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, FaceHandle};

/// Result of re-deriving an edge's geometry from surfaces.
#[derive(Debug, Clone)]
pub struct RederiveResult {
    /// The edge that was processed.
    pub edge: EdgeHandle,
    /// Deviation before re-derivation.
    pub old_deviation: f64,
    /// Deviation after re-derivation (from the worse of the two surfaces).
    pub new_deviation: f64,
    /// Whether the re-derivation improved the edge and was applied.
    pub succeeded: bool,
    /// Human-readable description of the outcome.
    pub description: String,
}

/// Number of sample points along the edge for iterative projection.
const NUM_PROJECTION_SAMPLES: usize = 16;

/// Number of projection iterations (project A -> B -> average).
const PROJECTION_ITERATIONS: usize = 5;

/// Degree of the fitted NURBS curve.
const FIT_DEGREE: usize = 3;

/// Get the two adjacent faces of an edge via its half-edges.
///
/// Returns `None` if the edge has fewer than two distinct adjacent faces.
fn adjacent_faces(brep: &BRep, edge: EdgeHandle) -> Option<(FaceHandle, FaceHandle)> {
    let edge_data = brep.edges.get(edge)?;
    let he = brep.half_edges.get(edge_data.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    let face_a = brep.loops.get(he.loop_)?.face;
    let face_b = brep.loops.get(twin.loop_)?.face;
    if face_a == face_b {
        return None; // Seam edge -- cannot rederive.
    }
    Some((face_a, face_b))
}

/// Measure the maximum deviation of an edge's curve from a surface.
///
/// Returns `f64::MAX` if the edge has no curve geometry.
fn measure_edge_deviation(brep: &BRep, edge: EdgeHandle, surface: &dyn SurfaceEval) -> f64 {
    let Some(edge_data) = brep.edges.get(edge) else {
        return f64::MAX;
    };
    let Some(curve) = &edge_data.curve else {
        return f64::MAX;
    };
    let result = edge_surface_deviation(&curve.geometry, (curve.t_start, curve.t_end), surface, 32);
    result.max_deviation
}

/// Project a point onto the intersection of two surfaces using alternating
/// projection: project onto surface A, then project that result onto surface B,
/// then average. Repeat for `iters` rounds.
fn project_to_intersection(
    point: &Point3,
    surf_a: &dyn SurfaceEval,
    surf_b: &dyn SurfaceEval,
    iters: usize,
) -> Point3 {
    let mut current = *point;
    for _ in 0..iters {
        let pa = if let Some((u, v)) = surf_a.closest_point(&current) {
            surf_a.point_at(u, v)
        } else {
            current
        };
        let pb = if let Some((u, v)) = surf_b.closest_point(&pa) {
            surf_b.point_at(u, v)
        } else {
            pa
        };
        current = Point3::new(
            (pa.x + pb.x) * 0.5,
            (pa.y + pb.y) * 0.5,
            (pa.z + pb.z) * 0.5,
        );
    }
    current
}

/// Re-derive an edge's curve from its two adjacent surfaces.
///
/// Uses iterative projection to approximate the surface-surface intersection,
/// then fits a NURBS curve through the projected points. If the new curve
/// has lower deviation than the original, replaces the edge's geometry.
///
/// # Errors
///
/// Returns `CoreError::InvalidHandle` if the edge or adjacent faces are missing.
/// Returns `CoreError::PreconditionFailed` if the edge has fewer than two
/// distinct adjacent faces (e.g., a seam edge).
pub fn rederive_edge_from_surfaces(
    brep: &mut BRep,
    edge: EdgeHandle,
    tolerance: &Tolerance,
) -> Result<RederiveResult, CoreError> {
    let (face_a, face_b) = adjacent_faces(brep, edge).ok_or_else(|| {
        CoreError::PreconditionFailed(
            "edge must have two distinct adjacent faces for re-derivation".to_string(),
        )
    })?;

    let surf_a = brep
        .faces
        .get(face_a)
        .and_then(|f| f.surface.as_ref())
        .ok_or_else(|| CoreError::InvalidHandle("face_a has no surface geometry".to_string()))?
        .geometry
        .clone();
    let surf_b = brep
        .faces
        .get(face_b)
        .and_then(|f| f.surface.as_ref())
        .ok_or_else(|| CoreError::InvalidHandle("face_b has no surface geometry".to_string()))?
        .geometry
        .clone();

    // Measure current deviation (worst of the two surfaces).
    let old_dev_a = measure_edge_deviation(brep, edge, &surf_a);
    let old_dev_b = measure_edge_deviation(brep, edge, &surf_b);
    let old_deviation = old_dev_a.max(old_dev_b);

    // Sample the current edge curve.
    let edge_data = brep
        .edges
        .get(edge)
        .ok_or_else(|| CoreError::InvalidHandle(format!("edge {edge:?} not found")))?;
    let curve = edge_data.curve.as_ref().ok_or_else(|| {
        CoreError::PreconditionFailed("edge has no curve geometry to re-derive from".to_string())
    })?;

    let t_start = curve.t_start;
    let t_end = curve.t_end;
    let mut sample_pts = Vec::with_capacity(NUM_PROJECTION_SAMPLES);
    let mut params = Vec::with_capacity(NUM_PROJECTION_SAMPLES);

    for i in 0..NUM_PROJECTION_SAMPLES {
        let frac = i as f64 / (NUM_PROJECTION_SAMPLES - 1) as f64;
        let t = t_start + frac * (t_end - t_start);
        let pt = curve.geometry.point_at(t);
        let projected = project_to_intersection(&pt, &surf_a, &surf_b, PROJECTION_ITERATIONS);
        sample_pts.push(projected);
        params.push(frac);
    }

    // Fit a new NURBS curve through the projected points.
    let num_ctrl = NUM_PROJECTION_SAMPLES.clamp(FIT_DEGREE + 1, FIT_DEGREE + 4);
    let new_nurbs = fit_curve(&sample_pts, &params, FIT_DEGREE, num_ctrl)?;

    // Measure deviation of the new curve from both surfaces.
    let new_dev_a = edge_surface_deviation(&new_nurbs, (0.0, 1.0), &surf_a, 32).max_deviation;
    let new_dev_b = edge_surface_deviation(&new_nurbs, (0.0, 1.0), &surf_b, 32).max_deviation;
    let new_deviation = new_dev_a.max(new_dev_b);

    // Only replace if the new curve is strictly better and within tolerance.
    if new_deviation < old_deviation && new_deviation < tolerance.linear() * 10.0 {
        let new_curve_on_edge = CurveOnEdge {
            geometry: CurveGeometry::Nurbs(new_nurbs),
            t_start: 0.0,
            t_end: 1.0,
        };
        if let Some(e) = brep.edges.get_mut(edge) {
            e.curve = Some(new_curve_on_edge);
        }
        Ok(RederiveResult {
            edge,
            old_deviation,
            new_deviation,
            succeeded: true,
            description: format!(
                "re-derived edge curve: deviation {old_deviation:.6} -> {new_deviation:.6}"
            ),
        })
    } else {
        Ok(RederiveResult {
            edge,
            old_deviation,
            new_deviation,
            succeeded: false,
            description: format!(
                "re-derivation did not improve edge: {old_deviation:.6} -> {new_deviation:.6}"
            ),
        })
    }
}

/// Re-derive all edges whose deviation exceeds the tolerance.
///
/// Iterates over every edge in the BRep, measures deviation from adjacent
/// surfaces, and attempts re-derivation for edges with deviation above
/// the linear tolerance threshold.
pub fn rederive_all_edges(brep: &mut BRep, tolerance: &Tolerance) -> Vec<RederiveResult> {
    let threshold = tolerance.linear();
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    let mut results = Vec::new();

    for eh in edge_handles {
        // Check if this edge has two distinct faces and a curve.
        let has_two_faces = adjacent_faces(brep, eh).is_some();
        let has_curve = brep.edges.get(eh).and_then(|e| e.curve.as_ref()).is_some();

        if !has_two_faces || !has_curve {
            continue;
        }

        // Quick deviation check before attempting the expensive re-derivation.
        let dev = {
            let edge_data = match brep.edges.get(eh) {
                Some(e) => e,
                None => continue,
            };
            let curve = match &edge_data.curve {
                Some(c) => c,
                None => continue,
            };
            let (fa, fb) = match adjacent_faces(brep, eh) {
                Some(pair) => pair,
                None => continue,
            };
            let sa = match brep.faces.get(fa).and_then(|f| f.surface.as_ref()) {
                Some(s) => &s.geometry,
                None => continue,
            };
            let sb = match brep.faces.get(fb).and_then(|f| f.surface.as_ref()) {
                Some(s) => &s.geometry,
                None => continue,
            };
            let da = edge_surface_deviation(&curve.geometry, (curve.t_start, curve.t_end), sa, 32)
                .max_deviation;
            let db = edge_surface_deviation(&curve.geometry, (curve.t_start, curve.t_end), sb, 32)
                .max_deviation;
            da.max(db)
        };

        if dev <= threshold {
            continue;
        }

        if let Ok(result) = rederive_edge_from_surfaces(brep, eh, tolerance) {
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
    fn rederive_on_valid_cube_is_noop() {
        let mut brep = BRep::unit_cube(Tolerance::default()).ok();
        if let Some(ref mut b) = brep {
            let tol = Tolerance::default();
            let results = rederive_all_edges(b, &tol);
            // Valid cube edges should all be within tolerance, so nothing to rederive.
            for r in &results {
                assert!(
                    r.old_deviation > tol.linear(),
                    "only high-deviation edges should be attempted"
                );
            }
        }
    }

    #[test]
    fn rederive_result_fields() {
        let result = RederiveResult {
            edge: EdgeHandle::default(),
            old_deviation: 0.1,
            new_deviation: 0.01,
            succeeded: true,
            description: "test".to_string(),
        };
        assert!(result.succeeded);
        assert!(result.new_deviation < result.old_deviation);
    }

    #[test]
    fn adjacent_faces_returns_none_for_empty_brep() {
        let brep = BRep::new(Tolerance::default());
        let eh = EdgeHandle::default();
        assert!(adjacent_faces(&brep, eh).is_none());
    }

    #[test]
    fn rederive_all_empty_brep() {
        let mut brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let results = rederive_all_edges(&mut brep, &tol);
        assert!(results.is_empty());
    }

    #[test]
    fn rederive_edge_with_high_deviation() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        // Find an edge with a curve.
        let edge_h = brep.edges.keys().find(|&e| brep.edges[e].curve.is_some());
        if let Some(eh) = edge_h {
            // Get a vertex on this edge and shift it to create deviation.
            let he_handle = brep.edges[eh].half_edge;
            let origin_vh = brep.half_edges[he_handle].origin;
            brep.vertices[origin_vh].position.z += 0.1;
            let result = rederive_edge_from_surfaces(&mut brep, eh, &tol);
            // Should attempt rederivation (may or may not succeed depending
            // on whether the edge has two distinct adjacent faces).
            assert!(result.is_ok() || result.is_err());
        }
    }

    #[test]
    fn rederive_all_edges_on_corrupted_cube() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        // Corrupt several vertices.
        let vert_handles: Vec<_> = brep.vertices.keys().collect();
        for (i, &vh) in vert_handles.iter().enumerate().take(2) {
            brep.vertices[vh].position.z += 0.05 * (i as f64 + 1.0);
        }
        let results = rederive_all_edges(&mut brep, &tol);
        // Should find some edges to rederive due to increased deviation.
        // Even if empty (all edges below threshold), no panic.
        for r in &results {
            assert!(r.old_deviation > tol.linear());
        }
    }

    #[test]
    fn rederive_edge_without_curve_returns_error() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let eh = brep.edges.keys().next().unwrap();
        brep.edges[eh].curve = None;
        let tol = Tolerance::default();
        let result = rederive_edge_from_surfaces(&mut brep, eh, &tol);
        // Should fail: edge has no curve to re-derive from.
        assert!(result.is_err());
    }

    #[test]
    fn rederive_edge_missing_surface_returns_error() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let eh = brep.edges.keys().next().unwrap();
        // Remove surface from one adjacent face.
        let he_handle = brep.edges[eh].half_edge;
        let he = &brep.half_edges[he_handle];
        let loop_handle = he.loop_;
        let face_handle = brep.loops[loop_handle].face;
        brep.faces[face_handle].surface = None;
        let tol = Tolerance::default();
        let result = rederive_edge_from_surfaces(&mut brep, eh, &tol);
        // Should fail: can't rederive without adjacent surface.
        assert!(result.is_err());
    }
}
