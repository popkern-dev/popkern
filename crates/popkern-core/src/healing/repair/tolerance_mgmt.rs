//! Tolerance hierarchy management.
//!
//! B-Rep tolerances follow a hierarchy: vertex tolerances must be at
//! least as large as the tolerances of their incident edges, and edge
//! tolerances must be at least as large as the tolerances of their
//! parent faces.  This module provides functions to suggest, enforce,
//! and auto-adjust tolerances across the model.

use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::healing::diagnostics::EntityRef;
use crate::healing::options::ToleranceBounds;
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Record of a single tolerance adjustment.
#[derive(Debug, Clone)]
pub struct ToleranceAdjustment {
    /// Which entity was adjusted.
    pub entity: EntityRef,
    /// Previous tolerance value (`None` if the entity used the BRep default).
    pub old_tolerance: Option<f64>,
    /// New tolerance value after adjustment.
    pub new_tolerance: f64,
}

/// Suggest a tolerance for an edge based on its curve-to-surface deviation.
///
/// Samples 8 points along the edge curve and measures the distance from
/// each to the adjacent face surfaces.  Returns the maximum deviation,
/// which is a reasonable tolerance for the edge.
///
/// Returns `None` if the edge lacks curve geometry or adjacent faces
/// lack surface geometry.
pub fn suggest_tolerance(brep: &BRep, edge: EdgeHandle) -> Option<f64> {
    let edge_data = brep.edges.get(edge)?;
    let curve_on_edge = edge_data.curve.as_ref()?;

    // Find the two faces adjacent to this edge via its half-edges.
    let he_handle = edge_data.half_edge;
    let he = brep.half_edges.get(he_handle)?;
    let twin = brep.half_edges.get(he.twin)?;

    let face_a_handle = brep.loops.get(he.loop_)?.face;
    let face_b_handle = brep.loops.get(twin.loop_)?.face;

    let surface_a = brep.faces.get(face_a_handle)?.surface.as_ref()?;
    let surface_b = brep.faces.get(face_b_handle)?.surface.as_ref()?;

    let n_samples = 8;
    let t0 = curve_on_edge.t_start;
    let t1 = curve_on_edge.t_end;
    let mut max_dev = 0.0_f64;

    for i in 0..n_samples {
        let t = t0 + (t1 - t0) * (i as f64) / ((n_samples - 1).max(1) as f64);
        let pt = curve_on_edge.geometry.point_at(t);

        // Distance to surface A.
        if let Some((u, v)) = surface_a.geometry.closest_point(&pt) {
            let on_surface = surface_a.geometry.point_at(u, v);
            max_dev = max_dev.max(pt.distance_to(&on_surface));
        }

        // Distance to surface B.
        if let Some((u, v)) = surface_b.geometry.closest_point(&pt) {
            let on_surface = surface_b.geometry.point_at(u, v);
            max_dev = max_dev.max(pt.distance_to(&on_surface));
        }
    }

    if max_dev > 0.0 { Some(max_dev) } else { None }
}

/// Enforce the tolerance hierarchy: vertex >= edge >= face.
///
/// Walks every edge and its endpoints.  If a vertex tolerance is
/// smaller than its incident edge tolerance, the vertex tolerance is
/// raised.  Similarly, if an edge tolerance is smaller than its parent
/// face tolerance, the edge tolerance is raised.
///
/// Returns a list of all adjustments made.
pub fn enforce_tolerance_hierarchy(brep: &mut BRep) -> Vec<ToleranceAdjustment> {
    let mut adjustments = Vec::new();

    // Phase 1: edge tolerance >= face tolerance.
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    for fh in face_handles {
        let face_tol = brep.face_tolerance(fh);
        let Some(face) = brep.faces.get(fh) else {
            continue;
        };
        let outer_loop = face.outer_loop;
        let Some(loop_) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start_he = loop_.half_edge;
        if !brep.half_edges.contains_key(start_he) {
            continue;
        }

        let mut current = start_he;
        loop {
            let Some(he) = brep.half_edges.get(current) else {
                break;
            };
            let edge_handle = he.edge;
            let next = he.next;

            let edge_tol = brep.edge_tolerance(edge_handle);
            if edge_tol < face_tol {
                let old = brep.edges.get(edge_handle).and_then(|e| e.tolerance);
                if let Some(e) = brep.edges.get_mut(edge_handle) {
                    e.tolerance = Some(face_tol);
                }
                adjustments.push(ToleranceAdjustment {
                    entity: EntityRef::Edge(edge_handle),
                    old_tolerance: old,
                    new_tolerance: face_tol,
                });
            }

            current = next;
            if current == start_he {
                break;
            }
        }
    }

    // Phase 2: vertex tolerance >= edge tolerance.
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    for eh in edge_handles {
        let edge_tol = brep.edge_tolerance(eh);
        let Some(edge_data) = brep.edges.get(eh) else {
            continue;
        };
        let he_handle = edge_data.half_edge;

        // Collect vertex handles from both half-edges.
        let mut vertex_handles: Vec<VertexHandle> = Vec::new();
        if let Some(he) = brep.half_edges.get(he_handle) {
            vertex_handles.push(he.origin);
            if let Some(twin) = brep.half_edges.get(he.twin) {
                vertex_handles.push(twin.origin);
            }
        }

        for vh in vertex_handles {
            let vtx_tol = brep.vertex_tolerance(vh);
            if vtx_tol < edge_tol {
                let old = brep.vertices.get(vh).and_then(|v| v.tolerance);
                if let Some(v) = brep.vertices.get_mut(vh) {
                    v.tolerance = Some(edge_tol);
                }
                adjustments.push(ToleranceAdjustment {
                    entity: EntityRef::Vertex(vh),
                    old_tolerance: old,
                    new_tolerance: edge_tol,
                });
            }
        }
    }

    adjustments
}

/// Clamp all entity tolerances to the given bounds, then enforce hierarchy.
///
/// First passes over vertices, edges, and faces, clamping any explicit
/// tolerance to `[bounds.min_linear, bounds.max_linear]`.  Then calls
/// [`enforce_tolerance_hierarchy`] to restore monotonicity.
///
/// Returns all adjustments made (clamping + hierarchy enforcement).
pub fn auto_adjust_tolerances(
    brep: &mut BRep,
    bounds: &ToleranceBounds,
) -> Vec<ToleranceAdjustment> {
    let mut adjustments = Vec::new();

    // Clamp vertex tolerances.
    let vhs: Vec<VertexHandle> = brep.vertices.keys().collect();
    for vh in vhs {
        if let Some(v) = brep.vertices.get(vh)
            && let Some(tol) = v.tolerance
        {
            let clamped = tol.clamp(bounds.min_linear, bounds.max_linear);
            if (clamped - tol).abs() > f64::EPSILON {
                if let Some(v) = brep.vertices.get_mut(vh) {
                    v.tolerance = Some(clamped);
                }
                adjustments.push(ToleranceAdjustment {
                    entity: EntityRef::Vertex(vh),
                    old_tolerance: Some(tol),
                    new_tolerance: clamped,
                });
            }
        }
    }

    // Clamp edge tolerances.
    let ehs: Vec<EdgeHandle> = brep.edges.keys().collect();
    for eh in ehs {
        if let Some(e) = brep.edges.get(eh)
            && let Some(tol) = e.tolerance
        {
            let clamped = tol.clamp(bounds.min_linear, bounds.max_linear);
            if (clamped - tol).abs() > f64::EPSILON {
                if let Some(e) = brep.edges.get_mut(eh) {
                    e.tolerance = Some(clamped);
                }
                adjustments.push(ToleranceAdjustment {
                    entity: EntityRef::Edge(eh),
                    old_tolerance: Some(tol),
                    new_tolerance: clamped,
                });
            }
        }
    }

    // Clamp face tolerances.
    let fhs: Vec<FaceHandle> = brep.faces.keys().collect();
    for fh in fhs {
        if let Some(f) = brep.faces.get(fh)
            && let Some(tol) = f.tolerance
        {
            let clamped = tol.clamp(bounds.min_linear, bounds.max_linear);
            if (clamped - tol).abs() > f64::EPSILON {
                if let Some(f) = brep.faces.get_mut(fh) {
                    f.tolerance = Some(clamped);
                }
                adjustments.push(ToleranceAdjustment {
                    entity: EntityRef::Face(fh),
                    old_tolerance: Some(tol),
                    new_tolerance: clamped,
                });
            }
        }
    }

    // Enforce hierarchy after clamping.
    let hierarchy_adj = enforce_tolerance_hierarchy(brep);
    adjustments.extend(hierarchy_adj);

    adjustments
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn enforce_hierarchy_on_valid_cube_is_noop() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let adjustments = enforce_tolerance_hierarchy(&mut brep);
        // A fresh cube has no per-entity tolerances, so the BRep default
        // applies uniformly -- no hierarchy violations.
        assert!(
            adjustments.is_empty(),
            "expected no adjustments on a fresh cube, got {}",
            adjustments.len()
        );
    }

    #[test]
    fn enforce_hierarchy_raises_vertex_tolerance() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // Set a high edge tolerance on the first edge.
        let edge = brep.edges.keys().next().unwrap();
        brep.edges[edge].tolerance = Some(0.1);
        // Vertices use the BRep default (1e-6), which is < 0.1.
        let adjustments = enforce_tolerance_hierarchy(&mut brep);
        assert!(
            !adjustments.is_empty(),
            "expected vertex tolerance to be raised"
        );
        // Check that affected vertices now have tolerance >= 0.1.
        let he = brep.edges[edge].half_edge;
        let origin = brep.half_edges[he].origin;
        assert!(
            brep.vertex_tolerance(origin) >= 0.1,
            "vertex tolerance should be >= edge tolerance"
        );
    }

    #[test]
    fn auto_adjust_clamps_and_enforces() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // Set an out-of-bounds edge tolerance.
        let edge = brep.edges.keys().next().unwrap();
        brep.edges[edge].tolerance = Some(1.0); // way above default max 1e-3

        let bounds = ToleranceBounds::default();
        let adjustments = auto_adjust_tolerances(&mut brep, &bounds);
        assert!(!adjustments.is_empty(), "expected clamping adjustments");
        // Edge tolerance should now be clamped to max_linear.
        let edge_tol = brep.edges[edge].tolerance.unwrap();
        assert!(
            edge_tol <= bounds.max_linear + f64::EPSILON,
            "edge tolerance {edge_tol} should be <= {}",
            bounds.max_linear
        );
    }

    #[test]
    fn suggest_tolerance_on_cube_edge() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        // The cube's edge curves lie exactly on the adjacent plane surfaces,
        // so the suggested tolerance should be very small (or None if no curve).
        if let Some(suggested) = suggest_tolerance(&brep, edge) {
            assert!(
                suggested < 1e-3,
                "suggested tolerance {suggested} should be small for a perfect cube"
            );
        }
    }
}
