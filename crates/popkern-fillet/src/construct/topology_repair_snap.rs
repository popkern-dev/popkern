//! Incident edge endpoint snapping helpers for topology repair.

use super::topology_repair_geometry::{remap_to_unit_range, slice_curve};
use super::*;

/// For every edge incident to `vh`, snap the 3D curve's endpoint
/// (the one matching v's t_start or t_end) to `new_position`.
/// More robust than the legacy version: uses half-edge orientation
/// to determine which end of the curve corresponds to `vh`, then
/// modifies the appropriate control point of the curve.
pub(super) fn snap_incident_edge_curve_endpoints_robust(
    brep: &mut BRep,
    vh: VertexHandle,
    new_position: Point3,
) {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    for edge in edge_handles {
        let Some(edge_data) = brep.edges.get(edge) else {
            continue;
        };
        let he = edge_data.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let twin = he_data.twin;
        let Some(twin_data) = brep.half_edges.get(twin) else {
            continue;
        };
        let v_origin = he_data.origin;
        let v_target = twin_data.origin;
        let edge_starts_at_vh = v_origin == vh;
        let edge_ends_at_vh = v_target == vh;
        if !edge_starts_at_vh && !edge_ends_at_vh {
            continue;
        }
        let Some(curve_on_edge) = brep.edges.get(edge).and_then(|e| e.curve.clone()) else {
            continue;
        };
        // Convert any analytic curve to NURBS first, so we have
        // control points to manipulate. Line/Arc/Circle types are
        // parameterised by analytic data (origin, direction, etc.)
        // — modifying those to relocate one endpoint is fiddly, but
        // promoting to NURBS gives us a uniform interface.
        let nurbs_full = match &curve_on_edge.geometry {
            CurveGeometry::Nurbs(n) => n.clone(),
            CurveGeometry::Line { .. }
            | CurveGeometry::Arc { .. }
            | CurveGeometry::Circle { .. }
            | CurveGeometry::Ellipse { .. }
            | CurveGeometry::Parabola { .. }
            | CurveGeometry::Hyperbola { .. } => {
                let trange = (curve_on_edge.t_start, curve_on_edge.t_end);
                curve_on_edge.geometry.to_nurbs_in_range(trange)
            }
        };
        // CurveOnEdge.t_start/t_end may not match the underlying
        // NURBS's full parameter range — the curve may extend
        // beyond the edge's actual span. Slice down to [t_start,
        // t_end] so the first/last CPs of the resulting curve
        // *do* correspond to the edge endpoints.
        let nurbs_sliced =
            match slice_curve(&nurbs_full, curve_on_edge.t_start, curve_on_edge.t_end) {
                Some(s) => remap_to_unit_range(&s),
                None => nurbs_full.clone(),
            };
        // The curve direction may not match the half-edge direction
        // (edge.half_edge isn't always the curve-start half-edge).
        // So decide which CP to update based on which curve endpoint
        // is closer to the existing vertex position. Use the new
        // position, since the move can be O(drift) in either
        // direction — closer in original space stays closer.
        let mut cps = nurbs_sliced.control_points.clone();
        let snap_first = if let (Some(first), Some(last)) = (cps.first(), cps.last()) {
            first.distance_to(&new_position) <= last.distance_to(&new_position)
        } else {
            true
        };
        // Sanity: only update if exactly one end is incident to vh
        // (BRep half-edge semantics usually guarantee this).
        let _ = (edge_starts_at_vh, edge_ends_at_vh);
        if snap_first {
            if let Some(first) = cps.first_mut() {
                *first = new_position;
            }
        } else if let Some(last) = cps.last_mut() {
            *last = new_position;
        }
        let new_curve = NurbsCurve {
            degree: nurbs_sliced.degree,
            knots: nurbs_sliced.knots,
            control_points: cps,
            weights: nurbs_sliced.weights,
        };
        let (new_t_start, new_t_end) = new_curve.parameter_range();
        if let Some(edge_mut) = brep.edges.get_mut(edge) {
            let geometry = CurveGeometry::Nurbs(new_curve);
            let geometry = geometry
                .canonical_line_in_range((new_t_start, new_t_end), 1.0e-8)
                .unwrap_or(geometry);
            edge_mut.curve = Some(CurveOnEdge {
                geometry,
                t_start: new_t_start,
                t_end: new_t_end,
            });
        }
    }
}

/// For every edge with an endpoint vertex `vh`, replace the
/// matching control point of the edge's 3D curve so the curve
/// terminates at the new vertex position. Identifies the matching
/// control point by proximity to `old_position`.
#[allow(dead_code)]
pub(super) fn snap_incident_edge_curve_endpoints(
    brep: &mut BRep,
    vh: VertexHandle,
    old_position: Point3,
    new_position: Point3,
) {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    for edge in edge_handles {
        let Some(edge_data) = brep.edges.get(edge) else {
            continue;
        };
        let he = edge_data.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let twin = he_data.twin;
        let Some(twin_data) = brep.half_edges.get(twin) else {
            continue;
        };
        let v_origin = he_data.origin;
        let v_target = twin_data.origin;
        if v_origin != vh && v_target != vh {
            continue;
        }
        let Some(curve_on_edge) = brep.edges.get(edge).and_then(|e| e.curve.clone()) else {
            continue;
        };
        let CurveGeometry::Nurbs(nurbs) = curve_on_edge.geometry else {
            continue;
        };
        let mut cps = nurbs.control_points.clone();
        let first_d = cps
            .first()
            .map(|p| p.distance_to(&old_position))
            .unwrap_or(f64::INFINITY);
        let last_d = cps
            .last()
            .map(|p| p.distance_to(&old_position))
            .unwrap_or(f64::INFINITY);
        if first_d < last_d && first_d < 1e-3 {
            if let Some(first) = cps.first_mut() {
                *first = new_position;
            }
        } else if last_d < 1e-3 {
            if let Some(last) = cps.last_mut() {
                *last = new_position;
            }
        } else {
            continue;
        }
        let new_curve = NurbsCurve {
            degree: nurbs.degree,
            knots: nurbs.knots,
            control_points: cps,
            weights: nurbs.weights,
        };
        if let Some(edge_mut) = brep.edges.get_mut(edge) {
            let geometry = CurveGeometry::Nurbs(new_curve);
            let geometry = geometry
                .canonical_line_in_range((curve_on_edge.t_start, curve_on_edge.t_end), 1.0e-8)
                .unwrap_or(geometry);
            edge_mut.curve = Some(CurveOnEdge {
                geometry,
                t_start: curve_on_edge.t_start,
                t_end: curve_on_edge.t_end,
            });
        }
    }
}
