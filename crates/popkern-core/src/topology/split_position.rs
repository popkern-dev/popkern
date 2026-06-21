//! Split position + curve propagation for [`BRep::split_edge`].
//!
//! Factored out of [`super::split`] so the rule "split.rs ≤ 500 lines"
//! holds after the β-1-A propagation and surface-derived fallback
//! were added (ADR `2026-06-21-split-edge-surface-derived-vertex-
//! placement.md`).

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::topology::{BRep, HalfEdgeHandle};

use super::split_surface_derived::surface_derived_split_position;

/// Compute the new vertex 3D position and the two sub-curves the
/// split should leave on the original and new edges.
///
/// Returns `(split_pos, edge_a_curve, edge_b_curve)`:
///   - `split_pos` is the 3D point for the new vertex.
///   - `edge_a_curve` is the replacement for the original edge's
///     curve (covers A→V of the parent curve).
///   - `edge_b_curve` is the curve assigned to the new V→B edge.
///
/// Preference order when the parent edge has no curve:
///   1. Surface-derived from `he_ab`'s face (pcurve image evaluated
///      on the face surface) — see
///      [`surface_derived_split_position`].
///   2. Surface-derived from `he_ba`'s face (twin).
///   3. Chord interpolation between `pos_a` and `pos_b` (legacy
///      fallback).
pub(super) fn split_position_and_propagation(
    brep: &BRep,
    curve_on_edge: &Option<CurveOnEdge>,
    he_ab: HalfEdgeHandle,
    he_ba: HalfEdgeHandle,
    pos_a: Point3,
    pos_b: Point3,
    t: f64,
) -> (Point3, Option<CurveOnEdge>, Option<CurveOnEdge>) {
    match curve_on_edge {
        Some(curve) => {
            let p_t_start = curve.geometry.point_at(curve.t_start);
            let p_t_end = curve.geometry.point_at(curve.t_end);
            let curve_aligned = pos_a.distance_to(&p_t_start) <= pos_a.distance_to(&p_t_end);
            let aligned_t = if curve_aligned { t } else { 1.0 - t };
            let t_param = curve.t_start + aligned_t * (curve.t_end - curve.t_start);
            let pos = curve.geometry.point_at(t_param);
            let (curve_at_a, curve_at_b) = if curve_aligned {
                (curve.t_start, curve.t_end)
            } else {
                (curve.t_end, curve.t_start)
            };
            let edge_a_curve = CurveOnEdge {
                geometry: curve.geometry.clone(),
                t_start: curve_at_a,
                t_end: t_param,
            };
            let edge_b_curve = CurveOnEdge {
                geometry: curve.geometry.clone(),
                t_start: t_param,
                t_end: curve_at_b,
            };
            (pos, Some(edge_a_curve), Some(edge_b_curve))
        }
        None => {
            let pos = surface_derived_split_position(brep, he_ab, t)
                .or_else(|| surface_derived_split_position(brep, he_ba, t))
                .unwrap_or_else(|| {
                    Point3::new(
                        pos_a.x + t * (pos_b.x - pos_a.x),
                        pos_a.y + t * (pos_b.y - pos_a.y),
                        pos_a.z + t * (pos_b.z - pos_a.z),
                    )
                });
            (pos, None, None)
        }
    }
}
