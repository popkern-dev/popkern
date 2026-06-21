//! Surface-derived split-edge vertex placement.
//!
//! Stage 1 of `docs/decisions/2026-06-21-split-edge-surface-derived-
//! vertex-placement.md`. When `BRep::split_edge` encounters an edge
//! with no curve representation, this module's helper derives the
//! new vertex 3D position from an incident face's surface via the
//! half-edge's pcurve. Falls back to `None` (caller chord-fallback)
//! when no surface or pcurve is available for the half-edge.
//!
//! The helper is `pub(super)` and called from `split.rs` in Stage 2.
//! Stage 1 ships the helper + unit tests; behavior is unchanged
//! until Stage 2 wires it in.

use crate::geometry::Point3;
use crate::geometry::surface::SurfaceEval;
use crate::topology::{BRep, HalfEdgeHandle};

/// Compute the 3D position for a split-edge vertex by evaluating the
/// face's surface at the pcurve's image of `t`.
///
/// `t` is the split parameter in `[0, 1]`, in the half-edge's
/// direction. Returns `None` when:
///
/// - the half-edge handle is stale,
/// - the half-edge has no pcurve,
/// - the half-edge's loop's face has no surface, or
/// - the pcurve's parameter range is degenerate.
///
/// Returns `Some(Point3)` otherwise — the 3D point lies on the
/// face's surface at the pcurve's image of `t`, which is by
/// construction the face's surface boundary at that parameter.
pub(super) fn surface_derived_split_position(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    t: f64,
) -> Option<Point3> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let pcurve = half_edge_data.pcurve.as_ref()?;
    let loop_data = brep.loops.get(half_edge_data.loop_)?;
    let face = brep.faces.get(loop_data.face)?;
    let surface_on_face = face.surface.as_ref()?;
    let (pc_start, pc_end) = pcurve.parameter_range();
    let pc_range = pc_end - pc_start;
    if pc_range.abs() < 1.0e-30 {
        return None;
    }
    let pc_t = pc_start + t * pc_range;
    let uv = pcurve.point_at(pc_t);
    Some(surface_on_face.geometry.point_at(uv.x, uv.y))
}

#[cfg(test)]
#[path = "split_surface_derived_tests.rs"]
mod tests;
