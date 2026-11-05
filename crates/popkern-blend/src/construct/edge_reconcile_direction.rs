//! Whether a half-edge runs its edge's curve backward over a range.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, HalfEdgeHandle};

/// Whether `half_edge` (on an open edge) starts at its edge curve's
/// `range.1` end: its origin vertex is clearly nearer that end.
pub(super) fn half_edge_runs_backward(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    range: (f64, f64),
) -> bool {
    let Some(h) = brep.half_edges.get(half_edge) else {
        return false;
    };
    let (Some(origin), Some(curve)) = (
        brep.vertices.get(h.origin).map(|v| v.position),
        brep.edges.get(h.edge).and_then(|e| e.curve.as_ref()),
    ) else {
        return false;
    };
    let (a, b) = (
        curve.geometry.point_at(range.0),
        curve.geometry.point_at(range.1),
    );
    let chord = a.distance_to(&b);
    let (da, db) = (origin.distance_to(&a), origin.distance_to(&b));
    chord > brep.tolerance.linear() && db + 1.0e-3 * chord < da
}
