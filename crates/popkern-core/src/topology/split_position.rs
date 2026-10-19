//! Split position + curve propagation for [`BRep::split_edge`].
//!
//! Factored out of [`super::split`] so the rule "split.rs ≤ 500 lines"
//! holds after the β-1-A propagation and surface-derived fallback
//! were added (ADR `2026-06-21-split-edge-surface-derived-vertex-
//! placement.md`).

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::nurbs_curve::NurbsCurve;
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
///   3. Chord interpolation between `pos_a` and `pos_b` (the original
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
            // Direction decision shared with the curve-routing site in
            // `split_edge` (they must agree). On CLOSED edges the
            // endpoint-distance witness is vacuous — the tie used to
            // break arbitrarily and could attach each sub-curve to the
            // WRONG half, flipping the loop's geometric winding while
            // the topology stayed consistent (measured: the cylinder
            // rim split in the PK-325 arc-join path served a CW outer
            // cap); ambiguous cases use a pcurve-image witness.
            let curve_aligned = curve_alignment_for_split(brep, he_ab, he_ba, curve, pos_a, pos_b);
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

/// Shared curve-direction decision for both `split_edge` sites (vertex
/// placement and NURBS half routing must agree).
pub(super) fn curve_alignment_for_split(
    brep: &BRep,
    he_ab: HalfEdgeHandle,
    he_ba: HalfEdgeHandle,
    curve: &CurveOnEdge,
    pos_a: Point3,
    pos_b: Point3,
) -> bool {
    let p_t_start = curve.geometry.point_at(curve.t_start);
    let p_t_end = curve.geometry.point_at(curve.t_end);
    let d_start = pos_a.distance_to(&p_t_start);
    let d_end = pos_a.distance_to(&p_t_end);
    let curve_scale = p_t_start.distance_to(&pos_b).max(d_start).max(d_end);
    let ambiguous = (d_start - d_end).abs() <= 1.0e-9 * (1.0 + curve_scale)
        || p_t_start.distance_to(&p_t_end) <= 1.0e-9 * (1.0 + curve_scale);
    if ambiguous {
        closed_edge_alignment_witness(brep, he_ab, curve)
            .or_else(|| closed_edge_alignment_witness(brep, he_ba, curve).map(|w| !w))
            .unwrap_or(true)
    } else {
        d_start <= d_end
    }
}

/// Expected 3D midpoints of the A→V and V→B spans, measured from the
/// half-edges' ORIGINAL (pre-split) pcurves mapped through their face
/// surfaces. `None` when no pcurve/surface pair resolves.
pub(super) fn closed_split_expected_midpoints(
    brep: &BRep,
    he_ab: HalfEdgeHandle,
    he_ba: HalfEdgeHandle,
    pcurve_ab: &Option<crate::geometry::pcurve::Pcurve>,
    pcurve_ba: &Option<crate::geometry::pcurve::Pcurve>,
    t: f64,
) -> Option<(Point3, Point3)> {
    // Fractions along he_ab of the two span midpoints.
    if let Some(mids) = pcurve_span_midpoints(brep, he_ab, pcurve_ab, 0.5 * t, 0.5 * (1.0 + t)) {
        return Some(mids);
    }
    // he_ba runs B→A: the A→V span sits at fraction 1 − t/2, the V→B
    // span at (1 − t)/2.
    pcurve_span_midpoints(brep, he_ba, pcurve_ba, 1.0 - 0.5 * t, 0.5 * (1.0 - t))
}

fn pcurve_span_midpoints(
    brep: &BRep,
    he: HalfEdgeHandle,
    pcurve: &Option<crate::geometry::pcurve::Pcurve>,
    frac_a: f64,
    frac_b: f64,
) -> Option<(Point3, Point3)> {
    use crate::geometry::surface::SurfaceEval;
    let pcurve = pcurve.as_ref()?;
    let he_data = brep.half_edges.get(he)?;
    let face = brep.loops.get(he_data.loop_)?.face;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let eval = |frac: f64| -> Point3 {
        let pt = pcurve.t_start + frac * (pcurve.t_end - pcurve.t_start);
        let uv = pcurve.point_at(pt);
        surface.point_at(uv.x, uv.y)
    };
    Some((eval(frac_a), eval(frac_b)))
}

/// Direction witness for closed/degenerate edges: does `he`'s pcurve
/// traverse the edge curve forward? Samples both slightly inside their
/// ranges and compares surface-mapped images — the offset common to
/// both candidates cancels, so the witness stays valid even when the
/// face surface has been swapped for an offset of itself (the shell
/// pipeline's pre-recompute state).
fn closed_edge_alignment_witness(
    brep: &BRep,
    he: HalfEdgeHandle,
    curve: &CurveOnEdge,
) -> Option<bool> {
    use crate::geometry::surface::SurfaceEval;
    let he_data = brep.half_edges.get(he)?;
    let pcurve = he_data.pcurve.as_ref()?;
    let face = brep.loops.get(he_data.loop_)?.face;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    const EPS: f64 = 0.05;
    let pt = pcurve.t_start + EPS * (pcurve.t_end - pcurve.t_start);
    let uv = pcurve.point_at(pt);
    let along_he = surface.point_at(uv.x, uv.y);
    let near_start = curve
        .geometry
        .point_at(curve.t_start + EPS * (curve.t_end - curve.t_start));
    let near_end = curve
        .geometry
        .point_at(curve.t_end - EPS * (curve.t_end - curve.t_start));
    let d_s = along_he.distance_to(&near_start);
    let d_e = along_he.distance_to(&near_end);
    if (d_s - d_e).abs() <= 1.0e-9 * (1.0 + d_s.max(d_e)) {
        return None; // still ambiguous — let the caller try the twin
    }
    Some(d_s <= d_e)
}

// Split curve geometry if present. Use the same direction-
// alignment logic as `split_pos` above so the curve is
// sampled consistently with the half-edge direction.
/// Route the split NURBS halves onto the two edges, agreeing with the
/// vertex-placement direction decision (shared witness) and using
/// direction-agnostic endpoint costs plus pcurve-image midpoint terms
/// so CLOSED edges (degenerate endpoint costs) route correctly.
#[allow(clippy::too_many_arguments)]
pub(super) fn route_split_curves(
    brep: &mut BRep,
    curve_on_edge: &Option<CurveOnEdge>,
    he_ab: HalfEdgeHandle,
    he_ba: HalfEdgeHandle,
    pcurve_ab: &Option<crate::geometry::pcurve::Pcurve>,
    pcurve_ba: &Option<crate::geometry::pcurve::Pcurve>,
    edge: crate::topology::EdgeHandle,
    edge_b: crate::topology::EdgeHandle,
    pos_a: Point3,
    pos_b: Point3,
    split_pos: Point3,
    t: f64,
) {
    if let Some(curve) = curve_on_edge {
        // Circles and ellipses are parameterised by angle without bound:
        // their two sub-ranges (set by `split_position_and_propagation`
        // with the same direction decision) are already exact. Their
        // rational form is one turn on [0, 2π] whatever the trim, with
        // parameter ≠ angle inside each quarter, so a NURBS split left an
        // edge trimmed past one turn with the whole turn (measured: a
        // closed ellipse rim on jack_in_the_box trimmed [2.98, 6.34]
        // split at 2π kept [0, 2π], a 3.44 shell gap).
        if matches!(
            curve.geometry,
            crate::geometry::curve::CurveGeometry::Circle { .. }
                | crate::geometry::curve::CurveGeometry::Ellipse { .. }
        ) {
            return;
        }
        let curve_aligned = curve_alignment_for_split(brep, he_ab, he_ba, curve, pos_a, pos_b);
        let aligned_t = if curve_aligned { t } else { 1.0 - t };
        let t_param = curve.t_start + aligned_t * (curve.t_end - curve.t_start);

        // Try to split NURBS; if the curve is analytic, convert to NURBS first.
        let nurbs = curve
            .geometry
            .to_nurbs_in_range((curve.t_start, curve.t_end));
        // An arc's rational form agrees with its angle only at quarter
        // points: split it where its POINT is the split vertex. Read
        // at the angle, a tee's branch-cap rim split where a fillet
        // strip ran out on the cap left its cylinder by 5.2e-2 (PK-766).
        let t_param = match curve.geometry {
            crate::geometry::curve::CurveGeometry::Arc { .. } => nurbs
                .closest_point(&curve.geometry.point_at(t_param))
                .unwrap_or(t_param),
            _ => t_param,
        };
        match nurbs.split_at(t_param) {
            Ok((left, right)) => {
                use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};

                // The original edge spans vertex_a → vertex_b in
                // half-edge convention, but the curve geometry
                // may have been authored in the *opposite*
                // direction (e.g. by a boolean engine that
                // captured intersection curves without enforcing
                // a canonical orientation). After splitting we
                // therefore inspect each half's endpoints and
                // route them to the half-edge whose vertex
                // positions they actually match. Without this
                // routing the loop validator sees a curve whose
                // span is the *other* half of the original edge,
                // producing an OpenLoop with a gap equal to the
                // length of the misrouted segment.
                let left_endpoints = {
                    let (t0, t1) = left.parameter_range();
                    (left.point_at(t0), left.point_at(t1))
                };
                let right_endpoints = {
                    let (t0, t1) = right.parameter_range();
                    (right.point_at(t0), right.point_at(t1))
                };
                let dist2 = |a: &Point3, b: &Point3| -> f64 {
                    (a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)
                };
                let half_mid = |half: &NurbsCurve| -> Point3 {
                    let (t0, t1) = half.parameter_range();
                    half.point_at(0.5 * (t0 + t1))
                };
                // Score "LEFT to old, RIGHT to new" vs "RIGHT to old, LEFT to new".
                // Each option pays the squared distance from the
                // expected endpoint to the actual vertex positions.
                // On CLOSED edges endpoint terms are DEGENERATE
                // (both halves connect the same two positions), so
                // the score also pays midpoint terms against the
                // half-edge's own traversal (original pcurve image
                // on the face surface): the half that covers the
                // A→V span of the traversal must land on the
                // original edge, else each sub-edge carries the
                // OTHER half's geometry and the loop's measured
                // winding flips while topology validates clean
                // (found by the PK-325 arc-join rim splits).
                let expected_mids =
                    closed_split_expected_midpoints(brep, he_ab, he_ba, pcurve_ab, pcurve_ba, t);
                let (mid_natural, mid_swapped) = match &expected_mids {
                    Some((exp_a, exp_b)) => (
                        dist2(&half_mid(&left), exp_a) + dist2(&half_mid(&right), exp_b),
                        dist2(&half_mid(&right), exp_a) + dist2(&half_mid(&left), exp_b),
                    ),
                    None => (0.0, 0.0),
                };
                // Edge curves are stored direction-free (the half-
                // edge convention resolves direction), so endpoint
                // matching must not privilege one traversal: pair
                // each half's endpoints with the edge's vertices
                // whichever way round is closer. A directional
                // pairing charged the correct closed-edge routing a
                // full diameter and always outvoted the midpoint
                // witness.
                let endpoint_cost = |ends: &(Point3, Point3), va: &Point3, vb: &Point3| -> f64 {
                    let fwd = dist2(&ends.0, va) + dist2(&ends.1, vb);
                    let rev = dist2(&ends.0, vb) + dist2(&ends.1, va);
                    fwd.min(rev)
                };
                let cost_natural = endpoint_cost(&left_endpoints, &pos_a, &split_pos)
                    + endpoint_cost(&right_endpoints, &split_pos, &pos_b)
                    + mid_natural;
                let cost_swapped = endpoint_cost(&right_endpoints, &pos_a, &split_pos)
                    + endpoint_cost(&left_endpoints, &split_pos, &pos_b)
                    + mid_swapped;
                let (old_half, new_half) = if cost_natural <= cost_swapped {
                    (left, right)
                } else {
                    (right, left)
                };

                let old_range = old_half.parameter_range();
                let old_geometry = CurveGeometry::Nurbs(old_half);
                let old_geometry = old_geometry
                    .canonical_line_in_range(old_range, brep.tolerance.linear().max(1.0e-9))
                    .unwrap_or(old_geometry);
                brep.edges[edge].curve = Some(CurveOnEdge {
                    geometry: old_geometry,
                    t_start: old_range.0,
                    t_end: old_range.1,
                });

                let new_range = new_half.parameter_range();
                let new_geometry = CurveGeometry::Nurbs(new_half);
                let new_geometry = new_geometry
                    .canonical_line_in_range(new_range, brep.tolerance.linear().max(1.0e-9))
                    .unwrap_or(new_geometry);
                brep.edges[edge_b].curve = Some(CurveOnEdge {
                    geometry: new_geometry,
                    t_start: new_range.0,
                    t_end: new_range.1,
                });
            }
            Err(_) => {
                // If split fails (degenerate), keep original on edge_a,
                // leave edge_b without a curve.
            }
        }
    }
}
