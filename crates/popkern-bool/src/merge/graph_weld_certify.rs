//! Weld certification helpers (split from `graph_weld.rs` under the
//! 500-line file budget).

/// FAR-TIER weld certificate: does `point` lie on `surface` within
/// `gate`? A cold NURBS projection can succeed into the wrong basin
/// of a multi-turn chart (the PK-463 class), so an off-carrier first
/// answer retries through the knot-span-seeded multi-start.
pub(super) fn on_carrier(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    point: &popkern_core::geometry::Point3,
    gate: f64,
) -> bool {
    use popkern_core::geometry::surface::{SurfaceEval as _, SurfaceGeometry};
    if let Some((u, v)) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, point)
        && surface.point_at(u, v).distance_to(point) <= gate
    {
        return true;
    }
    if let SurfaceGeometry::Nurbs(nurbs) = surface
        && let Some(seeds) = nurbs.projection_seeds()
        && let Some((u, v)) = popkern_core::geometry::surface::closest_point_from_seeds(
            nurbs,
            point,
            &seeds,
            (false, false),
            6,
        )
        && surface.point_at(u, v).distance_to(point) <= gate
    {
        return true;
    }
    false
}

use super::super::graph::IntersectionGraph;

// Beyond the tuned unconditional band, a VERIFIED annulus (up to
// 10× further): weld only when the candidate junction lies on the
// dead-end record curve's extension — perpendicular distance from
// the node to curve samples just past the endpoint. A shallow
// record × boundary crossing AMPLIFIES pcurve slop into the
// along-curve endpoint miss (measured, PK-464: the mouth-chamfer
// circle crossed the nut's crest-band helix boundary at slope
// 0.16, a sub-1e-3 boundary-pcurve wobble moved the clip window
// end 6.2e-4 from the true junction node, the ring segment
// dead-ended and the graze cascade unravelled the whole chamfer
// ring). The perpendicular gate keeps wrong-junction welds out:
// an off-curve node at the same spatial distance stays unwelded.
pub(super) fn on_extension(
    graph: &IntersectionGraph,
    incident: &[Vec<usize>],
    band: f64,
    n: usize,
    m_pos: &popkern_core::geometry::Point3,
) -> bool {
    {
        use popkern_core::geometry::curve::CurveEval as _;
        // Tangent-ray tolerance: the fit-slop class displaces the
        // endpoint mostly ALONG the curve; the perpendicular residual
        // of the true junction from the endpoint tangent stays within
        // curvature sag ((3.7e-2)²/2r ≈ 2.6e-4 on the nut's crest
        // circle) plus lateral slop.
        let perp_band = band;
        for &h in &incident[n] {
            let he = &graph.half_edges[h];
            let Some(c) = he
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[he.twin].curve.as_ref())
            else {
                // No stored curve on either side: fall back to the
                // CHORD tangent (twin origin -> this node). At weld
                // scales the chord direction is the outward tangent
                // to curvature sag (measured: the crest-arc piece
                // arrived curveless after its splits and the weld
                // silently skipped its only extension test).
                let p_n = graph.nodes[n].position;
                let q = graph.nodes[graph.half_edges[he.twin].origin].position;
                let Some(outward) = (p_n - q).normalized() else {
                    continue;
                };
                let to_m = *m_pos - p_n;
                let along = to_m.dot(&outward);
                let perp = (to_m - outward * along).length();
                if along >= -perp_band && perp < perp_band {
                    return true;
                }
                continue;
            };
            let span = c.t_end - c.t_start;
            if !span.is_finite() || span.abs() < 1e-30 {
                continue;
            }
            // Which curve end sits at `n`; outward tangent by an
            // inside finite difference (a clamped NURBS does not
            // extrapolate past its knots, so sampling beyond the end
            // is not an option).
            let p_n = graph.nodes[n].position;
            let d_start = c.geometry.point_at(c.t_start).distance_to(&p_n);
            let d_end = c.geometry.point_at(c.t_end).distance_to(&p_n);
            let (t_at, t_in) = if d_start <= d_end {
                (c.t_start, c.t_start + span * 1e-3)
            } else {
                (c.t_end, c.t_end - span * 1e-3)
            };
            let p_end = c.geometry.point_at(t_at);
            let inward = c.geometry.point_at(t_in) - p_end;
            let Some(outward) = (inward * -1.0).normalized() else {
                continue;
            };
            let to_m = *m_pos - p_end;
            let along = to_m.dot(&outward);
            let perp = (to_m - outward * along).length();
            if along < -perp_band {
                continue; // candidate sits back along the curve, not past the end
            }
            if perp < perp_band {
                return true;
            }
        }
        false
    }
}

/// Post-weld self-loop stub retirement (see the call site).
pub(super) fn retire_self_loop_stubs(graph: &mut IntersectionGraph, max_weld_d: f64, band: f64) {
    // POINT-LOOP STUB CLEANUP: welding two record endpoints a sliver
    // apart can close a micro record piece onto ITSELF (both ends the
    // welded node). Left registered, its pair sits in the wedge fan
    // as a zero-angle out-and-back and the rotational walk U-TURNS
    // through it, fusing the two regions meeting at the junction into
    // one fragment with one verdict (measured, PK-465: a 7.0e-5
    // mouth-junction stub on the angled tap's coil flank fused the
    // bore-void strip with the in-material flank; the served body
    // kept the strip as a masked ghost weave). Unregister such stubs
    // from faces and nodes; the gate is the weld scale itself.
    if max_weld_d > 0.0 {
        let gate = max_weld_d * 2.0 + band;
        let n_hes = graph.half_edges.len();
        for i in (0..n_hes).step_by(2) {
            let twin = graph.half_edges[i].twin;
            if !matches!(
                graph.half_edges[i].source,
                crate::merge::graph::EdgeSource::Shared { .. }
            ) {
                continue;
            }
            if graph.half_edges[i].origin != graph.half_edges[twin].origin {
                continue;
            }
            let node = graph.half_edges[i].origin;
            let anchor = graph.nodes[node].position;
            let extent = {
                use popkern_core::geometry::curve::CurveEval as _;
                match graph.half_edges[i]
                    .curve
                    .as_ref()
                    .or_else(|| graph.half_edges[twin].curve.as_ref())
                {
                    Some(c) => (0..=6)
                        .map(|k| {
                            let t = c.t_start + (c.t_end - c.t_start) * k as f64 / 6.0;
                            c.geometry.point_at(t).distance_to(&anchor)
                        })
                        .fold(0.0_f64, f64::max),
                    None => 0.0,
                }
            };
            if extent > gate {
                continue;
            }
            if popkern_core::trace::trace_enabled("boolean") {
                eprintln!(
                    "  [weld-stub] retiring self-looped record pair {i}/{twin} at n{node} (extent {extent:.2e})"
                );
            }
            for h in [i, twin] {
                let he_face = graph.half_edges[h].face;
                if let Some(f) = he_face
                    && let Some(edges) = graph.face_edges.get_mut(&f)
                {
                    edges.retain(|&e| e != h);
                }
                graph.half_edges[h].face = None;
            }
            // Records register their twin in face_a's list and their
            // fwd in face_b's list too — purge both lists wholesale.
            for edges in graph.face_edges.values_mut() {
                edges.retain(|&e| e != i && e != twin);
            }
            graph.nodes[node].outgoing.retain(|&e| e != i && e != twin);
        }
    }
}
