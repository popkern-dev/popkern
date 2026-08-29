//! Check 8: Loop closure — each loop's half-edge chain should be geometrically
//! closed within tolerance.
//!
//! Check 10: Shell watertight — for each edge, the two adjacent faces should
//! share the same geometric edge; no gaps.

use crate::error::ValidationError;
use crate::geometry::CurveEval;
use crate::topology::{BRep, TopologyRole};

impl BRep {
    /// Check 8: Loop closure.
    ///
    /// Walks each loop and verifies that the end of the last edge coincides
    /// with the start of the first edge, within `self.tolerance.linear()`.
    ///
    /// For edges without curve geometry, the check is based on vertex positions.
    /// For edges with [`CurveOnEdge`], the curve is evaluated at `t_end` and
    /// compared to the next edge's start.
    pub fn validate_loop_closure(&self) -> Result<(), ValidationError> {
        for (loop_handle, loop_data) in &self.loops {
            let start_he = loop_data.half_edge;

            // Collect vertex positions by walking the loop.
            // We use a bounded walk (at most half_edges.len() steps) to guard
            // against corrupted topology.
            let bound = self.half_edges.len().max(1);
            let mut cursor = start_he;
            let mut first_pos = None;
            let mut prev_end_pos = None;

            for _ in 0..=bound {
                let he = match self.half_edges.get(cursor) {
                    Some(h) => h,
                    None => break,
                };

                let origin_vertex = match self.vertices.get(he.origin) {
                    Some(v) => v,
                    None => break,
                };

                let start_pos = origin_vertex.position;
                if first_pos.is_none() {
                    first_pos = Some(start_pos);
                }

                // Determine the "end" position of this half-edge.
                // If there is curve geometry on the parent edge, evaluate at t_end.
                // Otherwise, the end is the next half-edge's origin vertex.
                let edge = match self.edges.get(he.edge) {
                    Some(e) => e,
                    None => break,
                };

                let end_pos = if let Some(curve_on_edge) = &edge.curve {
                    // The curve has two endpoints (t_start, t_end). One
                    // matches this half-edge's *origin* vertex; the other
                    // matches its destination. We can't assume the curve
                    // direction matches the canonical half-edge direction
                    // — boolean-engine outputs sometimes carry curves
                    // authored in the opposite direction. Pick the
                    // endpoint *farther* from `start_pos` (the half-edge's
                    // origin position) as this half-edge's "end".
                    let p_t_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
                    let p_t_end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);
                    let d_to_start = start_pos.distance_to(&p_t_start);
                    let d_to_end = start_pos.distance_to(&p_t_end);
                    if d_to_start <= d_to_end {
                        p_t_end
                    } else {
                        p_t_start
                    }
                } else {
                    // No curve: read the next half-edge's origin vertex.
                    let next_he = match self.half_edges.get(he.next) {
                        Some(n) => n,
                        None => break,
                    };
                    match self.vertices.get(next_he.origin) {
                        Some(v) => v.position,
                        None => break,
                    }
                };

                prev_end_pos = Some(end_pos);

                let next = he.next;
                cursor = next;
                if cursor == start_he {
                    break;
                }
            }

            // After the walk, the last end_pos must coincide with the first
            // start_pos (within linear tolerance).
            if let (Some(first), Some(last_end)) = (first_pos, prev_end_pos) {
                let gap = first.distance_to(&last_end);
                if gap > self.tolerance.linear() {
                    if crate::trace::trace_enabled("closure") {
                        eprintln!("  [closure] loop {loop_handle:?} gap {gap:.4}:");
                        let start_he = self.loops[loop_handle].half_edge;
                        let mut he_h = start_he;
                        for _ in 0..24 {
                            let he = &self.half_edges[he_h];
                            let v = self.vertices[he.origin].position;
                            let span = self.edges[he.edge].curve.as_ref().map(|c| {
                                use crate::geometry::curve::CurveEval as _;
                                let p0 = c.geometry.point_at(c.t_start);
                                let p1 = c.geometry.point_at(c.t_end);
                                format!(
                                    "({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) t[{:.4},{:.4}]",
                                    p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, c.t_start, c.t_end
                                )
                            });
                            eprintln!(
                                "    he edge={:?} origin=({:.3},{:.3},{:.3}) span={span:?}",
                                he.edge, v.x, v.y, v.z
                            );
                            he_h = he.next;
                            if he_h == start_he {
                                break;
                            }
                        }
                    }
                    return Err(ValidationError::OpenLoop {
                        loop_: loop_handle,
                        gap,
                    });
                }
            }
        }
        Ok(())
    }

    /// Check 10: Shell watertight.
    ///
    /// For each edge, the two adjacent half-edges belong to different faces.
    /// When curve geometry is present, evaluate the curve at its endpoints and
    /// verify both half-edges see the same positions (no gap).
    ///
    /// For v1a, gracefully skip edges with no geometry.
    pub fn validate_watertight(&self) -> Result<(), ValidationError> {
        for (edge_handle, edge_data) in &self.edges {
            // Retired scaffolding is not part of the production shell:
            // partition legs park displaced sides in ConstructionResidual
            // holders whose curves legitimately go stale, matching the
            // role skip in the naked-edge / degenerate / tessellation
            // detectors.
            if edge_data.role == TopologyRole::ConstructionResidual {
                continue;
            }
            let Some(curve_on_edge) = &edge_data.curve else {
                // No geometry — skip.
                continue;
            };

            // Evaluate the curve at its two endpoints.
            let p_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
            let p_end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);

            // The two half-edges for this edge must have origins that match
            // p_start and p_end respectively (or vice versa).
            let he = match self.half_edges.get(edge_data.half_edge) {
                Some(h) => h,
                None => continue,
            };
            let twin = match self.half_edges.get(he.twin) {
                Some(t) => t,
                None => continue,
            };

            let v_he = match self.vertices.get(he.origin) {
                Some(v) => v.position,
                None => continue,
            };
            let v_twin = match self.vertices.get(twin.origin) {
                Some(v) => v.position,
                None => continue,
            };

            // Edge curves are allowed to be authored in either direction.
            // Watertightness requires the two curve endpoints to match the two
            // half-edge origin vertices as an unordered pair.
            let natural_start = v_he.distance_to(&p_start);
            let natural_end = v_twin.distance_to(&p_end);
            let swapped_start = v_he.distance_to(&p_end);
            let swapped_end = v_twin.distance_to(&p_start);

            let natural_gap = natural_start.max(natural_end);
            let swapped_gap = swapped_start.max(swapped_end);
            let gap = natural_gap.min(swapped_gap);

            if gap > self.tolerance.linear() {
                return Err(ValidationError::ShellGap {
                    edge: edge_handle,
                    gap,
                });
            }
        }
        Ok(())
    }

    /// Find edges whose surface evaluations on either side disagree at
    /// interior points along the curve. `validate_watertight` only
    /// checks the two endpoints; OCCT-strict mesh slicers (Bambu
    /// Studio, FreeCAD) tessellate each pcurve to a polyline and
    /// require the resulting 3-D points to coincide between adjacent
    /// faces — interior disagreement appears in the imported mesh as
    /// "non-manifold edges."
    ///
    /// For each edge, samples `samples` points along the parameter
    /// range, evaluates `surface_a(pcurve_a(t))` and
    /// `surface_b(pcurve_b(t))`, and reports the edge with the max
    /// divergence if it exceeds `threshold`. The returned vector lists
    /// `(edge, max_divergence, midpoint_3d)` in descending order of
    /// divergence — useful for surfacing the worst offenders first.
    ///
    /// Edges without curve geometry, or with one or both half-edges
    /// missing pcurves, are skipped (we can't compute the divergence
    /// without both pcurves). Edges where both faces share the same
    /// surface (e.g. seam edges on cylinders) are also skipped.
    pub fn find_tessellation_non_manifold_edges(
        &self,
        samples: usize,
        threshold: f64,
    ) -> Vec<NonManifoldEdgeReport> {
        use crate::geometry::surface::SurfaceEval;
        let n_samples = samples.max(2);
        let mut reports: Vec<NonManifoldEdgeReport> = Vec::new();
        for (edge_handle, edge_data) in &self.edges {
            if edge_data.role == TopologyRole::ConstructionResidual {
                continue;
            }
            let Some(curve_on_edge) = edge_data.curve.as_ref() else {
                continue;
            };
            let Some(he) = self.half_edges.get(edge_data.half_edge) else {
                continue;
            };
            let Some(twin) = self.half_edges.get(he.twin) else {
                continue;
            };
            let Some(loop_a) = self.loops.get(he.loop_) else {
                continue;
            };
            let Some(loop_b) = self.loops.get(twin.loop_) else {
                continue;
            };
            let face_a = loop_a.face;
            let face_b = loop_b.face;
            if face_a == face_b {
                continue; // self-loop / seam edge
            }
            let face_a_is_residual = self
                .faces
                .get(face_a)
                .is_some_and(|face| face.role == TopologyRole::ConstructionResidual);
            let face_b_is_residual = self
                .faces
                .get(face_b)
                .is_some_and(|face| face.role == TopologyRole::ConstructionResidual);
            if face_a_is_residual || face_b_is_residual {
                continue;
            }
            let Some(pcurve_a) = he.pcurve.as_ref() else {
                continue;
            };
            let Some(pcurve_b) = twin.pcurve.as_ref() else {
                continue;
            };
            let Some(surf_a) = self.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let Some(surf_b) = self.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
                continue;
            };

            let (curve_t_min, curve_t_max) = (curve_on_edge.t_start, curve_on_edge.t_end);
            let (pa_min, pa_max) = pcurve_a.parameter_range();
            let (pb_min, pb_max) = pcurve_b.parameter_range();
            // The cyl-side and sliver-side pcurves run in opposite
            // parameter directions: one goes pcurve_a(t_min)→pcurve_a(t_max)
            // matching curve(t_min)→curve(t_max), the other reversed.
            // Detect direction by comparing 3-D positions at both ends.
            let curve_p_min = curve_on_edge.geometry.point_at(curve_t_min);
            let curve_p_max = curve_on_edge.geometry.point_at(curve_t_max);
            let p_a_at_min = surf_a
                .geometry
                .point_at(pcurve_a.point_at(pa_min).x, pcurve_a.point_at(pa_min).y);
            let p_a_at_max = surf_a
                .geometry
                .point_at(pcurve_a.point_at(pa_max).x, pcurve_a.point_at(pa_max).y);
            let pcurve_a_forward = curve_p_min.distance_to(&p_a_at_min)
                + curve_p_max.distance_to(&p_a_at_max)
                <= curve_p_min.distance_to(&p_a_at_max) + curve_p_max.distance_to(&p_a_at_min);
            let p_b_at_min = surf_b
                .geometry
                .point_at(pcurve_b.point_at(pb_min).x, pcurve_b.point_at(pb_min).y);
            let p_b_at_max = surf_b
                .geometry
                .point_at(pcurve_b.point_at(pb_max).x, pcurve_b.point_at(pb_max).y);
            let pcurve_b_forward = curve_p_min.distance_to(&p_b_at_min)
                + curve_p_max.distance_to(&p_b_at_max)
                <= curve_p_min.distance_to(&p_b_at_max) + curve_p_max.distance_to(&p_b_at_min);

            let mut max_div = 0.0_f64;
            let mut max_at = curve_p_min;
            for i in 0..n_samples {
                let s = i as f64 / (n_samples - 1) as f64;
                let t_curve = curve_t_min + s * (curve_t_max - curve_t_min);
                let s_a = if pcurve_a_forward { s } else { 1.0 - s };
                let s_b = if pcurve_b_forward { s } else { 1.0 - s };
                let t_a = pa_min + s_a * (pa_max - pa_min);
                let t_b = pb_min + s_b * (pb_max - pb_min);
                let uv_a = pcurve_a.point_at(t_a);
                let uv_b = pcurve_b.point_at(t_b);
                let p_a = surf_a.geometry.point_at(uv_a.x, uv_a.y);
                let p_b = surf_b.geometry.point_at(uv_b.x, uv_b.y);
                let div = p_a.distance_to(&p_b);
                if div > max_div {
                    max_div = div;
                    max_at = curve_on_edge.geometry.point_at(t_curve);
                }
            }
            if max_div > threshold {
                reports.push(NonManifoldEdgeReport {
                    edge: edge_handle,
                    face_a,
                    face_b,
                    max_divergence: max_div,
                    location: max_at,
                });
            }
        }
        reports.sort_by(|a, b| {
            b.max_divergence
                .partial_cmp(&a.max_divergence)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        reports
    }
}

/// One entry in `BRep::find_tessellation_non_manifold_edges`. The
/// `max_divergence` is the largest 3-D distance, sampled along the
/// edge curve, between `surface_a(pcurve_a(t))` and `surface_b(pcurve_b(t))`.
#[derive(Debug, Clone)]
pub struct NonManifoldEdgeReport {
    /// The edge being reported.
    pub edge: crate::topology::EdgeHandle,
    /// Face on the edge's primary half-edge side.
    pub face_a: crate::topology::FaceHandle,
    /// Face on the twin half-edge side.
    pub face_b: crate::topology::FaceHandle,
    /// The maximum 3-D distance between the two surface evaluations,
    /// sampled along the edge's parameter range.
    pub max_divergence: f64,
    /// 3-D location (on the edge's curve) where the maximum
    /// divergence was observed — useful for visualisation.
    pub location: crate::geometry::Point3,
}

#[cfg(test)]
#[path = "closure_tests.rs"]
mod tests;
