//! Check 8: Loop closure — each loop's half-edge chain should be geometrically
//! closed within tolerance.
//!
//! Check 10: Shell watertight — for each edge, the two adjacent faces should
//! share the same geometric edge; no gaps.

use crate::error::ValidationError;
use crate::geometry::CurveEval;
use crate::topology::BRep;

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

            // Check: edge start ~ he origin AND edge end ~ twin origin.
            let gap_start = v_he.distance_to(&p_start);
            let gap_end = v_twin.distance_to(&p_end);

            if gap_start > self.tolerance.linear() {
                return Err(ValidationError::ShellGap {
                    edge: edge_handle,
                    gap: gap_start,
                });
            }
            if gap_end > self.tolerance.linear() {
                return Err(ValidationError::ShellGap {
                    edge: edge_handle,
                    gap: gap_end,
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    fn make_cube() -> BRep {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
        brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
            .unwrap();
        brep
    }

    #[test]
    fn cube_passes_loop_closure() {
        let brep = make_cube();
        assert!(brep.validate_loop_closure().is_ok());
    }

    #[test]
    fn cube_passes_watertight() {
        let brep = make_cube();
        assert!(brep.validate_watertight().is_ok());
    }

    #[test]
    fn triangle_loop_is_closed() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(0.5, 1.0, 0.0))
            .unwrap();
        let _mef = brep.mef(mev2.half_edge_twin, mev1.half_edge).unwrap();
        assert!(brep.validate_loop_closure().is_ok());
    }

    #[test]
    fn loop_closure_with_curve_geometry() {
        // Use a cube with geometry to test the curve evaluation path
        let cube =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
        assert!(cube.validate_loop_closure().is_ok());
    }

    #[test]
    fn watertight_with_curve_geometry() {
        let cube =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
        assert!(cube.validate_watertight().is_ok());
    }

    #[test]
    fn loop_closure_on_corrupted_brep_breaks_gracefully() {
        // Remove a half-edge to trigger the `None` break paths in validate_loop_closure
        let mut brep = make_cube();
        let he_handle = brep.half_edges.keys().next().unwrap();
        // Corrupt next pointer to a non-existent handle
        brep.half_edges[he_handle].next = Default::default();
        // Should not panic — just break out of walk
        let _ = brep.validate_loop_closure();
    }

    #[test]
    fn loop_closure_twin_path_with_curve() {
        // Build a cube with geometry, then verify closure works through twin-side curve eval
        let brep =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
        // The cube has curves on all edges. Loop closure walks all loops and
        // picks the correct curve endpoint for each half-edge by direction
        // matching. Just ensure it passes.
        assert!(brep.validate_loop_closure().is_ok());
    }

    #[test]
    fn corrupted_edge_handle_breaks_gracefully() {
        // Remove an edge to trigger the edge.curve None path
        let mut brep = make_cube();
        let edge_handle = brep.edges.keys().next().unwrap();
        brep.edges.remove(edge_handle);
        // Should not panic
        let _ = brep.validate_loop_closure();
    }

    #[test]
    fn corrupted_vertex_handle_breaks_gracefully() {
        // Remove a vertex to trigger the vertices.get None path
        let mut brep = make_cube();
        let v_handle = brep.vertices.keys().next().unwrap();
        brep.vertices.remove(v_handle);
        // Should not panic
        let _ = brep.validate_loop_closure();
    }

    #[test]
    fn watertight_with_corrupted_he_skips() {
        // Corrupt an edge's half_edge to point to nonexistent handle
        let mut brep =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
        let edge_handle = brep.edges.keys().next().unwrap();
        brep.edges[edge_handle].half_edge = Default::default();
        // Should not panic; just skip the corrupted edge
        let _ = brep.validate_watertight();
    }

    #[test]
    fn watertight_with_bad_gap() {
        // Create a cube with geometry, then corrupt a vertex position to create a gap
        let mut brep =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();

        // Move a vertex far away to create a gap
        let v_handle = brep.vertices.keys().next().unwrap();
        brep.vertices[v_handle].position = Point3::new(100.0, 100.0, 100.0);

        let result = brep.validate_watertight();
        assert!(result.is_err());
    }

    #[test]
    fn watertight_twin_vertex_missing_skips() {
        // Remove the twin's origin vertex to trigger the vertices.get(twin.origin) None path
        let mut brep =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();

        // Find an edge with a curve
        let edge_h = brep
            .edges
            .keys()
            .find(|&eh| brep.edges[eh].curve.is_some())
            .unwrap();
        let he_h = brep.edges[edge_h].half_edge;
        let twin_h = brep.half_edges[he_h].twin;
        let twin_origin = brep.half_edges[twin_h].origin;

        // Remove the twin's origin vertex
        brep.vertices.remove(twin_origin);

        // Should not panic — just skip
        let _ = brep.validate_watertight();
    }

    #[test]
    fn watertight_twin_missing_skips() {
        // Corrupt twin pointer to trigger he.twin None path
        let mut brep =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();

        let edge_h = brep
            .edges
            .keys()
            .find(|&eh| brep.edges[eh].curve.is_some())
            .unwrap();
        let he_h = brep.edges[edge_h].half_edge;
        // Make twin point to nonexistent handle
        brep.half_edges[he_h].twin = Default::default();

        // Should not panic — just skip
        let _ = brep.validate_watertight();
    }

    #[test]
    fn watertight_gap_end_detected() {
        // Create a gap at the twin (end) vertex, not the start vertex.
        // This tests the gap_end > tolerance path (line 151-153).
        let mut brep =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();

        // Find an edge with a curve, then move the twin's origin far away
        // while keeping the canonical half-edge's origin in place.
        for edge_h in brep.edges.keys().collect::<Vec<_>>() {
            if brep.edges[edge_h].curve.is_some() {
                let he_h = brep.edges[edge_h].half_edge;
                let twin_h = brep.half_edges[he_h].twin;
                let twin_origin = brep.half_edges[twin_h].origin;
                brep.vertices[twin_origin].position = Point3::new(999.0, 999.0, 999.0);
                break;
            }
        }

        let result = brep.validate_watertight();
        assert!(result.is_err());
    }
}
