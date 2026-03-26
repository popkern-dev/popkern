//! Pcurve recomputation for BRep faces.
//!
//! Projects edge curve endpoints (or vertex positions) to surface UV space
//! and creates degree-1 NURBS pcurves for half-edges that lack them.

use super::FaceHandle;
use super::brep::BRep;

impl BRep {
    /// Recompute pcurves for all half-edges on a face by projecting
    /// edge curve endpoints to the face surface's UV space.
    ///
    /// Skips half-edges that already have pcurves. Falls back to vertex
    /// positions for edges that lack curve geometry. Clamps UV coordinates
    /// to the surface's finite parameter range.
    pub fn recompute_face_pcurves(
        &mut self,
        face_handle: FaceHandle,
    ) -> Result<(), crate::error::CoreError> {
        use crate::geometry::curve::CurveEval;
        use crate::geometry::surface::SurfaceEval;
        use crate::geometry::uv_projection::project_point_to_uv;
        use crate::geometry::{NurbsCurve, Point3};

        let face = self
            .faces
            .get(face_handle)
            .ok_or_else(|| crate::error::CoreError::InvalidHandle("face not found".into()))?;
        let surface = match &face.surface {
            Some(sof) => sof.geometry.clone(),
            None => return Ok(()),
        };

        let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

        // Collect all loop handles for this face.
        let outer_loop = face.outer_loop;
        let inner_loops = face.inner_loops.clone();

        let mut all_loop_handles = vec![outer_loop];
        all_loop_handles.extend(inner_loops);

        // Collect half-edge handles from all loops before mutating.
        let mut he_handles = Vec::new();
        for lh in &all_loop_handles {
            let loop_data = match self.loops.get(*lh) {
                Some(l) => l,
                None => continue,
            };
            let start_he = loop_data.half_edge;
            if !self.half_edges.contains_key(start_he) {
                continue;
            }
            he_handles.push(start_he);
            let mut current = self.half_edges[start_he].next;
            while current != start_he {
                he_handles.push(current);
                current = self.half_edges[current].next;
            }
        }

        // For each half-edge without a pcurve, compute one.
        for heh in he_handles {
            if self.half_edges[heh].pcurve.is_some() {
                continue;
            }

            let he_origin_pos = self.vertices[self.half_edges[heh].origin].position;
            let next_heh = self.half_edges[heh].next;
            let he_dest_pos = self.vertices[self.half_edges[next_heh].origin].position;

            let edge_handle = self.half_edges[heh].edge;
            let (start_3d, end_3d) = if let Some(ref curve_on_edge) = self.edges[edge_handle].curve
            {
                // Determine direction: compare he.origin position with curve start.
                let curve_start_pos = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
                let curve_end_pos = curve_on_edge.geometry.point_at(curve_on_edge.t_end);

                if he_origin_pos.distance_to(&curve_start_pos)
                    <= he_origin_pos.distance_to(&curve_end_pos)
                {
                    (curve_start_pos, curve_end_pos)
                } else {
                    (curve_end_pos, curve_start_pos)
                }
            } else {
                // No edge curve: fall back to vertex positions.
                (he_origin_pos, he_dest_pos)
            };

            // Project 3D start/end to UV.
            let (u_start, v_start) = match project_point_to_uv(&surface, &start_3d) {
                Some(uv) => uv,
                None => continue,
            };
            let (u_end, v_end) = match project_point_to_uv(&surface, &end_3d) {
                Some(uv) => uv,
                None => continue,
            };

            // Clamp UV to surface parameter range (if finite).
            let clamp_if_finite = |val: f64, lo: f64, hi: f64| -> f64 {
                if lo.is_finite() && hi.is_finite() {
                    val.clamp(lo, hi)
                } else {
                    val
                }
            };
            let u_s = clamp_if_finite(u_start, u_min, u_max);
            let v_s = clamp_if_finite(v_start, v_min, v_max);
            let u_e = clamp_if_finite(u_end, u_min, u_max);
            let v_e = clamp_if_finite(v_end, v_min, v_max);

            // Create degree-1 NURBS line pcurve.
            let pcurve = NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(u_s, v_s, 0.0), Point3::new(u_e, v_e, 0.0)],
                weights: vec![1.0, 1.0],
            };

            self.half_edges[heh].pcurve = Some(pcurve);
        }

        Ok(())
    }

    /// Recompute pcurves for ALL faces in the BRep.
    pub fn recompute_all_pcurves(&mut self) -> Result<(), crate::error::CoreError> {
        let faces: Vec<_> = self.faces.keys().collect();
        for fh in faces {
            self.recompute_face_pcurves(fh)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn recompute_pcurves_on_cube_is_idempotent() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // All half-edges already have pcurves from unit_cube.
        brep.recompute_all_pcurves().unwrap();
        for (_, he) in &brep.half_edges {
            assert!(he.pcurve.is_some(), "pcurve should still be present");
        }
    }

    #[test]
    fn recompute_fills_missing_pcurves() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // Strip all pcurves.
        let he_keys: Vec<_> = brep.half_edges.keys().collect();
        for heh in &he_keys {
            brep.half_edges[*heh].pcurve = None;
        }
        // Recompute should fill them back.
        brep.recompute_all_pcurves().unwrap();
        for (_, he) in &brep.half_edges {
            assert!(he.pcurve.is_some(), "pcurve should be recomputed");
        }
    }
}
