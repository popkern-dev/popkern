//! Self-intersection validation check (#11).
//!
//! Detects when non-adjacent faces in a solid intersect each other.

use crate::error::ValidationError;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::BRep;
use std::collections::HashSet;

impl BRep {
    /// Check #11: no face intersects another non-adjacent face.
    ///
    /// The check proceeds in three stages:
    ///
    /// 1. Build a set of adjacent face-index pairs by walking shared half-edges.
    /// 2. Sample each face surface on a 5×5 grid to compute an AABB.
    /// 3. For surviving non-adjacent pairs whose AABBs overlap, sample points on
    ///    each face and test whether any sample from face A lies within tolerance
    ///    of face B's surface (using point-to-surface distance).
    ///
    /// For v1f-2 the point-proximity test is used only for `Plane`-surface faces,
    /// where the distance from a point to the plane is analytic and exact.
    /// Other surface types skip the detailed test (the AABB pre-filter alone is not
    /// enough to report an error).
    pub fn validate_self_intersection(&self) -> Result<(), ValidationError> {
        let face_keys: Vec<_> = self.faces.keys().collect();
        let n = face_keys.len();
        if n < 2 {
            return Ok(());
        }

        // ------------------------------------------------------------------ //
        // 1. Build adjacency set (by face index pairs, not handles)           //
        // ------------------------------------------------------------------ //
        let mut adjacent: HashSet<(usize, usize)> = HashSet::new();

        for (i, &fk) in face_keys.iter().enumerate() {
            let face = &self.faces[fk];
            let start_he_handle = self.loops[face.outer_loop].half_edge;
            // Empty loop — skip
            if !self.half_edges.contains_key(start_he_handle) {
                continue;
            }
            let mut he_handle = start_he_handle;
            // Walk the outer loop; guard against cycles larger than the mesh.
            let bound = self.half_edges.len().max(1);
            for _ in 0..=bound {
                let he = match self.half_edges.get(he_handle) {
                    Some(h) => h,
                    None => break,
                };
                let twin_handle = he.twin;
                if let Some(twin) = self.half_edges.get(twin_handle) {
                    let twin_loop_handle = twin.loop_;
                    if let Some(twin_loop) = self.loops.get(twin_loop_handle) {
                        let twin_face = twin_loop.face;
                        if let Some(j) =
                            face_keys.iter().position(|&k| k == twin_face)
                        {
                            if i != j {
                                let pair = if i < j { (i, j) } else { (j, i) };
                                adjacent.insert(pair);
                            }
                        }
                    }
                }
                he_handle = he.next;
                if he_handle == start_he_handle {
                    break;
                }
            }
        }

        // ------------------------------------------------------------------ //
        // 2. AABB pre-filter: compute bounding boxes from surface samples     //
        // ------------------------------------------------------------------ //
        let sample_params = [0.0, 0.25, 0.5, 0.75, 1.0];

        let aabbs: Vec<Option<Aabb>> = face_keys
            .iter()
            .map(|&fk| {
                let face = &self.faces[fk];
                face.surface.as_ref().and_then(|s| {
                    let ((u0, u1), (v0, v1)) = s.geometry.parameter_range();
                    let mut pts = Vec::with_capacity(25);
                    for &tp in &sample_params {
                        let u = u0 + tp * (u1 - u0);
                        for &tq in &sample_params {
                            let v = v0 + tq * (v1 - v0);
                            pts.push(s.geometry.point_at(u, v));
                        }
                    }
                    Aabb::from_points(&pts)
                })
            })
            .collect();

        // ------------------------------------------------------------------ //
        // 3. Test non-adjacent, AABB-overlapping pairs                       //
        // ------------------------------------------------------------------ //
        let linear_tol = self.tolerance.linear;

        for i in 0..n {
            for j in (i + 1)..n {
                let pair = (i, j);
                if adjacent.contains(&pair) {
                    continue;
                }

                // AABB gate: skip if either face has no geometry or no overlap.
                let (aabb_i, aabb_j) = match (&aabbs[i], &aabbs[j]) {
                    (Some(a), Some(b)) => (a, b),
                    _ => continue,
                };
                if !aabb_i.overlaps(aabb_j, linear_tol) {
                    continue;
                }

                // Detailed analytic check for Plane–Plane pairs.
                // For other surface types the AABB overlap is a conservative
                // approximation; we skip the expensive check and rely on the
                // caller to run a dedicated intersect_surfaces step later.
                let face_i = &self.faces[face_keys[i]];
                let face_j = &self.faces[face_keys[j]];

                let surf_i = match face_i.surface.as_ref() {
                    Some(s) => s,
                    None => continue,
                };
                let surf_j = match face_j.surface.as_ref() {
                    Some(s) => s,
                    None => continue,
                };

                if let (
                    SurfaceGeometry::Plane {
                        origin: orig_j,
                        normal: nrm_j,
                        ..
                    },
                    SurfaceGeometry::Plane {
                        origin: orig_i,
                        normal: nrm_i,
                        ..
                    },
                ) = (&surf_j.geometry, &surf_i.geometry)
                {
                    // Sample points on face i; check if any lie within tolerance
                    // of the plane of face j.
                    let ((u0i, u1i), (v0i, v1i)) = surf_i.geometry.parameter_range();
                    for &tp in &sample_params {
                        let u = u0i + tp * (u1i - u0i);
                        for &tq in &sample_params {
                            let v = v0i + tq * (v1i - v0i);
                            let p = surf_i.geometry.point_at(u, v);
                            // Signed distance from p to plane j.
                            let d = (p - *orig_j).dot(nrm_j).abs();
                            if d <= linear_tol {
                                // p lies on (within tol) face j's plane.
                                // Report self-intersection.
                                return Err(ValidationError::SelfIntersection {
                                    face_a: face_keys[i],
                                    face_b: face_keys[j],
                                });
                            }
                        }
                    }

                    // Also sample face j and check against plane i.
                    let ((u0j, u1j), (v0j, v1j)) = surf_j.geometry.parameter_range();
                    for &tp in &sample_params {
                        let u = u0j + tp * (u1j - u0j);
                        for &tq in &sample_params {
                            let v = v0j + tq * (v1j - v0j);
                            let p = surf_j.geometry.point_at(u, v);
                            let d = (p - *orig_i).dot(nrm_i).abs();
                            if d <= linear_tol {
                                return Err(ValidationError::SelfIntersection {
                                    face_a: face_keys[i],
                                    face_b: face_keys[j],
                                });
                            }
                        }
                    }
                }
                // Non-Plane pairs: AABB overlap is not sufficient to assert
                // intersection — skip to avoid false positives.
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

    #[test]
    fn cube_no_self_intersection() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        assert!(cube.validate_self_intersection().is_ok());
    }

    #[test]
    fn prism_no_self_intersection() {
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
        ];
        let prism = BRep::prism(
            &profile,
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            Tolerance::default(),
        )
        .unwrap();
        assert!(prism.validate_self_intersection().is_ok());
    }

    #[test]
    fn empty_brep_passes() {
        let brep = BRep::new(Tolerance::default());
        assert!(brep.validate_self_intersection().is_ok());
    }

    #[test]
    fn single_face_brep_passes() {
        // A brep with only one face cannot self-intersect.
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
        let _mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
        assert!(brep.validate_self_intersection().is_ok());
    }
}
