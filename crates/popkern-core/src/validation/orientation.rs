//! Check 7: Face orientation — surface normals should not be degenerate.

use crate::error::ValidationError;
use crate::geometry::SurfaceEval;
use crate::topology::BRep;

impl BRep {
    /// Check 7: Face orientation.
    ///
    /// For faces with no surface geometry, this check is skipped (returns Ok).
    /// When a surface IS present, verify that the surface normal at the
    /// parameter-space origin is non-zero (degenerate normal = bad geometry).
    ///
    /// Full consistency checking between adjacent faces is deferred to a future
    /// version; the stub here returns Ok after the per-face check.
    pub fn validate_face_orientation(&self) -> Result<(), ValidationError> {
        for (face_handle, face_data) in &self.faces {
            let Some(surface_on_face) = &face_data.surface else {
                // No geometry attached — skip gracefully.
                continue;
            };

            // Sample the normal at the midpoint of the parameter range.
            // We avoid (0, 0) because some surfaces (e.g. spheres) have
            // degenerate normals at parameter-space corners (poles).
            let ((u0, u1), (v0, v1)) = surface_on_face.geometry.parameter_range();
            let u_mid = if u0.is_finite() && u1.is_finite() {
                0.5 * (u0 + u1)
            } else {
                0.0
            };
            let v_mid = if v0.is_finite() && v1.is_finite() {
                0.5 * (v0 + v1)
            } else {
                0.0
            };
            let normal = surface_on_face.geometry.normal_at(u_mid, v_mid);
            if normal.length_squared() < f64::EPSILON {
                return Err(ValidationError::InconsistentOrientation { face: face_handle });
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
        brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0)
            .unwrap();
        brep
    }

    #[test]
    fn topology_only_cube_passes_orientation() {
        let brep = make_cube();
        assert!(brep.validate_face_orientation().is_ok());
    }

    #[test]
    fn cube_with_surface_geometry_passes_orientation() {
        let cube =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
        assert!(cube.validate_face_orientation().is_ok());
    }

    #[test]
    fn degenerate_surface_normal_fails_orientation() {
        use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        let mut brep = make_cube();
        // Attach a degenerate sphere surface (radius 0) to a face.
        // Normal at (0,0) will be zero for a zero-radius sphere.
        let face_handle = brep.faces.keys().next().unwrap();
        brep.faces[face_handle].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 0.0,
            },
            orientation: true,
        });
        let result = brep.validate_face_orientation();
        assert!(result.is_err());
    }
}
