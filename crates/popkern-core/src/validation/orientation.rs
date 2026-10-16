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

            // Sample the normal at a representative interior parameter point.
            // We avoid (0, 0) because some surfaces (e.g. spheres) have
            // degenerate normals at parameter-space corners (poles).
            // For semi-infinite ranges (e.g. a cone's v in [0, ∞)), using
            // v=0 hits the apex (degenerate), so we fall back to v=1.0.
            let ((u0, u1), (v0, v1)) = surface_on_face.geometry.parameter_range();
            let u_mid = if u0.is_finite() && u1.is_finite() {
                0.5 * (u0 + u1)
            } else {
                0.0
            };
            let v_mid = if v0.is_finite() && v1.is_finite() {
                0.5 * (v0 + v1)
            } else if v0.is_finite() {
                // Semi-infinite upward (e.g. [0, ∞)): use v0 + 1.0 to stay off
                // any degenerate apex at v0.
                v0 + 1.0
            } else {
                0.0
            };
            // Degenerate when the chart's partials are parallel or
            // vanish, read scale-free: the squared length of the raw
            // normal against f64::EPSILON flagged every face under
            // ~1e-4 across (a 1.6 mm rounded cube's 70 um corner
            // patches: |Su x Sv|^2 ~ 2e-17), so a valid mm-scale part
            // failed validation.
            let geometry = &surface_on_face.geometry;
            let (su, sv) = (
                geometry.partial_u(u_mid, v_mid),
                geometry.partial_v(u_mid, v_mid),
            );
            let normal = geometry.normal_at(u_mid, v_mid);
            let scale = su.length() * sv.length();
            let degenerate = if scale > 0.0 && scale.is_finite() {
                su.cross(&sv).length() <= 1e-12 * scale
            } else {
                normal.length_squared() < f64::EPSILON
            };
            if degenerate || !normal.length().is_finite() {
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
        brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
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
                axis: crate::geometry::Vec3::Z,
            },
            orientation: true,
        });
        let result = brep.validate_face_orientation();
        assert!(result.is_err());
    }

    /// Check 7 reads degeneracy scale-free: a 10 um sphere's chart
    /// normal has squared length ~1e-20, which the old absolute
    /// f64::EPSILON test called degenerate.
    #[test]
    fn a_tiny_sphere_face_is_not_degenerate() {
        use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        let mut brep = make_cube();
        let face_handle = brep.faces.keys().next().unwrap();
        brep.faces[face_handle].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1e-5,
                axis: crate::geometry::Vec3::Z,
            },
            orientation: true,
        });
        assert!(brep.validate_face_orientation().is_ok());
    }
}
