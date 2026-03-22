//! Geometric and topological validation for BRep models.

mod closure;
mod curve_surface;
mod euler_poincare;
mod half_edge_checks;
mod manifold;
mod orientation;

use crate::error::ValidationError;
use crate::topology::BRep;

impl BRep {
    /// Run all validation checks (structural + geometric). Returns on the first failure.
    pub fn validate(&self) -> Result<(), ValidationError> {
        // Structural checks (1-6)
        self.validate_euler_poincare()?;
        self.validate_half_edge_pairing()?;
        self.validate_half_edge_cycles()?;
        self.validate_vertex_consistency()?;
        self.validate_manifold_edges()?;
        self.validate_vertex_valence()?;
        // Geometric checks (7-10)
        self.validate_face_orientation()?;
        self.validate_loop_closure()?;
        self.validate_curve_surface()?;
        self.validate_watertight()?;
        Ok(())
    }

    /// Run all validation checks and collect all errors.
    pub fn validate_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();

        // Structural checks (1-6)
        if let Err(e) = self.validate_euler_poincare() {
            errors.push(e);
        }
        if let Err(e) = self.validate_half_edge_pairing() {
            errors.push(e);
        }
        if let Err(e) = self.validate_half_edge_cycles() {
            errors.push(e);
        }
        if let Err(e) = self.validate_vertex_consistency() {
            errors.push(e);
        }
        if let Err(e) = self.validate_manifold_edges() {
            errors.push(e);
        }
        if let Err(e) = self.validate_vertex_valence() {
            errors.push(e);
        }
        // Geometric checks (7-10)
        if let Err(e) = self.validate_face_orientation() {
            errors.push(e);
        }
        if let Err(e) = self.validate_loop_closure() {
            errors.push(e);
        }
        if let Err(e) = self.validate_curve_surface() {
            errors.push(e);
        }
        if let Err(e) = self.validate_watertight() {
            errors.push(e);
        }

        errors
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
    fn cube_passes_full_validation() {
        let brep = make_cube();
        assert!(brep.validate().is_ok());
    }

    #[test]
    fn cube_validate_all_returns_no_errors() {
        let brep = make_cube();
        assert!(brep.validate_all().is_empty());
    }

    #[test]
    fn corrupted_brep_validate_all_collects_errors() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Corrupt the twin pointer.
        brep.half_edges[mev.half_edge].twin = Default::default();
        let errors = brep.validate_all();
        // Should have at least 1 error (broken twin).
        assert!(!errors.is_empty());
    }

    #[test]
    fn mvfs_only_passes_validation() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        assert!(brep.validate().is_ok());
    }

    // --- Task 20: geometric checks ---

    #[test]
    fn geometric_checks_pass_on_topology_only_brep() {
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

        assert!(brep.validate_face_orientation().is_ok());
        assert!(brep.validate_loop_closure().is_ok());
        assert!(brep.validate_curve_surface().is_ok());
        assert!(brep.validate_watertight().is_ok());
    }

    #[test]
    fn full_validate_passes_for_cube() {
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

        assert!(brep.validate().is_ok());
    }

    #[test]
    fn validate_all_collects_multiple_errors() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Corrupt twin AND next to trigger multiple validation errors
        brep.half_edges[mev.half_edge].twin = Default::default();
        brep.half_edges[mev.half_edge_twin].twin = Default::default();
        let errors = brep.validate_all();
        // Should have multiple errors (broken twin, broken cycle, manifold, etc.)
        assert!(
            errors.len() >= 2,
            "Expected at least 2 errors, got {}",
            errors.len()
        );
    }

    #[test]
    fn validate_all_hits_all_check_branches() {
        use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        // Build a cube with geometry, then corrupt it to fail many checks
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();

        // Remove a face to break Euler-Poincare (check 1)
        let face_h = brep.faces.keys().next().unwrap();
        brep.faces.remove(face_h);

        // Corrupt a twin pointer to break pairing (check 2)
        let he_h = brep.half_edges.keys().next().unwrap();
        brep.half_edges[he_h].twin = Default::default();

        // Attach degenerate surface to trigger orientation failure (check 7)
        let face2_h = brep.faces.keys().next().unwrap();
        brep.faces[face2_h].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 0.0,
            },
            orientation: true,
        });

        let errors = brep.validate_all();
        // Should have errors from multiple check categories
        assert!(
            errors.len() >= 3,
            "Expected at least 3 errors from different checks, got {}",
            errors.len()
        );
    }

    #[test]
    fn validate_all_on_cube_with_surface_geometry() {
        let cube =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
        let errors = cube.validate_all();
        assert!(
            errors.is_empty(),
            "Cube with geometry should pass all validation: {:?}",
            errors
        );
    }

    #[test]
    fn loop_closure_from_vertex_positions() {
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
}
