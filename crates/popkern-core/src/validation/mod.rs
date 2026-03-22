//! Geometric and topological validation for BRep models.

mod euler_poincare;
mod half_edge_checks;
mod manifold;

use crate::error::ValidationError;
use crate::topology::BRep;

impl BRep {
    /// Run all structural validation checks. Returns on the first failure.
    pub fn validate(&self) -> Result<(), ValidationError> {
        self.validate_euler_poincare()?;
        self.validate_half_edge_pairing()?;
        self.validate_half_edge_cycles()?;
        self.validate_vertex_consistency()?;
        self.validate_manifold_edges()?;
        self.validate_vertex_valence()?;
        Ok(())
    }

    /// Run all structural validation checks and collect all errors.
    pub fn validate_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();

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
}
