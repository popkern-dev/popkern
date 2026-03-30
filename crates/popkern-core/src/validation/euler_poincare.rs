//! Euler-Poincaré invariant check.

use crate::error::ValidationError;
use crate::topology::BRep;

impl BRep {
    /// Validate the Euler-Poincaré formula for this BRep.
    ///
    /// For genus-0 models: `V - E + F - (L - F) = 2 * S`
    ///
    /// This simplifies from the full formula `V - E + F - (L - F) - 2(S - G) = 0`
    /// with genus G = 0.
    pub fn validate_euler_poincare(&self) -> Result<(), ValidationError> {
        let v = self.vertices.len() as i64;
        let e = self.edges.len() as i64;
        let f = self.faces.len() as i64;
        let l = self.loops.len() as i64;
        let s = self.shells.len() as i64;

        // V - E + F - (L - F) = 2 * S
        let lhs = v - e + f - (l - f);
        let rhs = 2 * s;

        if lhs != rhs {
            // Use the first solid handle for the error, or default if none.
            let solid = self.solids.keys().next().unwrap_or_default();

            return Err(ValidationError::EulerPoincare {
                solid,
                expected: rhs,
                actual: lhs,
            });
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// Helper: build a cube via mvfs + mev + mef + extrude.
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
    fn cube_passes_euler_poincare() {
        let brep = make_cube();
        assert!(brep.validate_euler_poincare().is_ok());
    }

    #[test]
    fn mvfs_only_passes_euler_poincare() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // V=1, E=0, F=1, L=1, S=1 => 1 - 0 + 1 - (1 - 1) = 2 = 2*1
        assert!(brep.validate_euler_poincare().is_ok());
    }

    #[test]
    fn mvfs_plus_mev_passes_euler_poincare() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let _mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // V=2, E=1, F=1, L=1, S=1 => 2 - 1 + 1 - (1 - 1) = 2 = 2*1
        assert!(brep.validate_euler_poincare().is_ok());
    }

    #[test]
    fn corrupted_counts_fail_euler_poincare() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // Remove the face but keep vertex — breaks Euler-Poincaré
        brep.faces.remove(mvfs.face);
        assert!(brep.validate_euler_poincare().is_err());
    }
}
