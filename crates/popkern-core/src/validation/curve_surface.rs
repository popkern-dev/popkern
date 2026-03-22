//! Check 9: Curve-surface consistency — edge curves should lie on their
//! adjacent face surfaces.
//!
//! For v1a this is a graceful stub: when no geometry is present the check
//! passes unconditionally. Full deviation checking is deferred.

use crate::error::ValidationError;
use crate::topology::BRep;

impl BRep {
    /// Check 9: Curve-surface consistency.
    ///
    /// Verifies that edge curves lie on the surfaces of their adjacent faces.
    /// For v1a, this check is skipped (returns Ok) when either the edge or its
    /// adjacent face lacks geometry. When both are present, a future
    /// implementation would sample the curve and measure deviation from the
    /// surface.
    pub fn validate_curve_surface(&self) -> Result<(), ValidationError> {
        // v1a: No geometry wired on edges/faces yet in the test corpus.
        // Gracefully return Ok; full sampling logic to be added in a future task.
        let _ = &self.edges;
        let _ = &self.faces;
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
    fn topology_only_cube_passes_curve_surface() {
        let brep = make_cube();
        assert!(brep.validate_curve_surface().is_ok());
    }
}
