//! Unit cube constructor.

use crate::error::CoreError;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

impl BRep {
    /// Build a unit cube [0,1]^3 using Euler operators with analytic geometry.
    ///
    /// # Errors
    ///
    /// Returns [`CoreError`] if any Euler operation or validation fails.
    pub fn unit_cube(tolerance: Tolerance) -> Result<Self, CoreError> {
        let mut brep = BRep::new(tolerance);

        // 1. mvfs at (0,0,0)
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0))?;

        // 2. mev_at_vertex to (1,0,0)
        let mev1 = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))?;

        // 3. mev to (1,1,0)
        let mev2 = brep.mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))?;

        // 4. mev to (0,1,0)
        let mev3 = brep.mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))?;

        // 5. mef to close the bottom face (connect v3 back to v0)
        let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge)?;

        // 6. extrude along (0,0,1) by 1.0
        brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0)?;

        // Attach geometry to all edges and faces.
        brep.attach_geometry()?;

        // Validate the result.
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}
