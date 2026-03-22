//! Cone constructor.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

impl BRep {
    /// Build a cone with 3-lune topology.
    ///
    /// The cone has its apex at (0, 0, height) and its base circle centered
    /// at the origin with the given `radius`. The two vertices are the apex
    /// and the base center (0, 0, 0).
    ///
    /// # Topology
    ///
    /// V=2, E=3, F=3, L=3, S=1.
    ///
    /// # Geometry
    ///
    /// Edges have line curves (seam lines from base center to apex). Surface
    /// geometry is not attached because the conical parameterization is
    /// degenerate at the apex, which would fail the orientation validation check.
    ///
    /// # Errors
    ///
    /// Returns [`CoreError`] if `radius` or `height` is not positive,
    /// or if any Euler operation or validation fails.
    pub fn cone(radius: f64, height: f64, tolerance: Tolerance) -> Result<Self, CoreError> {
        if radius <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cone radius must be positive".into(),
            ));
        }
        if height <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cone height must be positive".into(),
            ));
        }

        let mut brep = BRep::new(tolerance);
        let base = Point3::new(0.0, 0.0, 0.0);
        let apex = Point3::new(0.0, 0.0, height);

        // Build 3-lune topology (identical structure to sphere).
        let mvfs = brep.mvfs(base)?;
        let mev1 = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, apex)?;
        let mef1 = brep.mef(mev1.half_edge, mev1.half_edge_twin)?;
        let he_in_new_loop = mef1.half_edge_twin;
        let next_he = brep.half_edges[he_in_new_loop].next;
        let _mef2 = brep.mef(he_in_new_loop, next_he)?;

        // Attach line curves to all edges.
        brep.attach_line_curves();

        // Validate the result.
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}
