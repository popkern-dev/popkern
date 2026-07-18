//! Prism constructor.

use crate::error::CoreError;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

impl BRep {
    /// Build a prism by extruding a polygon profile along a direction.
    ///
    /// The profile must have at least 3 non-collinear points. The prism
    /// is created by extruding the profile face along `direction * distance`.
    ///
    /// # Errors
    ///
    /// Returns [`CoreError`] if the profile has fewer than 3 points,
    /// or if any Euler operation or validation fails.
    pub fn prism(
        profile: &[Point3],
        direction: Vec3,
        distance: f64,
        tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        if profile.len() < 3 {
            return Err(CoreError::PreconditionFailed(
                "prism profile must have at least 3 vertices".into(),
            ));
        }

        let mut brep = BRep::new(tolerance);

        // 1. mvfs at profile[0]
        let mvfs = brep.mvfs(profile[0])?;

        // 2. mev_at_vertex to profile[1]
        let mev1 = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, profile[1])?;

        // 3. mev for each subsequent vertex
        // We need to track the last twin half-edge to chain mev calls,
        // and track the first half-edge for the closing mef.
        let mut last_twin = mev1.half_edge_twin;
        for &point in &profile[2..] {
            let mev = brep.mev(last_twin, point)?;
            last_twin = mev.half_edge_twin;
        }

        // 4. mef to close the profile face
        let mef = brep.mef(last_twin, mev1.half_edge)?;

        // 5. extrude along direction * distance
        brep.extrude(mef.face, direction, distance, None)?;

        // Attach geometry to all edges and faces.
        brep.attach_geometry()?;

        // Attach pcurves on the planar faces. Without them every
        // downstream trim-domain consumer flies blind: the boolean
        // UV clipper cannot classify records against prism faces
        // (measured: a prism through a cylinder wall kept the full
        // unbounded plane×cylinder intersection lines and the union
        // was silent garbage — PK-150), and the tessellator samples
        // faces at unrelated densities (358 T-vertex welded opens on
        // a pristine prism box — PK-149).
        brep.attach_planar_pcurves()?;

        // Validate the result.
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}
