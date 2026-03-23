//! Sphere constructor.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

impl BRep {
    /// Build a sphere with 3-lune topology.
    ///
    /// The sphere is centered at the origin with the given `radius`.
    /// The two vertices (poles) are at (0, 0, radius) and (0, 0, -radius).
    /// Three edges run between the poles, dividing the sphere into 3 lune faces.
    ///
    /// # Topology
    ///
    /// V=2, E=3, F=3, L=3, S=1.
    ///
    /// # Geometry
    ///
    /// All three faces are assigned `SurfaceGeometry::Sphere` with the same
    /// center and radius. Edges have line curves (seam lines from pole to pole);
    /// proper meridian arc curves are deferred to a future version.
    ///
    /// # Errors
    ///
    /// Returns [`CoreError`] if `radius` is not positive,
    /// or if any Euler operation or validation fails.
    pub fn sphere(radius: f64, tolerance: Tolerance) -> Result<Self, CoreError> {
        if radius <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "sphere radius must be positive".into(),
            ));
        }

        let mut brep = BRep::new(tolerance);
        let south = Point3::new(0.0, 0.0, -radius);
        let north = Point3::new(0.0, 0.0, radius);

        // 1. mvfs at south pole
        let mvfs = brep.mvfs(south)?;

        // 2. mev_at_vertex to north pole — edge 1
        let mev1 = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, north)?;

        // 3. mef to create edge 2 and face 2
        // Loop has: he1 (south->north), he1_twin (north->south)
        // mef(he1, he1_twin) connects origin(he1)=south to origin(he1_twin)=north
        let mef1 = brep.mef(mev1.half_edge, mev1.half_edge_twin)?;

        // 4. mef to create edge 3 and face 3
        // One of the loops now has 2 half-edges. We need to find two half-edges
        // in the same loop that originate at different vertices.
        // mef1.half_edge_twin is in the new loop and originates at north.
        // Following next from mef1.half_edge_twin should give us a he originating at south.
        let he_in_new_loop = mef1.half_edge_twin;
        let next_he = brep.half_edges[he_in_new_loop].next;
        let _mef2 = brep.mef(he_in_new_loop, next_he)?;

        // Attach line curves to all edges.
        brep.attach_line_curves();

        // Attach sphere surface geometry to every face.
        let sphere_geom = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius,
        };
        let face_handles: Vec<_> = brep.faces.keys().collect();
        for fh in face_handles {
            brep.faces[fh].surface = Some(SurfaceOnFace {
                geometry: sphere_geom.clone(),
                orientation: true,
            });
        }

        // Validate the result.
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}
