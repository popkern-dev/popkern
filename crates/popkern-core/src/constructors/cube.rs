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

        // Attach pcurves to half-edges on planar faces.
        brep.attach_planar_pcurves()?;

        // Validate the result.
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::curve::CurveEval;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn cube_half_edges_have_pcurves() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        for (_, he) in &cube.half_edges {
            assert!(
                he.pcurve.is_some(),
                "Every cube half-edge should have a pcurve"
            );
        }
    }

    #[test]
    fn cube_pcurve_endpoints_match_uv_projection() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        for (_, face) in &cube.faces {
            let sof = face.surface.as_ref().unwrap();
            let (plane_origin, u_axis, v_axis) = match &sof.geometry {
                SurfaceGeometry::Plane {
                    origin,
                    u_axis,
                    v_axis,
                    ..
                } => (*origin, *u_axis, *v_axis),
                _ => continue,
            };

            let u_dot_u = u_axis.dot(&u_axis);
            let v_dot_v = v_axis.dot(&v_axis);

            let outer_loop = face.outer_loop;
            let start_he = cube.loops[outer_loop].half_edge;
            let mut current = start_he;
            loop {
                let he = &cube.half_edges[current];
                let pcurve = he.pcurve.as_ref().unwrap();

                // Check start point.
                let origin_pos = cube.vertices[he.origin].position;
                let d_start = origin_pos - plane_origin;
                let u_start = d_start.dot(&u_axis) / u_dot_u;
                let v_start = d_start.dot(&v_axis) / v_dot_v;

                let pc_start = pcurve.point_at(0.0);
                assert!(
                    (pc_start.x - u_start).abs() < 1e-10,
                    "pcurve start u mismatch"
                );
                assert!(
                    (pc_start.y - v_start).abs() < 1e-10,
                    "pcurve start v mismatch"
                );

                // Check end point.
                let next_he = &cube.half_edges[he.next];
                let dest_pos = cube.vertices[next_he.origin].position;
                let d_end = dest_pos - plane_origin;
                let u_end = d_end.dot(&u_axis) / u_dot_u;
                let v_end = d_end.dot(&v_axis) / v_dot_v;

                let pc_end = pcurve.point_at(1.0);
                assert!((pc_end.x - u_end).abs() < 1e-10, "pcurve end u mismatch");
                assert!((pc_end.y - v_end).abs() < 1e-10, "pcurve end v mismatch");

                current = he.next;
                if current == start_he {
                    break;
                }
            }
        }
    }
}
