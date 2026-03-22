//! Convenience constructors for common solids.
//!
//! Build topology via Euler operators, then attach analytic geometry
//! (line curves and planar surfaces) to all edges and faces.

use crate::error::CoreError;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Compute a plane surface from 3 non-collinear points.
fn plane_from_points(p0: Point3, p1: Point3, p2: Point3) -> Result<SurfaceGeometry, CoreError> {
    let u_axis = p1 - p0;
    let v_axis = p2 - p0;
    let normal = u_axis.cross(&v_axis);
    let normal = normal
        .normalized()
        .ok_or_else(|| CoreError::DegenerateSurface("collinear points for plane".into()))?;
    Ok(SurfaceGeometry::Plane {
        origin: p0,
        normal,
        u_axis,
        v_axis,
    })
}

/// Create a line `CurveOnEdge` between two points.
fn line_between(start: Point3, end: Point3) -> CurveOnEdge {
    let direction = end - start;
    CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: start,
            direction,
        },
        t_start: 0.0,
        t_end: 1.0,
    }
}

impl BRep {
    /// Attach line curves to all edges and planar surfaces to all faces.
    ///
    /// For each edge, computes the line from the origin of one half-edge
    /// to the origin of its twin. For each face, computes a plane from
    /// three vertices of the face's outer loop.
    fn attach_geometry(&mut self) -> Result<(), CoreError> {
        // Attach curves to edges.
        let edge_handles: Vec<_> = self.edges.keys().collect();
        for eh in edge_handles {
            let he_handle = self.edges[eh].half_edge;
            let he = &self.half_edges[he_handle];
            let start_pos = self.vertices[he.origin].position;
            let twin = &self.half_edges[he.twin];
            let end_pos = self.vertices[twin.origin].position;
            self.edges[eh].curve = Some(line_between(start_pos, end_pos));
        }

        // Attach surfaces to faces.
        let face_handles: Vec<_> = self.faces.keys().collect();
        for fh in face_handles {
            let outer_loop = self.faces[fh].outer_loop;
            let loop_data = self
                .loops
                .get(outer_loop)
                .ok_or_else(|| CoreError::InvalidHandle("face outer loop invalid".into()))?;
            let start_he = loop_data.half_edge;

            // Skip degenerate faces (no half-edges, e.g. the initial mvfs face).
            if !self.half_edges.contains_key(start_he) {
                continue;
            }

            // Collect 3 distinct vertex positions from the loop.
            let he0 = &self.half_edges[start_he];
            let p0 = self.vertices[he0.origin].position;

            let he1 = &self.half_edges[he0.next];
            let p1 = self.vertices[he1.origin].position;

            let he2 = &self.half_edges[he1.next];
            let p2 = self.vertices[he2.origin].position;

            let surface = plane_from_points(p0, p1, p2)?;
            self.faces[fh].surface = Some(SurfaceOnFace {
                geometry: surface,
                orientation: true,
            });
        }

        Ok(())
    }

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
        brep.extrude(mef.face, direction, distance)?;

        // Attach geometry to all edges and faces.
        brep.attach_geometry()?;

        // Validate the result.
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_cube_topology() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        assert_eq!(cube.vertices.len(), 8);
        assert_eq!(cube.edges.len(), 12);
        assert_eq!(cube.faces.len(), 6);
    }

    #[test]
    fn unit_cube_passes_validation() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        assert!(cube.validate().is_ok(), "Cube should pass full validation");
    }

    #[test]
    fn unit_cube_has_geometry() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        // Every edge should have a curve
        for edge in cube.edges.values() {
            assert!(edge.curve.is_some(), "Every edge should have geometry");
        }
        // Every face should have a surface
        for face in cube.faces.values() {
            assert!(face.surface.is_some(), "Every face should have geometry");
        }
    }

    #[test]
    fn prism_triangle_topology() {
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
        ];
        let prism = BRep::prism(
            &profile,
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            Tolerance::default(),
        )
        .unwrap();
        // Triangular prism: V=6, E=9, F=5
        assert_eq!(prism.vertices.len(), 6);
        assert_eq!(prism.edges.len(), 9);
        assert_eq!(prism.faces.len(), 5);
    }

    #[test]
    fn prism_passes_validation() {
        let profile = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
        ];
        let prism = BRep::prism(
            &profile,
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            Tolerance::default(),
        )
        .unwrap();
        assert!(prism.validate().is_ok(), "Prism should pass validation");
    }

    #[test]
    fn prism_too_few_vertices_fails() {
        let profile = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let result = BRep::prism(
            &profile,
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            Tolerance::default(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn prism_pentagon_topology() {
        let profile = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.309, 0.951, 0.0),
            Point3::new(-0.809, 0.588, 0.0),
            Point3::new(-0.809, -0.588, 0.0),
            Point3::new(0.309, -0.951, 0.0),
        ];
        let prism = BRep::prism(
            &profile,
            Vec3::new(0.0, 0.0, 1.0),
            2.0,
            Tolerance::default(),
        )
        .unwrap();
        // Pentagonal prism: V=10, E=15, F=7
        assert_eq!(prism.vertices.len(), 10);
        assert_eq!(prism.edges.len(), 15);
        assert_eq!(prism.faces.len(), 7);
    }
}
