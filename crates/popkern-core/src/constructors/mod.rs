//! Convenience constructors for common solids.
//!
//! Build topology via Euler operators, then attach analytic geometry
//! (line curves and planar surfaces) to all edges and faces.

mod cone;
mod cube;
mod cylinder;
mod prism;
mod sphere;
mod torus;

use crate::error::CoreError;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::topology::BRep;

/// Compute a plane surface from 3 non-collinear points.
pub(super) fn plane_from_points(
    p0: Point3,
    p1: Point3,
    p2: Point3,
) -> Result<SurfaceGeometry, CoreError> {
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
pub(super) fn line_between(start: Point3, end: Point3) -> CurveOnEdge {
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
    /// Attach line curves to all edges that lack geometry.
    pub(super) fn attach_line_curves(&mut self) {
        let edge_handles: Vec<_> = self.edges.keys().collect();
        for eh in edge_handles {
            if self.edges[eh].curve.is_some() {
                continue;
            }
            let he_handle = self.edges[eh].half_edge;
            let he = &self.half_edges[he_handle];
            let start_pos = self.vertices[he.origin].position;
            let twin = &self.half_edges[he.twin];
            let end_pos = self.vertices[twin.origin].position;
            self.edges[eh].curve = Some(line_between(start_pos, end_pos));
        }
    }

    /// Attach line curves to all edges and planar surfaces to all faces.
    ///
    /// For each edge, computes the line from the origin of one half-edge
    /// to the origin of its twin. For each face, computes a plane from
    /// three vertices of the face's outer loop.
    pub(super) fn attach_geometry(&mut self) -> Result<(), CoreError> {
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

    /// Build a prism without attaching geometry (internal helper for cylinder).
    pub(super) fn prism_without_geometry(
        profile: &[Point3],
        direction: Vec3,
        distance: f64,
        tolerance: crate::tolerance::Tolerance,
    ) -> Result<Self, CoreError> {
        if profile.len() < 3 {
            return Err(CoreError::PreconditionFailed(
                "prism profile must have at least 3 vertices".into(),
            ));
        }

        let mut brep = BRep::new(tolerance);

        let mvfs = brep.mvfs(profile[0])?;
        let mev1 = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, profile[1])?;

        let mut last_twin = mev1.half_edge_twin;
        for &point in &profile[2..] {
            let mev = brep.mev(last_twin, point)?;
            last_twin = mev.half_edge_twin;
        }

        let mef = brep.mef(last_twin, mev1.half_edge)?;
        brep.extrude(mef.face, direction, distance)?;

        Ok(brep)
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

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

    // --- Cylinder tests ---

    #[test]
    fn cylinder_topology() {
        let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
        // Hexagonal prism: V=12, E=18, F=8
        assert_eq!(cyl.vertices.len(), 12);
        assert_eq!(cyl.edges.len(), 18);
        assert_eq!(cyl.faces.len(), 8);
    }

    #[test]
    fn cylinder_passes_validation() {
        let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
        assert!(cyl.validate().is_ok(), "Cylinder should pass validation");
    }

    #[test]
    fn cylinder_has_geometry() {
        let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
        for edge in cyl.edges.values() {
            assert!(
                edge.curve.is_some(),
                "Every edge should have curve geometry"
            );
        }
        for face in cyl.faces.values() {
            assert!(
                face.surface.is_some(),
                "Every face should have surface geometry"
            );
        }
    }

    #[test]
    fn cylinder_has_cylindrical_side_faces() {
        let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
        let mut cylinder_faces = 0;
        let mut plane_faces = 0;
        for face in cyl.faces.values() {
            if let Some(ref sof) = face.surface {
                match sof.geometry {
                    SurfaceGeometry::Cylinder { .. } => cylinder_faces += 1,
                    SurfaceGeometry::Plane { .. } => plane_faces += 1,
                    _ => {}
                }
            }
        }
        assert_eq!(cylinder_faces, 6, "Should have 6 cylindrical side faces");
        assert_eq!(plane_faces, 2, "Should have 2 planar cap faces");
    }

    #[test]
    fn cylinder_invalid_radius() {
        let result = BRep::cylinder(0.0, 1.0, Tolerance::default());
        assert!(result.is_err());
        let result = BRep::cylinder(-1.0, 1.0, Tolerance::default());
        assert!(result.is_err());
    }

    #[test]
    fn cylinder_invalid_height() {
        let result = BRep::cylinder(1.0, 0.0, Tolerance::default());
        assert!(result.is_err());
        let result = BRep::cylinder(1.0, -1.0, Tolerance::default());
        assert!(result.is_err());
    }

    // --- Sphere tests ---

    #[test]
    fn sphere_topology() {
        let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
        assert_eq!(sph.vertices.len(), 2);
        assert_eq!(sph.edges.len(), 3);
        assert_eq!(sph.faces.len(), 3);
    }

    #[test]
    fn sphere_passes_validation() {
        let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
        assert!(sph.validate().is_ok(), "Sphere should pass validation");
    }

    #[test]
    fn sphere_has_edge_geometry() {
        let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
        for edge in sph.edges.values() {
            assert!(
                edge.curve.is_some(),
                "Every edge should have curve geometry"
            );
        }
    }

    #[test]
    fn sphere_vertex_positions() {
        let sph = BRep::sphere(2.5, Tolerance::default()).unwrap();
        let positions: Vec<_> = sph.vertices.values().map(|v| v.position).collect();
        // One vertex at (0,0,-2.5) and one at (0,0,2.5)
        let has_south = positions
            .iter()
            .any(|p| (p.z + 2.5).abs() < 1e-10 && p.x.abs() < 1e-10 && p.y.abs() < 1e-10);
        let has_north = positions
            .iter()
            .any(|p| (p.z - 2.5).abs() < 1e-10 && p.x.abs() < 1e-10 && p.y.abs() < 1e-10);
        assert!(has_south, "Should have south pole vertex");
        assert!(has_north, "Should have north pole vertex");
    }

    #[test]
    fn sphere_invalid_radius() {
        let result = BRep::sphere(0.0, Tolerance::default());
        assert!(result.is_err());
        let result = BRep::sphere(-1.0, Tolerance::default());
        assert!(result.is_err());
    }

    // --- Cone tests ---

    #[test]
    fn cone_topology() {
        let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
        assert_eq!(cone.vertices.len(), 2);
        assert_eq!(cone.edges.len(), 3);
        assert_eq!(cone.faces.len(), 3);
    }

    #[test]
    fn cone_passes_validation() {
        let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
        assert!(cone.validate().is_ok(), "Cone should pass validation");
    }

    #[test]
    fn cone_has_edge_geometry() {
        let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
        for edge in cone.edges.values() {
            assert!(
                edge.curve.is_some(),
                "Every edge should have curve geometry"
            );
        }
    }

    #[test]
    fn cone_vertex_positions() {
        let cone = BRep::cone(1.0, 3.0, Tolerance::default()).unwrap();
        let positions: Vec<_> = cone.vertices.values().map(|v| v.position).collect();
        let has_base = positions
            .iter()
            .any(|p| p.z.abs() < 1e-10 && p.x.abs() < 1e-10 && p.y.abs() < 1e-10);
        let has_apex = positions
            .iter()
            .any(|p| (p.z - 3.0).abs() < 1e-10 && p.x.abs() < 1e-10 && p.y.abs() < 1e-10);
        assert!(has_base, "Should have base vertex");
        assert!(has_apex, "Should have apex vertex");
    }

    #[test]
    fn cone_invalid_params() {
        assert!(BRep::cone(0.0, 1.0, Tolerance::default()).is_err());
        assert!(BRep::cone(-1.0, 1.0, Tolerance::default()).is_err());
        assert!(BRep::cone(1.0, 0.0, Tolerance::default()).is_err());
        assert!(BRep::cone(1.0, -1.0, Tolerance::default()).is_err());
    }

    // --- Torus tests ---

    #[test]
    fn torus_basic_succeeds() {
        let result = BRep::torus(2.0, 0.5, Tolerance::default());
        assert!(result.is_ok(), "Torus constructor should succeed");
    }
}
