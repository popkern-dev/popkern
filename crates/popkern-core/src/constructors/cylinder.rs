//! Cylinder constructor.

use crate::error::CoreError;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

use super::plane_from_points;

impl BRep {
    /// Build a cylinder as a hexagonal prism with cylindrical surface geometry
    /// on the side faces and planar surfaces on the caps.
    ///
    /// The cylinder has its base centered at the origin in the XY plane and
    /// extends along the +Z axis by `height`. The hexagonal cross-section
    /// approximates a circle of the given `radius`.
    ///
    /// # Topology
    ///
    /// V=12, E=18, F=8 (6 side quads + 2 hexagonal caps).
    ///
    /// # Errors
    ///
    /// Returns [`CoreError`] if `radius` or `height` is not positive,
    /// or if any Euler operation or validation fails.
    pub fn cylinder(radius: f64, height: f64, tolerance: Tolerance) -> Result<Self, CoreError> {
        if radius <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cylinder radius must be positive".into(),
            ));
        }
        if height <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cylinder height must be positive".into(),
            ));
        }

        let n = 6;
        let profile: Vec<Point3> = (0..n)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                Point3::new(radius * angle.cos(), radius * angle.sin(), 0.0)
            })
            .collect();

        let mut brep =
            Self::prism_without_geometry(&profile, Vec3::new(0.0, 0.0, 1.0), height, tolerance)?;

        // Attach line curves to all edges.
        brep.attach_line_curves();

        // Attach planar surfaces to top and bottom cap faces, and cylindrical
        // surfaces to the 6 side faces.
        //
        // Strategy: for each face, check if all vertices share the same z
        // coordinate. If so, it's a cap face (planar). Otherwise it's a side
        // face (cylindrical).
        let face_handles: Vec<_> = brep.faces.keys().collect();
        for fh in face_handles {
            let outer_loop = brep.faces[fh].outer_loop;
            let loop_data = match brep.loops.get(outer_loop) {
                Some(l) => l,
                None => continue,
            };
            let start_he = loop_data.half_edge;
            if !brep.half_edges.contains_key(start_he) {
                continue;
            }

            // Collect vertex positions in this face's loop.
            let mut positions = Vec::new();
            let mut cursor = start_he;
            let bound = brep.half_edges.len();
            for _ in 0..=bound {
                let he = &brep.half_edges[cursor];
                positions.push(brep.vertices[he.origin].position);
                cursor = he.next;
                if cursor == start_he {
                    break;
                }
            }

            if positions.len() < 3 {
                continue;
            }

            // Check if all z-coordinates are the same (cap face).
            let z0 = positions[0].z;
            let is_cap = positions.iter().all(|p| (p.z - z0).abs() < 1e-10);

            let surface = if is_cap {
                plane_from_points(positions[0], positions[1], positions[2])?
            } else {
                SurfaceGeometry::Cylinder {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    axis: Vec3::new(0.0, 0.0, 1.0),
                    radius,
                }
            };

            brep.faces[fh].surface = Some(SurfaceOnFace {
                geometry: surface,
                orientation: true,
            });
        }

        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}
