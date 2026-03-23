//! Cylinder constructor.

use crate::error::CoreError;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

use super::plane_from_points;

/// Map a 3D point on a Z-axis cylinder to UV parameters.
///
/// For a cylinder with `origin=(0,0,0)`, `axis=(0,0,1)`, `radius=r`:
/// `u = atan2(y, x)` shifted to `[0, 2PI)`, `v = z`.
fn cylinder_uv(pos: Point3) -> (f64, f64) {
    let mut u = pos.y.atan2(pos.x);
    if u < -1e-14 {
        u += std::f64::consts::TAU;
    }
    (u, pos.z)
}

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
    /// # Geometry
    ///
    /// Side faces carry `SurfaceGeometry::Cylinder`; cap faces carry
    /// `SurfaceGeometry::Plane`. Every half-edge has a pcurve in the
    /// parent face's UV space.
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
        for fh in &face_handles {
            let outer_loop = brep.faces[*fh].outer_loop;
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

            brep.faces[*fh].surface = Some(SurfaceOnFace {
                geometry: surface,
                orientation: true,
            });
        }

        // Attach pcurves: planar caps via the shared helper, cylindrical sides manually.
        brep.attach_planar_pcurves()?;
        brep.attach_cylinder_side_pcurves()?;

        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }

    /// Attach pcurves to half-edges on faces with cylindrical surface geometry.
    ///
    /// For each cylindrical side face, projects vertex positions to the
    /// cylinder's UV space using `u = atan2(y, x)` and `v = z`.
    fn attach_cylinder_side_pcurves(&mut self) -> Result<(), CoreError> {
        let face_handles: Vec<_> = self.faces.keys().collect();
        for fh in face_handles {
            let is_cylinder = matches!(
                self.faces[fh].surface,
                Some(SurfaceOnFace {
                    geometry: SurfaceGeometry::Cylinder { .. },
                    ..
                })
            );
            if !is_cylinder {
                continue;
            }

            let outer_loop = self.faces[fh].outer_loop;
            let loop_data = self
                .loops
                .get(outer_loop)
                .ok_or_else(|| CoreError::InvalidHandle("face outer loop invalid".into()))?;
            let start_he = loop_data.half_edge;

            if !self.half_edges.contains_key(start_he) {
                continue;
            }

            // Collect half-edge handles for this loop.
            let mut he_handles = vec![start_he];
            let mut current = self.half_edges[start_he].next;
            while current != start_he {
                he_handles.push(current);
                current = self.half_edges[current].next;
            }

            // First pass: compute all UV coordinates for consistency.
            // We need to ensure u values are consistent within a face
            // (no jumps across the 0/2PI seam).
            let mut uvs: Vec<(f64, f64)> = Vec::new();
            for &heh in &he_handles {
                let pos = self.vertices[self.half_edges[heh].origin].position;
                uvs.push(cylinder_uv(pos));
            }

            // Ensure u values don't jump across seam within a single face.
            // All u values should be close to each other. If any pair differs
            // by more than PI, adjust by adding/subtracting TAU.
            if !uvs.is_empty() {
                let u_ref = uvs[0].0;
                for uv in &mut uvs[1..] {
                    while uv.0 - u_ref > std::f64::consts::PI {
                        uv.0 -= std::f64::consts::TAU;
                    }
                    while u_ref - uv.0 > std::f64::consts::PI {
                        uv.0 += std::f64::consts::TAU;
                    }
                }
            }

            // Assign pcurves as degree-1 NURBS lines in UV space.
            let n_he = he_handles.len();
            for (i, &heh) in he_handles.iter().enumerate() {
                let (u_start, v_start) = uvs[i];
                let (u_end, v_end) = uvs[(i + 1) % n_he];

                let pcurve = NurbsCurve {
                    degree: 1,
                    knots: vec![0.0, 0.0, 1.0, 1.0],
                    control_points: vec![
                        Point3::new(u_start, v_start, 0.0),
                        Point3::new(u_end, v_end, 0.0),
                    ],
                    weights: vec![1.0, 1.0],
                };

                self.half_edges[heh].pcurve = Some(pcurve);
            }
        }
        Ok(())
    }
}
