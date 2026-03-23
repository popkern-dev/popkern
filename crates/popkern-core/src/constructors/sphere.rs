//! Sphere constructor.

use std::f64::consts::{PI, TAU};

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;
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
    /// center and radius. Edges have line curves (seam lines from pole to pole).
    /// Each half-edge carries a pcurve (degree-1 NURBS line in UV space) representing
    /// a meridian at constant u, running from v=PI (south pole) to v=0 (north pole)
    /// or vice versa.
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
        let mef1 = brep.mef(mev1.half_edge, mev1.half_edge_twin)?;

        // 4. mef to create edge 3 and face 3
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
        for fh in &face_handles {
            brep.faces[*fh].surface = Some(SurfaceOnFace {
                geometry: sphere_geom.clone(),
                orientation: true,
            });
        }

        // Attach pcurves to all half-edges.
        // Each face is a lune bounded by two meridians. For face index i,
        // the left meridian is at u = i * TAU/3, the right at u = (i+1) * TAU/3.
        // In UV space: south pole = (u, PI), north pole = (u, 0).
        let south_vertex = mvfs.vertex;
        for (face_idx, &fh) in face_handles.iter().enumerate() {
            let u_left = face_idx as f64 * TAU / 3.0;
            let u_right = (face_idx as f64 + 1.0) * TAU / 3.0;

            let outer_loop = brep.faces[fh].outer_loop;
            let loop_data = brep
                .loops
                .get(outer_loop)
                .ok_or_else(|| CoreError::InvalidHandle("face outer loop invalid".into()))?;
            let start_he = loop_data.half_edge;

            if !brep.half_edges.contains_key(start_he) {
                continue;
            }

            // Walk the loop (exactly 2 half-edges per lune face).
            let mut he_handles = vec![start_he];
            let mut current = brep.half_edges[start_he].next;
            while current != start_he {
                he_handles.push(current);
                current = brep.half_edges[current].next;
            }

            for &heh in &he_handles {
                let origin = brep.half_edges[heh].origin;
                let is_from_south = origin == south_vertex;

                let pcurve = if is_from_south {
                    // South->North: (u_left, PI) -> (u_left, 0)
                    NurbsCurve {
                        degree: 1,
                        knots: vec![0.0, 0.0, 1.0, 1.0],
                        control_points: vec![
                            Point3::new(u_left, PI, 0.0),
                            Point3::new(u_left, 0.0, 0.0),
                        ],
                        weights: vec![1.0, 1.0],
                    }
                } else {
                    // North->South: (u_right, 0) -> (u_right, PI)
                    NurbsCurve {
                        degree: 1,
                        knots: vec![0.0, 0.0, 1.0, 1.0],
                        control_points: vec![
                            Point3::new(u_right, 0.0, 0.0),
                            Point3::new(u_right, PI, 0.0),
                        ],
                        weights: vec![1.0, 1.0],
                    }
                };

                brep.half_edges[heh].pcurve = Some(pcurve);
            }
        }

        // Validate the result.
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}
