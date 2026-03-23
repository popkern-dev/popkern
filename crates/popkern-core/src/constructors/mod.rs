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
use crate::geometry::Point3;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
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

    /// Attach pcurves to all half-edges on faces with planar surfaces.
    ///
    /// For each face with a `Plane` surface, walks the outer loop and computes
    /// a degree-1 NURBS line in UV space for each half-edge. The UV projection
    /// uses the inverse of `point_at(u,v) = origin + u_axis*u + v_axis*v`,
    /// which for orthogonal axes is:
    /// `u = (p - origin) . u_axis / |u_axis|^2`
    /// `v = (p - origin) . v_axis / |v_axis|^2`
    pub(super) fn attach_planar_pcurves(&mut self) -> Result<(), CoreError> {
        let face_handles: Vec<_> = self.faces.keys().collect();
        for fh in face_handles {
            let surface_data = match &self.faces[fh].surface {
                Some(sof) => match &sof.geometry {
                    SurfaceGeometry::Plane {
                        origin,
                        u_axis,
                        v_axis,
                        ..
                    } => Some((*origin, *u_axis, *v_axis)),
                    _ => None,
                },
                None => None,
            };

            let (plane_origin, u_axis, v_axis) = match surface_data {
                Some(data) => data,
                None => continue,
            };

            let u_dot_u = u_axis.dot(&u_axis);
            let v_dot_v = v_axis.dot(&v_axis);

            // Skip degenerate axes.
            if u_dot_u < 1e-30 || v_dot_v < 1e-30 {
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

            // Compute and assign pcurves.
            for &heh in &he_handles {
                let origin_pos = self.vertices[self.half_edges[heh].origin].position;
                let next_heh = self.half_edges[heh].next;
                let dest_pos = self.vertices[self.half_edges[next_heh].origin].position;

                let d_start = origin_pos - plane_origin;
                let u_start = d_start.dot(&u_axis) / u_dot_u;
                let v_start = d_start.dot(&v_axis) / v_dot_v;

                let d_end = dest_pos - plane_origin;
                let u_end = d_end.dot(&u_axis) / u_dot_u;
                let v_end = d_end.dot(&v_axis) / v_dot_v;

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

#[cfg(test)]
mod tests;
