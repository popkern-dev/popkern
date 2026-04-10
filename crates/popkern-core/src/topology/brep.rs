//! BRep — the top-level boundary representation container.

use crate::naming::{NamePool, OperationLabels};
use crate::tolerance::Tolerance;
use crate::topology::*;
use serde::{Deserialize, Serialize};
use slotmap::SlotMap;

/// A boundary representation model.
///
/// Owns all topological entities in typed arenas (SlotMaps).
/// All topology modifications go through Euler operators.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    /// All vertices.
    pub vertices: SlotMap<VertexHandle, Vertex>,
    /// All half-edges.
    pub half_edges: SlotMap<HalfEdgeHandle, HalfEdge>,
    /// All edges.
    pub edges: SlotMap<EdgeHandle, Edge>,
    /// All loops.
    pub loops: SlotMap<LoopHandle, Loop>,
    /// All faces.
    pub faces: SlotMap<FaceHandle, Face>,
    /// All shells.
    pub shells: SlotMap<ShellHandle, Shell>,
    /// All solids.
    pub solids: SlotMap<SolidHandle, Solid>,
    /// Geometric tolerance for this model.
    pub tolerance: Tolerance,
    /// Pool of stable topological names for entities in this model.
    #[serde(default)]
    pub name_pool: NamePool,
    /// Human-readable labels for modelling operations.
    #[serde(default)]
    pub operation_labels: OperationLabels,
}

impl BRep {
    /// Create a new empty BRep with the given tolerance.
    pub fn new(tolerance: Tolerance) -> Self {
        Self {
            vertices: SlotMap::with_key(),
            half_edges: SlotMap::with_key(),
            edges: SlotMap::with_key(),
            loops: SlotMap::with_key(),
            faces: SlotMap::with_key(),
            shells: SlotMap::with_key(),
            solids: SlotMap::with_key(),
            tolerance,
            name_pool: NamePool::new(),
            operation_labels: OperationLabels::default(),
        }
    }
}

use crate::geometry::transform::Transform3;

impl BRep {
    /// Collect vertex positions from the outer loop of a face by walking half-edges.
    ///
    /// Returns an empty `Vec` if the face handle is invalid, the loop is missing,
    /// or the loop has no half-edges.  This is the canonical implementation —
    /// downstream crates should call this instead of re-implementing the walk.
    pub fn face_outer_loop_positions(
        &self,
        face_handle: FaceHandle,
    ) -> Vec<crate::geometry::Point3> {
        let Some(face) = self.faces.get(face_handle) else {
            return Vec::new();
        };
        let Some(loop_) = self.loops.get(face.outer_loop) else {
            return Vec::new();
        };
        let start_he = loop_.half_edge;
        if !self.half_edges.contains_key(start_he) {
            return Vec::new();
        }
        let mut positions = Vec::new();
        let mut current = start_he;
        loop {
            let he = &self.half_edges[current];
            let vertex = &self.vertices[he.origin];
            positions.push(vertex.position);
            current = he.next;
            if current == start_he {
                break;
            }
        }
        positions
    }

    /// Collect sample points along the outer loop's edge curves.
    ///
    /// Walks the outer loop's half-edges. For each edge that has curve
    /// geometry, samples the curve at `n` equally-spaced parameter values
    /// and collects the 3D points. Includes vertex positions as well.
    ///
    /// This is more representative than [`face_outer_loop_positions`] for
    /// faces with self-loop edges (e.g., circular caps) where the vertex
    /// count is 1 but the curve traverses a full circle.
    pub fn face_boundary_sample_points(
        &self,
        face_handle: FaceHandle,
        samples_per_edge: usize,
    ) -> Vec<crate::geometry::Point3> {
        use crate::geometry::curve::CurveEval;

        let Some(face) = self.faces.get(face_handle) else {
            return Vec::new();
        };
        let Some(loop_) = self.loops.get(face.outer_loop) else {
            return Vec::new();
        };
        let start_he = loop_.half_edge;
        if !self.half_edges.contains_key(start_he) {
            return Vec::new();
        }

        let mut points = Vec::new();
        let mut current = start_he;
        loop {
            let he = &self.half_edges[current];
            // Always include the vertex position.
            points.push(self.vertices[he.origin].position);

            // If the edge has curve geometry, sample along it.
            let edge = &self.edges[he.edge];
            if let Some(ref curve) = edge.curve {
                let n = samples_per_edge.max(2);
                for i in 1..n {
                    let t = curve.t_start + (curve.t_end - curve.t_start) * (i as f64 / n as f64);
                    points.push(curve.geometry.point_at(t));
                }
            }

            current = he.next;
            if current == start_he {
                break;
            }
        }
        points
    }

    /// Get the effective tolerance for a vertex (entity-level or BRep default).
    pub fn vertex_tolerance(&self, vh: VertexHandle) -> f64 {
        self.vertices[vh]
            .tolerance
            .unwrap_or(self.tolerance.linear())
    }

    /// Get the effective tolerance for an edge (entity-level or BRep default).
    pub fn edge_tolerance(&self, eh: EdgeHandle) -> f64 {
        self.edges[eh].tolerance.unwrap_or(self.tolerance.linear())
    }

    /// Get the effective tolerance for a face (entity-level or BRep default).
    pub fn face_tolerance(&self, fh: FaceHandle) -> f64 {
        self.faces[fh].tolerance.unwrap_or(self.tolerance.linear())
    }

    /// Apply a rigid transformation to the entire BRep, returning a new copy.
    ///
    /// Transforms vertex positions, surface geometry on every face, edge
    /// curve geometry, and NURBS control points on pcurves (for non-planar
    /// faces; pcurves on planar faces live in UV space and are
    /// parameterization-invariant for rigid transforms that keep the plane
    /// canonical form, so they are regenerated via `attach_planar_pcurves`).
    pub fn transform(&self, t: &Transform3) -> Result<BRep, crate::error::CoreError> {
        let mut result = self.clone();

        // 1. Transform vertex positions.
        for (_, vertex) in &mut result.vertices {
            vertex.position = t.apply_point(&vertex.position);
        }

        // 2. Transform face surface geometry.
        for (_, face) in &mut result.faces {
            if let Some(ref mut sof) = face.surface {
                sof.geometry = transform_surface(&sof.geometry, t);
            }
        }

        // 3. Transform edge curve geometry.
        for (_, edge) in &mut result.edges {
            if let Some(ref mut curve_on_edge) = edge.curve {
                curve_on_edge.geometry = transform_curve(&curve_on_edge.geometry, t);
            }
        }

        // 4. Regenerate pcurves for planar faces (UV coordinates change
        //    when the plane origin/axes change).
        result.attach_planar_pcurves()?;

        Ok(result)
    }
}

/// Transform a `SurfaceGeometry` by a rigid transform.
fn transform_surface(
    surface: &crate::geometry::surface::SurfaceGeometry,
    t: &Transform3,
) -> crate::geometry::surface::SurfaceGeometry {
    use crate::geometry::surface::SurfaceGeometry;

    match surface {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => SurfaceGeometry::Plane {
            origin: t.apply_point(origin),
            normal: t.apply_direction(normal),
            u_axis: t.apply_direction(u_axis),
            v_axis: t.apply_direction(v_axis),
        },
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => SurfaceGeometry::Cylinder {
            origin: t.apply_point(origin),
            axis: t.apply_direction(axis),
            radius: *radius,
        },
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => SurfaceGeometry::Cone {
            apex: t.apply_point(apex),
            axis: t.apply_direction(axis),
            half_angle: *half_angle,
        },
        SurfaceGeometry::Sphere { center, radius } => SurfaceGeometry::Sphere {
            center: t.apply_point(center),
            radius: *radius,
        },
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => SurfaceGeometry::Torus {
            center: t.apply_point(center),
            axis: t.apply_direction(axis),
            major_radius: *major_radius,
            minor_radius: *minor_radius,
        },
        SurfaceGeometry::Nurbs(nurbs) => {
            let control_points = nurbs
                .control_points
                .iter()
                .map(|p| t.apply_point(p))
                .collect();
            SurfaceGeometry::Nurbs(crate::geometry::NurbsSurface {
                degree_u: nurbs.degree_u,
                degree_v: nurbs.degree_v,
                knots_u: nurbs.knots_u.clone(),
                knots_v: nurbs.knots_v.clone(),
                control_points,
                num_u: nurbs.num_u,
                num_v: nurbs.num_v,
                weights: nurbs.weights.clone(),
            })
        }
    }
}

/// Transform a `CurveGeometry` by a rigid transform.
fn transform_curve(
    curve: &crate::geometry::curve::CurveGeometry,
    t: &Transform3,
) -> crate::geometry::curve::CurveGeometry {
    use crate::geometry::curve::CurveGeometry;

    match curve {
        CurveGeometry::Line { origin, direction } => CurveGeometry::Line {
            origin: t.apply_point(origin),
            direction: t.apply_direction(direction),
        },
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => CurveGeometry::Circle {
            center: t.apply_point(center),
            normal: t.apply_direction(normal),
            radius: *radius,
        },
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        } => CurveGeometry::Arc {
            center: t.apply_point(center),
            normal: t.apply_direction(normal),
            radius: *radius,
            start_angle: *start_angle,
            end_angle: *end_angle,
        },
        CurveGeometry::Nurbs(nurbs) => {
            let control_points = nurbs
                .control_points
                .iter()
                .map(|p| t.apply_point(p))
                .collect();
            CurveGeometry::Nurbs(crate::geometry::NurbsCurve {
                degree: nurbs.degree,
                knots: nurbs.knots.clone(),
                control_points,
                weights: nurbs.weights.clone(),
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_brep_is_empty() {
        let brep = BRep::new(Tolerance::default());
        assert_eq!(brep.vertices.len(), 0);
        assert_eq!(brep.half_edges.len(), 0);
        assert_eq!(brep.edges.len(), 0);
        assert_eq!(brep.loops.len(), 0);
        assert_eq!(brep.faces.len(), 0);
        assert_eq!(brep.shells.len(), 0);
        assert_eq!(brep.solids.len(), 0);
    }

    #[test]
    fn brep_has_tolerance() {
        let tol = Tolerance::default().with_linear(1e-3);
        let brep = BRep::new(tol);
        assert_eq!(brep.tolerance.linear(), 1e-3);
    }

    #[test]
    fn brep_serializes_round_trip() {
        let brep = BRep::new(Tolerance::default());
        let json = serde_json::to_string(&brep).unwrap();
        let brep2: BRep = serde_json::from_str(&json).unwrap();
        assert_eq!(brep2.vertices.len(), 0);
    }

    #[test]
    fn brep_transform_translates_vertices() {
        use crate::geometry::transform::Transform3;
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let t = Transform3::translation(10.0, 0.0, 0.0);
        let transformed = brep.transform(&t).unwrap();
        for (_, v) in &transformed.vertices {
            assert!(v.position.x >= 10.0 - 1e-10);
        }
    }
}

#[cfg(test)]
mod face_loop_tests {
    use super::*;
    use crate::tolerance::Tolerance;

    #[test]
    fn face_outer_loop_positions_unit_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face_handle = brep.faces.keys().next().unwrap();
        let positions = brep.face_outer_loop_positions(face_handle);
        assert_eq!(positions.len(), 4, "each cube face should have 4 vertices");
        for pos in &positions {
            assert!(
                pos.x >= -1e-10 && pos.x <= 1.0 + 1e-10,
                "x={} out of [0,1]",
                pos.x
            );
            assert!(
                pos.y >= -1e-10 && pos.y <= 1.0 + 1e-10,
                "y={} out of [0,1]",
                pos.y
            );
            assert!(
                pos.z >= -1e-10 && pos.z <= 1.0 + 1e-10,
                "z={} out of [0,1]",
                pos.z
            );
        }
    }

    #[test]
    fn face_outer_loop_positions_invalid_face_returns_empty() {
        let brep = BRep::new(Tolerance::default());
        let bogus = FaceHandle::default();
        let positions = brep.face_outer_loop_positions(bogus);
        assert!(
            positions.is_empty(),
            "invalid handle should yield empty vec"
        );
    }
}
