//! BRep — the top-level boundary representation container.

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
        }
    }
}

use crate::geometry::transform::Transform3;

impl BRep {
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
        assert_eq!(brep.tolerance.linear, 1e-3);
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
