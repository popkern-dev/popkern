//! Geometry attachment methods for BRepBuilder.

use crate::constructors::builder::BRepBuilder;
use crate::geometry::curve::CurveOnEdge;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceOnFace;
use crate::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle};

impl BRepBuilder {
    /// Attach a curve to an edge.
    pub fn set_edge_curve(&mut self, edge: EdgeHandle, curve: CurveOnEdge) {
        self.brep.edges[edge].curve = Some(curve);
    }

    /// Attach a surface to a face.
    pub fn set_face_surface(&mut self, face: FaceHandle, surface: SurfaceOnFace) {
        self.brep.faces[face].surface = Some(surface);
    }

    /// Attach a pcurve (UV parameter curve) to a half-edge.
    pub fn set_half_edge_pcurve(&mut self, he: HalfEdgeHandle, pcurve: NurbsCurve) {
        self.brep.half_edges[he].pcurve = Some(pcurve);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;

    #[test]
    fn set_face_surface_attaches_surface() {
        let mut b = BRepBuilder::new(Tolerance::default());
        let solid = b.add_solid();
        let shell = b.add_shell(solid);
        let face = b.add_face(shell);
        assert!(b.brep.faces[face].surface.is_none());

        let surface = SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        };
        b.set_face_surface(face, surface);
        assert!(b.brep.faces[face].surface.is_some());
    }

    #[test]
    fn set_half_edge_pcurve_attaches_pcurve() {
        let mut b = BRepBuilder::new(Tolerance::default());
        let v1 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
        let v2 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
        let (_, he1, _) = b.add_edge(v1, v2);
        assert!(b.brep.half_edges[he1].pcurve.is_none());

        let pcurve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        b.set_half_edge_pcurve(he1, pcurve);
        assert!(b.brep.half_edges[he1].pcurve.is_some());
    }
}
