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
