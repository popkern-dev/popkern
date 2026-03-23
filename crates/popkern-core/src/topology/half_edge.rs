use crate::topology::{EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};
use serde::{Deserialize, Serialize};

/// A directed half-edge in the BRep model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HalfEdge {
    /// Vertex at the start of this half-edge.
    pub origin: VertexHandle,
    /// The opposite-direction half-edge.
    pub twin: HalfEdgeHandle,
    /// Next half-edge in the face loop.
    pub next: HalfEdgeHandle,
    /// Previous half-edge in the face loop.
    pub prev: HalfEdgeHandle,
    /// Parent edge.
    pub edge: EdgeHandle,
    /// Parent loop.
    pub loop_: LoopHandle,
    /// UV parameter curve on the parent face's surface.
    /// Uses Point3(u, v, 0) convention. `None` for edges without surface geometry.
    pub pcurve: Option<crate::geometry::nurbs_curve::NurbsCurve>,
}

#[cfg(test)]
mod tests {
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn half_edge_pcurve_is_none_by_default() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        for (_, he) in &brep.half_edges {
            assert!(he.pcurve.is_none());
        }
    }
}
