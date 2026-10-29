//! The edges blended in one fillet operation, with each one's support
//! rails — what the anchor-split stage needs to tell a fellow fillet
//! edge from a sharp neighbour and, at a reflex corner shared by two
//! fillet edges, where their support rails meet (PK-740).

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{EdgeHandle, FaceHandle};

use super::EdgeGeometry;

/// Fillet edges of the operation and their rails per support face.
#[derive(Debug, Default)]
pub(crate) struct BlendEdgeSet {
    edges: HashSet<EdgeHandle>,
    rails: HashMap<(EdgeHandle, FaceHandle), NurbsCurve>,
    faces: HashMap<EdgeHandle, (FaceHandle, FaceHandle)>,
    concave: HashSet<EdgeHandle>,
    /// Planar strips (chamfer bevels) as `(origin, normal)`.
    bevels: HashMap<EdgeHandle, (Point3, Vec3)>,
}

impl BlendEdgeSet {
    /// Every analysed edge with its two support rails.
    pub(crate) fn from_geometries(edge_geometries: &[EdgeGeometry]) -> Self {
        let mut set = Self::default();
        for eg in edge_geometries {
            let e = eg.analysis.edge;
            set.edges.insert(e);
            set.rails
                .insert((e, eg.analysis.face_1), eg.trim_3d_1.clone());
            set.rails
                .insert((e, eg.analysis.face_2), eg.trim_3d_2.clone());
            set.faces
                .insert(e, (eg.analysis.face_1, eg.analysis.face_2));
            if eg.analysis.convexity == crate::analysis::Convexity::Concave {
                set.concave.insert(e);
            }
            if let SurfaceGeometry::Plane { origin, normal, .. } = &eg.surface {
                set.bevels.insert(e, (*origin, *normal));
            }
        }
        set
    }

    /// Whether `edge` is blended by this operation.
    pub(crate) fn contains(&self, edge: &EdgeHandle) -> bool {
        self.edges.contains(edge)
    }

    /// The rail `edge`'s strip lays on `face`, when `edge` is blended.
    pub(crate) fn rail(&self, edge: EdgeHandle, face: FaceHandle) -> Option<&NurbsCurve> {
        self.rails.get(&(edge, face))
    }

    /// The neighbour's rail on `face` when `edge` and `neighbour` are
    /// both CONCAVE fillets — a boss-base reflex pair (PK-740). Saddle
    /// corners (a convex fillet meeting a concave one) belong to the
    /// setback machinery and never mitre here.
    pub(crate) fn concave_pair_rail(
        &self,
        edge: EdgeHandle,
        neighbour: EdgeHandle,
        face: FaceHandle,
    ) -> Option<&NurbsCurve> {
        (self.concave.contains(&edge) && self.concave.contains(&neighbour))
            .then(|| self.rail(neighbour, face))
            .flatten()
    }

    /// The bevel plane of `edge` when it is blended with a planar
    /// strip (a chamfer).
    pub(crate) fn bevel(&self, edge: EdgeHandle) -> Option<(Point3, Vec3)> {
        self.bevels.get(&edge).copied()
    }

    /// Whether every edge at `vertex` of `brep` is blended by this
    /// operation (a corner chamfered all round belongs to the mitre
    /// apex stage, not to the bevel-pair mitre).
    pub(crate) fn blends_every_edge_at(
        &self,
        brep: &popkern_core::topology::BRep,
        vertex: popkern_core::topology::VertexHandle,
    ) -> bool {
        brep.half_edges
            .values()
            .filter(|he| he.origin == vertex)
            .all(|he| self.edges.contains(&he.edge))
    }

    /// The support `edge`'s strip lays on other than `face`.
    pub(crate) fn other_face(&self, edge: EdgeHandle, face: FaceHandle) -> Option<FaceHandle> {
        let &(a, b) = self.faces.get(&edge)?;
        if a == face {
            Some(b)
        } else if b == face {
            Some(a)
        } else {
            None
        }
    }
}
