//! The near region of a certification pair (PK-900).
//!
//! A pair the base meshes cannot decide is re-meshed finer, but only
//! where it can matter: a triangle whose certified lower bound to the
//! other body (mesh distance less both triangles' deviations) exceeds the
//! band is farther than the band on the true faces too, so contact or
//! overlap can only involve the faces of the triangles inside it. A
//! screw in a 0.1 mm thread fit used to re-mesh both whole bodies, the
//! casting's fins and cycloid flanks with it, at every finer level.

use std::collections::BTreeSet;

use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::{FaceMesh, MeshError, TessellationOptions};

use super::{Aabb, LocalMesh, WorldMesh, build, deviation, triangle_distance};

impl WorldMesh {
    /// Local triangle indices of this mesh and of `other` whose
    /// certified lower bound to the other mesh is within `cutoff`.
    pub fn near_triangles(&self, other: &WorldMesh, cutoff: f64) -> (Vec<usize>, Vec<usize>) {
        let (mut here, mut there) = (BTreeSet::new(), BTreeSet::new());
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return (Vec::new(), Vec::new());
        }
        let mut stack = vec![(0usize, 0usize)];
        while let Some((i, j)) = stack.pop() {
            let (a, b) = (&self.nodes[i], &other.nodes[j]);
            if a.bounds.distance(&b.bounds) - a.eps - b.eps > cutoff {
                continue;
            }
            match (a.children, b.children) {
                (None, None) => {
                    for ta in a.first..a.first + a.count {
                        for tb in b.first..b.first + b.count {
                            let (d, _, _) = triangle_distance(&self.tris[ta], &other.tris[tb]);
                            if d - self.eps[ta] - other.eps[tb] <= cutoff {
                                here.insert(self.local[ta]);
                                there.insert(other.local[tb]);
                            }
                        }
                    }
                }
                (Some((l, r)), None) => stack.extend([(l, j), (r, j)]),
                (None, Some((l, r))) => stack.extend([(i, l), (i, r)]),
                (Some((al, ar)), Some((bl, br))) => {
                    // Split the larger box.
                    if a.bounds.diagonal() >= b.bounds.diagonal() {
                        stack.extend([(al, j), (ar, j)]);
                    } else {
                        stack.extend([(i, bl), (i, br)]);
                    }
                }
            }
        }
        (here.into_iter().collect(), there.into_iter().collect())
    }
}

impl LocalMesh {
    /// The fine meshes of `faces` alone, as one local mesh (a partial
    /// body: not closed, so containment is not asked of it).
    pub fn from_faces(
        brep: &BRep,
        faces: &[FaceHandle],
        options: &TessellationOptions,
    ) -> Result<LocalMesh, MeshError> {
        let meshes: Vec<FaceMesh> = popkern_mesh::tessellate_faces(brep, faces, options)?;
        let mut out = LocalMesh::empty();
        for fm in &meshes {
            for t in &fm.triangles {
                let corners = t.map(|i| fm.vertices[i].position);
                let uvs = t.map(|i| fm.vertices[i].uv);
                let eps = deviation(brep, fm.face, &corners, &uvs);
                out.max_eps = out.max_eps.max(eps);
                for p in &corners {
                    out.bounds.add(p);
                }
                out.eps.push(eps);
                out.tris.push(corners);
                out.faces.push(fm.face);
                out.uvs.push(uvs);
            }
        }
        out.index();
        Ok(out)
    }

    /// The meshes of `pieces` (each a face set of one body) as one.
    pub fn concat(pieces: &[std::sync::Arc<LocalMesh>]) -> LocalMesh {
        let mut out = LocalMesh::empty();
        for piece in pieces {
            out.tris.extend_from_slice(&piece.tris);
            out.faces.extend_from_slice(&piece.faces);
            out.eps.extend_from_slice(&piece.eps);
            out.uvs.extend_from_slice(&piece.uvs);
            for k in 0..3 {
                out.bounds.lo[k] = out.bounds.lo[k].min(piece.bounds.lo[k]);
                out.bounds.hi[k] = out.bounds.hi[k].max(piece.bounds.hi[k]);
            }
            out.max_eps = out.max_eps.max(piece.max_eps);
        }
        out.index();
        out
    }

    fn empty() -> LocalMesh {
        LocalMesh {
            tris: Vec::new(),
            faces: Vec::new(),
            eps: Vec::new(),
            uvs: Vec::new(),
            bounds: Aabb::empty(),
            max_eps: 0.0,
            order: Vec::new(),
            nodes: Vec::new(),
        }
    }

    /// Builds the hierarchy (none for an empty region: a body with no
    /// face near the other stays in the pair with no triangles).
    fn index(&mut self) {
        self.order = (0..self.tris.len()).collect();
        self.nodes.clear();
        if !self.tris.is_empty() {
            build(
                &self.tris,
                &self.eps,
                &mut self.order,
                0,
                self.tris.len(),
                &mut self.nodes,
            );
        }
    }
}
