//! Direct BRep construction for import scenarios.
//!
//! Bypasses Euler operators by directly populating the BRep slotmaps.
//! Validates the finished topology on `.build()`.

use crate::error::CoreError;
use crate::geometry::point::Point3;
use crate::tolerance::Tolerance;
use crate::topology::*;

/// Builder for constructing BRep topology directly (without Euler operators).
///
/// Used by importers (STEP, IGES) that have finished topology descriptions
/// and need to populate the BRep slotmaps in bulk. Call `.build()` when done
/// to validate and consume the builder.
pub struct BRepBuilder {
    pub(crate) brep: BRep,
}

impl BRepBuilder {
    /// Create a new empty builder with the given tolerance.
    pub fn new(tolerance: Tolerance) -> Self {
        Self {
            brep: BRep::new(tolerance),
        }
    }

    /// Insert a vertex at the given position.
    pub fn add_vertex(&mut self, position: Point3) -> VertexHandle {
        self.brep.vertices.insert(Vertex {
            position,
            outgoing: HalfEdgeHandle::default(),
            name: None,
        })
    }

    /// Insert an empty solid and return its handle.
    pub fn add_solid(&mut self) -> SolidHandle {
        self.brep.solids.insert(Solid {
            shell: ShellHandle::default(),
        })
    }

    /// Insert a shell belonging to the given solid.
    pub fn add_shell(&mut self, solid: SolidHandle) -> ShellHandle {
        let sh = self.brep.shells.insert(Shell {
            face: FaceHandle::default(),
            solid,
        });
        self.brep.solids[solid].shell = sh;
        sh
    }

    /// Insert an edge between two vertices, creating a twin half-edge pair.
    ///
    /// Returns `(EdgeHandle, HalfEdgeHandle, HalfEdgeHandle)` where the first
    /// half-edge goes from `v1` to `v2` and the second from `v2` to `v1`.
    ///
    /// The half-edges' `next`, `prev`, and `loop_` fields are left as defaults.
    /// Wire them via `add_loop`.
    pub fn add_edge(
        &mut self,
        v1: VertexHandle,
        v2: VertexHandle,
    ) -> (EdgeHandle, HalfEdgeHandle, HalfEdgeHandle) {
        let edge_handle = self.brep.edges.insert(Edge {
            half_edge: HalfEdgeHandle::default(),
            curve: None,
            name: None,
        });

        let he1 = self.brep.half_edges.insert(HalfEdge {
            origin: v1,
            twin: HalfEdgeHandle::default(),
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge: edge_handle,
            loop_: LoopHandle::default(),
            pcurve: None,
        });

        let he2 = self.brep.half_edges.insert(HalfEdge {
            origin: v2,
            twin: he1,
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge: edge_handle,
            loop_: LoopHandle::default(),
            pcurve: None,
        });

        self.brep.half_edges[he1].twin = he2;
        self.brep.edges[edge_handle].half_edge = he1;
        if self.brep.vertices[v1].outgoing == HalfEdgeHandle::default() {
            self.brep.vertices[v1].outgoing = he1;
        }
        if self.brep.vertices[v2].outgoing == HalfEdgeHandle::default() {
            self.brep.vertices[v2].outgoing = he2;
        }

        (edge_handle, he1, he2)
    }

    /// Insert a face belonging to the given shell.
    ///
    /// The face has no loops yet. Add them with `add_loop` and `add_inner_loop`.
    pub fn add_face(&mut self, shell: ShellHandle) -> FaceHandle {
        let fh = self.brep.faces.insert(Face {
            outer_loop: LoopHandle::default(),
            inner_loops: Vec::new(),
            shell,
            surface: None,
            name: None,
        });
        if self.brep.shells[shell].face == FaceHandle::default() {
            self.brep.shells[shell].face = fh;
        }
        fh
    }

    /// Wire a sequence of half-edges into a loop and assign it as
    /// the outer loop of the given face.
    ///
    /// Sets `next`, `prev`, and `loop_` on each half-edge.
    pub fn add_loop(
        &mut self,
        half_edges: &[HalfEdgeHandle],
        face: FaceHandle,
    ) -> Result<LoopHandle, CoreError> {
        if half_edges.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "loop requires at least one half-edge".into(),
            ));
        }

        let lh = self.brep.loops.insert(Loop {
            half_edge: half_edges[0],
            face,
        });

        let n = half_edges.len();
        for i in 0..n {
            let curr = half_edges[i];
            let next = half_edges[(i + 1) % n];
            let prev = half_edges[(i + n - 1) % n];
            self.brep.half_edges[curr].next = next;
            self.brep.half_edges[curr].prev = prev;
            self.brep.half_edges[curr].loop_ = lh;
        }

        if self.brep.faces[face].outer_loop == LoopHandle::default() {
            self.brep.faces[face].outer_loop = lh;
        }

        Ok(lh)
    }

    /// Add an inner loop (hole) to the given face.
    ///
    /// Same wiring as `add_loop`, but appends to `inner_loops` instead
    /// of setting `outer_loop`.
    pub fn add_inner_loop(
        &mut self,
        half_edges: &[HalfEdgeHandle],
        face: FaceHandle,
    ) -> Result<LoopHandle, CoreError> {
        if half_edges.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "loop requires at least one half-edge".into(),
            ));
        }

        let lh = self.brep.loops.insert(Loop {
            half_edge: half_edges[0],
            face,
        });

        let n = half_edges.len();
        for i in 0..n {
            let curr = half_edges[i];
            let next = half_edges[(i + 1) % n];
            let prev = half_edges[(i + n - 1) % n];
            self.brep.half_edges[curr].next = next;
            self.brep.half_edges[curr].prev = prev;
            self.brep.half_edges[curr].loop_ = lh;
        }

        self.brep.faces[face].inner_loops.push(lh);
        Ok(lh)
    }
}

impl BRepBuilder {
    /// Consume the builder, validate the topology, and return the BRep.
    ///
    /// Runs the full validation suite. Returns `CoreError::Validation`
    /// if any check fails.
    pub fn build(self) -> Result<BRep, CoreError> {
        self.brep.validate().map_err(CoreError::from)?;
        Ok(self.brep)
    }

    /// Consume the builder without validation.
    ///
    /// Use when you want to run validation yourself (e.g., with `validate_all()`
    /// to collect all errors instead of failing on the first).
    pub fn build_unchecked(self) -> BRep {
        self.brep
    }

    /// Get a shared reference to the in-progress BRep.
    pub fn brep(&self) -> &BRep {
        &self.brep
    }

    /// Get a mutable reference to the in-progress BRep.
    pub fn brep_mut(&mut self) -> &mut BRep {
        &mut self.brep
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;

    #[test]
    fn builder_add_vertices() {
        let mut b = BRepBuilder::new(Tolerance::default());
        let v1 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
        let v2 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
        assert_ne!(v1, v2);
        assert_eq!(b.brep.vertices.len(), 2);
    }

    #[test]
    fn builder_add_solid_and_shell() {
        let mut b = BRepBuilder::new(Tolerance::default());
        let solid = b.add_solid();
        let shell = b.add_shell(solid);
        assert_eq!(b.brep.solids[solid].shell, shell);
        assert_eq!(b.brep.shells[shell].solid, solid);
    }

    #[test]
    fn builder_add_edge_creates_twin_pair() {
        let mut b = BRepBuilder::new(Tolerance::default());
        let v1 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
        let v2 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
        let (edge, he1, he2) = b.add_edge(v1, v2);
        assert_eq!(b.brep.half_edges[he1].twin, he2);
        assert_eq!(b.brep.half_edges[he2].twin, he1);
        assert_eq!(b.brep.half_edges[he1].origin, v1);
        assert_eq!(b.brep.half_edges[he2].origin, v2);
        assert_eq!(b.brep.half_edges[he1].edge, edge);
        assert_eq!(b.brep.half_edges[he2].edge, edge);
        assert_eq!(b.brep.vertices[v1].outgoing, he1);
    }

    #[test]
    fn builder_add_loop_wires_chain() {
        let mut b = BRepBuilder::new(Tolerance::default());
        let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
        let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
        let v2 = b.add_vertex(Point3::new(0.5, 1.0, 0.0));

        let (_, he01, _) = b.add_edge(v0, v1);
        let (_, he12, _) = b.add_edge(v1, v2);
        let (_, he20, _) = b.add_edge(v2, v0);

        let solid = b.add_solid();
        let shell = b.add_shell(solid);
        let face = b.add_face(shell);
        let loop_h = b.add_loop(&[he01, he12, he20], face).unwrap();

        assert_eq!(b.brep.half_edges[he01].next, he12);
        assert_eq!(b.brep.half_edges[he12].next, he20);
        assert_eq!(b.brep.half_edges[he20].next, he01);
        assert_eq!(b.brep.half_edges[he01].prev, he20);
        assert_eq!(b.brep.half_edges[he12].prev, he01);
        assert_eq!(b.brep.half_edges[he20].prev, he12);
        assert_eq!(b.brep.half_edges[he01].loop_, loop_h);
        assert_eq!(b.brep.loops[loop_h].half_edge, he01);
        assert_eq!(b.brep.loops[loop_h].face, face);
        assert_eq!(b.brep.faces[face].outer_loop, loop_h);
        assert_eq!(b.brep.faces[face].shell, shell);
    }

    #[test]
    fn builder_set_edge_curve() {
        use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
        let mut b = BRepBuilder::new(Tolerance::default());
        let v1 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
        let v2 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
        let (edge, _, _) = b.add_edge(v1, v2);
        b.set_edge_curve(
            edge,
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    direction: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
                },
                t_start: 0.0,
                t_end: 1.0,
            },
        );
        assert!(b.brep.edges[edge].curve.is_some());
    }

    #[test]
    fn builder_constructs_valid_cube() {
        let mut b = BRepBuilder::new(Tolerance::default());

        // 8 vertices
        let v = [
            b.add_vertex(Point3::new(0.0, 0.0, 0.0)),
            b.add_vertex(Point3::new(1.0, 0.0, 0.0)),
            b.add_vertex(Point3::new(1.0, 1.0, 0.0)),
            b.add_vertex(Point3::new(0.0, 1.0, 0.0)),
            b.add_vertex(Point3::new(0.0, 0.0, 1.0)),
            b.add_vertex(Point3::new(1.0, 0.0, 1.0)),
            b.add_vertex(Point3::new(1.0, 1.0, 1.0)),
            b.add_vertex(Point3::new(0.0, 1.0, 1.0)),
        ];

        // 12 edges
        let (_, he01, he10) = b.add_edge(v[0], v[1]);
        let (_, he12, he21) = b.add_edge(v[1], v[2]);
        let (_, he23, he32) = b.add_edge(v[2], v[3]);
        let (_, he30, he03) = b.add_edge(v[3], v[0]);
        let (_, he45, he54) = b.add_edge(v[4], v[5]);
        let (_, he56, he65) = b.add_edge(v[5], v[6]);
        let (_, he67, he76) = b.add_edge(v[6], v[7]);
        let (_, he74, he47) = b.add_edge(v[7], v[4]);
        let (_, he04, he40) = b.add_edge(v[0], v[4]);
        let (_, he15, he51) = b.add_edge(v[1], v[5]);
        let (_, he26, he62) = b.add_edge(v[2], v[6]);
        let (_, he37, he73) = b.add_edge(v[3], v[7]);

        let solid = b.add_solid();
        let shell = b.add_shell(solid);

        // 6 faces (CCW from outside for outward normals)
        // Front (z=0): v0->v3->v2->v1
        let f = b.add_face(shell);
        b.add_loop(&[he03, he32, he21, he10], f).unwrap();

        // Back (z=1): v4->v5->v6->v7
        let f = b.add_face(shell);
        b.add_loop(&[he45, he56, he67, he74], f).unwrap();

        // Bottom (y=0): v0->v1->v5->v4
        let f = b.add_face(shell);
        b.add_loop(&[he01, he15, he54, he40], f).unwrap();

        // Top (y=1): v3->v7->v6->v2
        let f = b.add_face(shell);
        b.add_loop(&[he37, he76, he65, he23], f).unwrap();

        // Left (x=0): v0->v4->v7->v3
        let f = b.add_face(shell);
        b.add_loop(&[he04, he47, he73, he30], f).unwrap();

        // Right (x=1): v1->v2->v6->v5
        let f = b.add_face(shell);
        b.add_loop(&[he12, he26, he62, he51], f).unwrap();

        // Attach geometry
        b.brep_mut().attach_geometry().unwrap();
        b.brep_mut().attach_planar_pcurves().unwrap();

        // Build with validation
        let brep = b.build().unwrap();
        assert_eq!(brep.vertices.len(), 8);
        assert_eq!(brep.edges.len(), 12);
        assert_eq!(brep.faces.len(), 6);
        assert_eq!(brep.solids.len(), 1);
    }
}
