//! Topological entities: Vertex, Edge, Loop, Face, Shell, Solid.

use crate::geometry::Point3;
use serde::{Deserialize, Serialize};

/// A topological vertex — a point in the BRep model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Vertex {
    /// Unique identifier.
    pub id: usize,
    /// The geometric point this vertex represents.
    pub point: Point3,
}

impl Vertex {
    /// Creates a new vertex with the given id and position.
    pub fn new(id: usize, point: Point3) -> Self {
        Self { id, point }
    }
}

/// A topological edge — connects two vertices.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge {
    /// Unique identifier.
    pub id: usize,
    /// ID of the start vertex.
    pub start_vertex: usize,
    /// ID of the end vertex.
    pub end_vertex: usize,
}

impl Edge {
    /// Creates a new edge between two vertices.
    pub fn new(id: usize, start_vertex: usize, end_vertex: usize) -> Self {
        Self {
            id,
            start_vertex,
            end_vertex,
        }
    }
}

/// A topological face — bounded by a loop of edges.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Face {
    /// Unique identifier.
    pub id: usize,
    /// IDs of the edges forming this face's boundary loop.
    pub edge_ids: Vec<usize>,
}

impl Face {
    /// Creates a new face with the given boundary edges.
    pub fn new(id: usize, edge_ids: Vec<usize>) -> Self {
        Self { id, edge_ids }
    }
}

/// A topological shell — a connected set of faces forming a closed surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Shell {
    /// Unique identifier.
    pub id: usize,
    /// IDs of the faces in this shell.
    pub face_ids: Vec<usize>,
}

impl Shell {
    /// Creates a new shell containing the given faces.
    pub fn new(id: usize, face_ids: Vec<usize>) -> Self {
        Self { id, face_ids }
    }
}

/// A topological solid — bounded by a shell.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Solid {
    /// Unique identifier.
    pub id: usize,
    /// ID of the bounding shell.
    pub shell_id: usize,
}

impl Solid {
    /// Creates a new solid bounded by the given shell.
    pub fn new(id: usize, shell_id: usize) -> Self {
        Self { id, shell_id }
    }
}

/// A complete BRep model containing all topological entities.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRepModel {
    /// All vertices in the model.
    pub vertices: Vec<Vertex>,
    /// All edges in the model.
    pub edges: Vec<Edge>,
    /// All faces in the model.
    pub faces: Vec<Face>,
    /// All shells in the model.
    pub shells: Vec<Shell>,
    /// All solids in the model.
    pub solids: Vec<Solid>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;

    #[test]
    fn vertex_creation() {
        let v = Vertex::new(0, Point3::new(0.0, 0.0, 0.0));
        assert_eq!(v.id, 0);
        assert_eq!(v.point.x, 0.0);
    }

    #[test]
    fn edge_references_two_vertices() {
        let e = Edge::new(0, 0, 1);
        assert_eq!(e.start_vertex, 0);
        assert_eq!(e.end_vertex, 1);
    }

    #[test]
    fn face_has_edge_loop() {
        let f = Face::new(0, vec![0, 1, 2, 3]);
        assert_eq!(f.edge_ids.len(), 4);
    }

    #[test]
    fn shell_contains_faces() {
        let s = Shell::new(0, vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(s.face_ids.len(), 6);
    }

    #[test]
    fn solid_contains_shell() {
        let s = Solid::new(0, 0);
        assert_eq!(s.shell_id, 0);
    }

    #[test]
    fn brep_model_counts() {
        let model = make_unit_cube();
        assert_eq!(model.vertices.len(), 8);
        assert_eq!(model.edges.len(), 12);
        assert_eq!(model.faces.len(), 6);
    }

    /// Helper: build a unit cube BRep model for testing.
    fn make_unit_cube() -> BRepModel {
        let vertices = vec![
            Vertex::new(0, Point3::new(0.0, 0.0, 0.0)),
            Vertex::new(1, Point3::new(1.0, 0.0, 0.0)),
            Vertex::new(2, Point3::new(1.0, 1.0, 0.0)),
            Vertex::new(3, Point3::new(0.0, 1.0, 0.0)),
            Vertex::new(4, Point3::new(0.0, 0.0, 1.0)),
            Vertex::new(5, Point3::new(1.0, 0.0, 1.0)),
            Vertex::new(6, Point3::new(1.0, 1.0, 1.0)),
            Vertex::new(7, Point3::new(0.0, 1.0, 1.0)),
        ];
        let edges = vec![
            Edge::new(0, 0, 1),
            Edge::new(1, 1, 2),
            Edge::new(2, 2, 3),
            Edge::new(3, 3, 0),
            Edge::new(4, 4, 5),
            Edge::new(5, 5, 6),
            Edge::new(6, 6, 7),
            Edge::new(7, 7, 4),
            Edge::new(8, 0, 4),
            Edge::new(9, 1, 5),
            Edge::new(10, 2, 6),
            Edge::new(11, 3, 7),
        ];
        let faces = vec![
            Face::new(0, vec![0, 1, 2, 3]),
            Face::new(1, vec![4, 5, 6, 7]),
            Face::new(2, vec![0, 9, 4, 8]),
            Face::new(3, vec![1, 10, 5, 9]),
            Face::new(4, vec![2, 11, 6, 10]),
            Face::new(5, vec![3, 8, 7, 11]),
        ];
        let shell = Shell::new(0, vec![0, 1, 2, 3, 4, 5]);
        let solid = Solid::new(0, 0);
        BRepModel {
            vertices,
            edges,
            faces,
            shells: vec![shell],
            solids: vec![solid],
        }
    }
}
