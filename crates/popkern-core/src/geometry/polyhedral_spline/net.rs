use super::types::PolyhedralSplineError;
use crate::geometry::Point3;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ControlVertex {
    pub point: Point3,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ControlFace {
    pub vertices: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ControlNet {
    pub vertices: Vec<ControlVertex>,
    pub faces: Vec<ControlFace>,
}

#[derive(Debug, Clone)]
pub(super) struct HalfedgeNet<'a> {
    pub net: &'a ControlNet,
    halfedges: Vec<Halfedge>,
    outgoing: Vec<Vec<usize>>,
    face_halfedges: Vec<Vec<usize>>,
}

#[derive(Debug, Clone)]
struct Halfedge {
    from: usize,
    to: usize,
    face: usize,
    next: usize,
    prev: usize,
    opposite: Option<usize>,
}

impl ControlNet {
    pub fn new(
        vertices: Vec<Point3>,
        faces: Vec<Vec<usize>>,
    ) -> Result<Self, PolyhedralSplineError> {
        for (face, verts) in faces.iter().enumerate() {
            for &vertex in verts {
                if vertex >= vertices.len() {
                    return Err(PolyhedralSplineError::MissingVertex { face, vertex });
                }
            }
        }
        let vertices = vertices
            .into_iter()
            .map(|point| ControlVertex { point })
            .collect();
        let faces = faces
            .into_iter()
            .map(|vertices| ControlFace { vertices })
            .collect();
        Ok(Self { vertices, faces })
    }

    pub(super) fn halfedges(&self) -> Result<HalfedgeNet<'_>, PolyhedralSplineError> {
        HalfedgeNet::new(self)
    }

    pub fn from_obj_str(text: &str) -> Result<Self, PolyhedralSplineError> {
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        for (line_idx, line) in text.lines().enumerate() {
            let line_no = line_idx + 1;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let mut parts = trimmed.split_whitespace();
            match parts.next() {
                Some("v") => {
                    let coords: Vec<f64> = parts
                        .take(3)
                        .map(|part| {
                            part.parse::<f64>()
                                .map_err(|_| PolyhedralSplineError::ObjParse {
                                    line: line_no,
                                    reason: format!("invalid vertex coordinate '{part}'"),
                                })
                        })
                        .collect::<Result<_, _>>()?;
                    if coords.len() != 3 {
                        return Err(PolyhedralSplineError::ObjParse {
                            line: line_no,
                            reason: "vertex record needs three coordinates".to_string(),
                        });
                    }
                    vertices.push(Point3::new(coords[0], coords[1], coords[2]));
                }
                Some("f") => {
                    let face = parts
                        .map(|part| parse_obj_face_vertex(part, line_no))
                        .collect::<Result<Vec<_>, _>>()?;
                    if face.len() < 3 {
                        return Err(PolyhedralSplineError::ObjParse {
                            line: line_no,
                            reason: "face record needs at least three vertices".to_string(),
                        });
                    }
                    faces.push(face);
                }
                _ => {}
            }
        }
        Self::new(vertices, faces)
    }
}

fn parse_obj_face_vertex(part: &str, line: usize) -> Result<usize, PolyhedralSplineError> {
    let raw = part.split('/').next().unwrap_or_default();
    let index = raw
        .parse::<usize>()
        .map_err(|_| PolyhedralSplineError::ObjParse {
            line,
            reason: format!("invalid face vertex '{part}'"),
        })?;
    index
        .checked_sub(1)
        .ok_or_else(|| PolyhedralSplineError::ObjParse {
            line,
            reason: "OBJ face indices are 1-based".to_string(),
        })
}

impl<'a> HalfedgeNet<'a> {
    fn new(net: &'a ControlNet) -> Result<Self, PolyhedralSplineError> {
        let mut halfedges = Vec::new();
        let mut face_halfedges = Vec::with_capacity(net.faces.len());
        let mut directed = HashMap::<(usize, usize), usize>::new();

        for (face_id, face) in net.faces.iter().enumerate() {
            let mut ids = Vec::with_capacity(face.vertices.len());
            for index_in_face in 0..face.vertices.len() {
                let from = face.vertices[index_in_face];
                let to = face.vertices[(index_in_face + 1) % face.vertices.len()];
                if directed.contains_key(&(from, to)) {
                    return Err(PolyhedralSplineError::DuplicateDirectedEdge { from, to });
                }
                let id = halfedges.len();
                directed.insert((from, to), id);
                ids.push(id);
                halfedges.push(Halfedge {
                    from,
                    to,
                    face: face_id,
                    next: id,
                    prev: id,
                    opposite: None,
                });
            }
            face_halfedges.push(ids);
        }

        for ids in &face_halfedges {
            for (i, &h) in ids.iter().enumerate() {
                let n = ids.len();
                halfedges[h].next = ids[(i + 1) % n];
                halfedges[h].prev = ids[(i + n - 1) % n];
            }
        }

        for halfedge in &mut halfedges {
            let from = halfedge.from;
            let to = halfedge.to;
            halfedge.opposite = directed.get(&(to, from)).copied();
        }

        let mut outgoing = vec![Vec::new(); net.vertices.len()];
        for (h, halfedge) in halfedges.iter().enumerate() {
            outgoing[halfedge.from].push(h);
        }

        Ok(Self {
            net,
            halfedges,
            outgoing,
            face_halfedges,
        })
    }

    pub fn point(&self, vertex: usize) -> Point3 {
        self.net.vertices[vertex].point
    }

    pub fn face_valence(&self, face: usize) -> usize {
        self.net.faces[face].vertices.len()
    }

    pub fn vertex_valence(&self, vertex: usize) -> usize {
        self.outgoing[vertex].len()
    }

    pub fn is_boundary_vertex(&self, vertex: usize) -> bool {
        self.outgoing[vertex]
            .iter()
            .any(|&h| self.halfedges[h].opposite.is_none())
    }

    pub fn face_vertices(&self, face: usize) -> &[usize] {
        &self.net.faces[face].vertices
    }

    pub fn face_halfedges(&self, face: usize) -> &[usize] {
        &self.face_halfedges[face]
    }

    pub fn outgoing_halfedges(&self, vertex: usize) -> &[usize] {
        &self.outgoing[vertex]
    }

    pub fn halfedge_from_vertex(&self, halfedge: usize) -> usize {
        self.halfedges[halfedge].from
    }

    pub fn next_halfedge(&self, halfedge: usize) -> usize {
        self.halfedges[halfedge].next
    }

    pub fn ordered_outgoing(&self, vertex: usize) -> Vec<usize> {
        let outgoing = &self.outgoing[vertex];
        if outgoing.len() <= 1 {
            return outgoing.clone();
        }
        let mut ordered = Vec::with_capacity(outgoing.len());
        let mut seen = HashSet::new();
        let mut current = outgoing[0];
        for _ in 0..outgoing.len() {
            if !seen.insert(current) {
                break;
            }
            ordered.push(current);
            let prev = self.halfedges[current].prev;
            let Some(next) = self.halfedges[prev].opposite else {
                break;
            };
            if self.halfedges[next].from != vertex {
                break;
            }
            current = next;
        }
        for &h in outgoing {
            if seen.insert(h) {
                ordered.push(h);
            }
        }
        ordered
    }

    pub fn faces_around_vertex(&self, vertex: usize) -> Vec<usize> {
        self.ordered_outgoing(vertex)
            .into_iter()
            .map(|h| self.halfedges[h].face)
            .collect()
    }

    pub fn neighbor_faces_of_face(&self, face: usize) -> Vec<usize> {
        let mut faces = HashSet::new();
        for &vertex in self.face_vertices(face) {
            for h in &self.outgoing[vertex] {
                let f = self.halfedges[*h].face;
                if f != face {
                    faces.insert(f);
                }
            }
        }
        let mut out: Vec<usize> = faces.into_iter().collect();
        out.sort_unstable();
        out
    }

    pub fn all_faces_have_size(&self, faces: &[usize], size: usize) -> bool {
        faces.iter().all(|&face| self.face_valence(face) == size)
    }

    pub fn command_vertices(&self, start: usize, commands: &[u8]) -> Option<Vec<usize>> {
        self.command_vertices_with_end(start, commands)
            .map(|(vertices, _)| vertices)
    }

    pub fn command_vertices_with_end(
        &self,
        start: usize,
        commands: &[u8],
    ) -> Option<(Vec<usize>, usize)> {
        let mut h = start;
        let mut out = Vec::new();
        for command in commands {
            match command {
                1 => h = self.halfedges[h].next,
                2 => h = self.halfedges[h].prev,
                3 => h = self.halfedges[h].opposite?,
                4 => out.push(self.halfedges[h].from),
                _ => return None,
            }
        }
        Some((out, h))
    }
}
