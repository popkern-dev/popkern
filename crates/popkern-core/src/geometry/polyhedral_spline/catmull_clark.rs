//! One Catmull–Clark refinement step on a [`ControlNet`].
//!
//! The polyhedral-net spline family covers specific local
//! configurations (all-quad vertex neighborhoods, the T0/T1/T2/Ngon
//! transition faces). A general cage is not guaranteed to decompose
//! into them — but one Catmull–Clark step normalizes any manifold cage
//! into an all-quad cage whose interior vertices all have all-quad
//! neighborhoods: original vertices keep their valence, edge points
//! have valence 4, and an n-gon's face point has valence n. Provided
//! the cage's valences and face arities are within the family's 3..=8
//! window, every quarter-cell of the refined cage is then served by a
//! `RegularQuad`/`ExtraordinaryStar` vertex spline. Valences or
//! arities above 8 are *preserved* by refinement, so repeating the
//! step never rescues them — callers should decline rather than loop.
//!
//! Semantically the refined net represents the same Catmull–Clark
//! limit surface the cage implies, so refinement changes the
//! realization of the design intent, not the intent itself. See
//! `docs/architecture/smooth-mesh-cage-normalization.md`.

use std::collections::HashMap;

use super::net::ControlNet;
use super::types::PolyhedralSplineError;
use crate::geometry::{Point3, Vec3};

/// An undirected cage edge with its incident faces.
struct CageEdge {
    a: usize,
    b: usize,
    faces: Vec<usize>,
}

impl ControlNet {
    /// One Catmull–Clark refinement step.
    ///
    /// Standard rules: face points at face centroids; interior edge
    /// points at `(a + b + F1 + F2) / 4` (boundary edges at the
    /// midpoint); interior vertex points at `(Q + 2R + (n-3)S) / n`
    /// with `Q` the mean incident face point, `R` the mean incident
    /// edge midpoint and `n` the valence (boundary vertices at
    /// `(prev + 6S + next) / 8` along the boundary polyline). Each
    /// original face contributes one quad per corner, preserving the
    /// face's winding.
    ///
    /// # Errors
    ///
    /// [`PolyhedralSplineError`] when the net is not a manifold cage
    /// (duplicate directed edge — also produced by inconsistent
    /// winding), or when a boundary vertex does not have exactly two
    /// boundary edges (a non-manifold boundary pinch, reported as
    /// [`PolyhedralSplineError::BoundaryPinch`]).
    pub fn catmull_clark_refined(&self) -> Result<ControlNet, PolyhedralSplineError> {
        // Validates manifoldness (duplicate directed edges) up front.
        let _ = self.halfedges()?;

        let face_points: Vec<Point3> = self
            .faces
            .iter()
            .map(|face| centroid(face.vertices.iter().map(|&v| self.vertices[v].point)))
            .collect();

        let (edges, edge_index) = undirected_edges(self);

        let edge_points: Vec<Point3> = edges
            .iter()
            .map(|edge| {
                let a = self.vertices[edge.a].point;
                let b = self.vertices[edge.b].point;
                if edge.faces.len() == 2 {
                    average(&[a, b, face_points[edge.faces[0]], face_points[edge.faces[1]]])
                } else {
                    average(&[a, b])
                }
            })
            .collect();

        let vertex_points = refined_vertex_points(self, &face_points, &edges)?;

        // New vertex layout: [vertex points][face points][edge points].
        let face_base = vertex_points.len();
        let edge_base = face_base + face_points.len();
        let mut vertices = vertex_points;
        vertices.extend(face_points.iter().copied());
        vertices.extend(edge_points.iter().copied());

        let mut faces = Vec::new();
        for (face_id, face) in self.faces.iter().enumerate() {
            let n = face.vertices.len();
            for corner in 0..n {
                let v = face.vertices[corner];
                let v_next = face.vertices[(corner + 1) % n];
                let v_prev = face.vertices[(corner + n - 1) % n];
                let edge_out = edge_base + edge_index[&undirected_key(v, v_next)];
                let edge_in = edge_base + edge_index[&undirected_key(v_prev, v)];
                faces.push(vec![v, edge_out, face_base + face_id, edge_in]);
            }
        }

        ControlNet::new(vertices, faces)
    }
}

fn undirected_key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

fn undirected_edges(net: &ControlNet) -> (Vec<CageEdge>, HashMap<(usize, usize), usize>) {
    let mut edges: Vec<CageEdge> = Vec::new();
    let mut index = HashMap::new();
    for (face_id, face) in net.faces.iter().enumerate() {
        let n = face.vertices.len();
        for corner in 0..n {
            let a = face.vertices[corner];
            let b = face.vertices[(corner + 1) % n];
            let key = undirected_key(a, b);
            let slot = *index.entry(key).or_insert_with(|| {
                edges.push(CageEdge {
                    a: key.0,
                    b: key.1,
                    faces: Vec::new(),
                });
                edges.len() - 1
            });
            edges[slot].faces.push(face_id);
        }
    }
    (edges, index)
}

fn refined_vertex_points(
    net: &ControlNet,
    face_points: &[Point3],
    edges: &[CageEdge],
) -> Result<Vec<Point3>, PolyhedralSplineError> {
    let mut incident_faces: Vec<Vec<usize>> = vec![Vec::new(); net.vertices.len()];
    for (face_id, face) in net.faces.iter().enumerate() {
        for &v in &face.vertices {
            incident_faces[v].push(face_id);
        }
    }
    let mut incident_edges: Vec<Vec<usize>> = vec![Vec::new(); net.vertices.len()];
    for (edge_id, edge) in edges.iter().enumerate() {
        incident_edges[edge.a].push(edge_id);
        incident_edges[edge.b].push(edge_id);
    }

    let mut out = Vec::with_capacity(net.vertices.len());
    for (v, vertex) in net.vertices.iter().enumerate() {
        let s = vertex.point;
        let boundary: Vec<usize> = incident_edges[v]
            .iter()
            .copied()
            .filter(|&e| edges[e].faces.len() < 2)
            .collect();
        if boundary.is_empty() {
            let n = incident_edges[v].len() as f64;
            let q = centroid(incident_faces[v].iter().map(|&f| face_points[f]));
            let r = centroid(incident_edges[v].iter().map(|&e| {
                average(&[
                    net.vertices[edges[e].a].point,
                    net.vertices[edges[e].b].point,
                ])
            }));
            let blended = (q - Point3::new(0.0, 0.0, 0.0))
                + (r - Point3::new(0.0, 0.0, 0.0)) * 2.0
                + (s - Point3::new(0.0, 0.0, 0.0)) * (n - 3.0);
            out.push(scale(blended, 1.0 / n));
        } else if boundary.len() == 2 {
            let other = |e: usize| {
                if edges[e].a == v {
                    edges[e].b
                } else {
                    edges[e].a
                }
            };
            let p0 = net.vertices[other(boundary[0])].point;
            let p1 = net.vertices[other(boundary[1])].point;
            let blended = (p0 - Point3::new(0.0, 0.0, 0.0))
                + (s - Point3::new(0.0, 0.0, 0.0)) * 6.0
                + (p1 - Point3::new(0.0, 0.0, 0.0));
            out.push(scale(blended, 1.0 / 8.0));
        } else {
            return Err(PolyhedralSplineError::BoundaryPinch {
                vertex: v,
                boundary_edges: boundary.len(),
            });
        }
    }
    Ok(out)
}

fn centroid(points: impl Iterator<Item = Point3>) -> Point3 {
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    let mut count = 0usize;
    for point in points {
        sum = sum + (point - Point3::new(0.0, 0.0, 0.0));
        count += 1;
    }
    scale(sum, 1.0 / count.max(1) as f64)
}

fn average(points: &[Point3]) -> Point3 {
    centroid(points.iter().copied())
}

fn scale(v: Vec3, factor: f64) -> Point3 {
    Point3::new(v.x * factor, v.y * factor, v.z * factor)
}
