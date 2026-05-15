use super::construct::{
    try_extraordinary_vertex, try_ngon_face, try_polar_vertex, try_regular_vertex, try_t0_face,
    try_t1_face, try_t2_face, try_two_triangles_two_quads_vertex,
};
use super::net::ControlNet;
use super::types::{PolyhedralSpline, PolyhedralSplineError};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ControlNetSplineOrigin {
    Face { face: usize },
    Vertex { vertex: usize },
}

#[derive(Debug, Clone, PartialEq)]
pub struct ControlNetSpline {
    pub origin: ControlNetSplineOrigin,
    pub spline: PolyhedralSpline,
}

pub fn construct_from_control_net(
    net: &ControlNet,
) -> Result<Vec<PolyhedralSpline>, PolyhedralSplineError> {
    Ok(construct_from_control_net_with_origins(net)?
        .into_iter()
        .map(|item| item.spline)
        .collect())
}

pub fn construct_from_control_net_with_origins(
    net: &ControlNet,
) -> Result<Vec<ControlNetSpline>, PolyhedralSplineError> {
    let halfedges = net.halfedges()?;
    let mut splines = Vec::new();

    for face in 0..net.faces.len() {
        if let Some(spline) = try_t0_face(&halfedges, face)? {
            splines.push(face_spline(face, spline));
            continue;
        }
        if let Some(spline) = try_t1_face(&halfedges, face)? {
            splines.push(face_spline(face, spline));
            continue;
        }
        if let Some(spline) = try_t2_face(&halfedges, face)? {
            splines.push(face_spline(face, spline));
            continue;
        }
        if let Some(spline) = try_ngon_face(&halfedges, face)? {
            splines.push(face_spline(face, spline));
        }
    }

    for vertex in 0..net.vertices.len() {
        if let Some(spline) = try_regular_vertex(&halfedges, vertex)? {
            splines.push(vertex_spline(vertex, spline));
            continue;
        }
        if let Some(spline) = try_extraordinary_vertex(&halfedges, vertex)? {
            splines.push(vertex_spline(vertex, spline));
            continue;
        }
        if let Some(spline) = try_two_triangles_two_quads_vertex(&halfedges, vertex)? {
            splines.push(vertex_spline(vertex, spline));
            continue;
        }
        if let Some(spline) = try_polar_vertex(&halfedges, vertex)? {
            splines.push(vertex_spline(vertex, spline));
        }
    }

    Ok(splines)
}

fn face_spline(face: usize, spline: PolyhedralSpline) -> ControlNetSpline {
    ControlNetSpline {
        origin: ControlNetSplineOrigin::Face { face },
        spline,
    }
}

fn vertex_spline(vertex: usize, spline: PolyhedralSpline) -> ControlNetSpline {
    ControlNetSpline {
        origin: ControlNetSplineOrigin::Vertex { vertex },
        spline,
    }
}
