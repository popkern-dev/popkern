use super::mask::{apply_inline_mask_to_patch, apply_mask_to_patches};
use super::net::HalfedgeNet;
use super::patch_subdivision::quarter_subdivided;
use super::tables;
use super::types::{
    PolyhedralPatch, PolyhedralSpline, PolyhedralSplineError, PolyhedralSplineFamily,
};
use crate::geometry::Point3;

pub fn construct_regular_quad(
    stencil: [Point3; 9],
) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    let family = PolyhedralSplineFamily::RegularQuad;
    let patch = apply_inline_mask_to_patch(&tables::REGULAR_MASK, family, &stencil)?;
    // The biquadratic mask spans the vertex's whole dual cell (its
    // corners are the four face centroids), while extraordinary-star
    // sectors tile per quarter-cell. Splitting at the parametric
    // midlines is exact (de Casteljau) and gives every component the
    // same quarter-cell boundary granularity, so mixed regular/star
    // complexes weld 1:1. See
    // docs/architecture/smooth-mesh-cage-normalization.md.
    spline(family, quarter_subdivided(&patch)?)
}

pub fn construct_two_triangles_two_quads(
    stencil: [Point3; 9],
) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    let family = PolyhedralSplineFamily::TwoTrianglesTwoQuads;
    let patch = apply_inline_mask_to_patch(&tables::REGULAR_MASK, family, &stencil)?;
    // Same dual-cell mask as RegularQuad, same quarter split (see
    // construct_regular_quad): the quarter sides then pair 1:1 with
    // regular neighbours AND with polar sector edges (whose endpoints
    // are spoke midpoints and triangle centroids — points the coarse
    // dual-cell side skipped over).
    spline(family, quarter_subdivided(&patch)?)
}

pub fn construct_extraordinary_star(
    valence: usize,
    stencil: &[Point3],
) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    let family = PolyhedralSplineFamily::ExtraordinaryStar { valence };
    let patches = apply_mask_to_patches(tables::extraordinary(valence)?, family, 3, 3, stencil)?;
    spline(family, patches)
}

pub fn construct_ngon(
    sides: usize,
    stencil: &[Point3],
) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    let family = PolyhedralSplineFamily::Ngon { sides };
    let patches = apply_mask_to_patches(tables::ngon(sides)?, family, 3, 3, stencil)?;
    spline(family, patches)
}

pub fn construct_polar(
    sectors: usize,
    stencil: &[Point3],
) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    let family = PolyhedralSplineFamily::Polar { sectors };
    let patches = apply_mask_to_patches(tables::polar(sectors)?, family, 3, 2, stencil)?;
    spline(family, patches)
}

pub fn construct_t0(stencil: [Point3; 14]) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    let family = PolyhedralSplineFamily::T0;
    let patches = apply_mask_to_patches(tables::t0()?, family, 3, 3, &stencil)?;
    spline(family, patches)
}

pub fn construct_t1(stencil: [Point3; 18]) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    let family = PolyhedralSplineFamily::T1;
    let patches = apply_mask_to_patches(tables::t1()?, family, 3, 3, &stencil)?;
    spline(family, patches)
}

pub fn construct_t2(stencil: [Point3; 20]) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    let family = PolyhedralSplineFamily::T2;
    let patches = apply_mask_to_patches(tables::t2()?, family, 3, 3, &stencil)?;
    spline(family, patches)
}

fn spline(
    family: PolyhedralSplineFamily,
    patches: Vec<PolyhedralPatch>,
) -> Result<PolyhedralSpline, PolyhedralSplineError> {
    for patch in &patches {
        for (index, point) in patch.control_points.iter().enumerate() {
            if !(point.x.is_finite() && point.y.is_finite() && point.z.is_finite()) {
                return Err(PolyhedralSplineError::NonFinitePoint { family, index });
            }
        }
    }
    Ok(PolyhedralSpline { family, patches })
}

pub(super) fn try_ngon_face(
    net: &HalfedgeNet<'_>,
    face: usize,
) -> Result<Option<PolyhedralSpline>, PolyhedralSplineError> {
    let sides = net.face_valence(face);
    if !(3..=8).contains(&sides) || sides == 4 {
        return Ok(None);
    }
    if !net
        .face_vertices(face)
        .iter()
        .all(|&vertex| net.vertex_valence(vertex) == 4)
    {
        return Ok(None);
    }
    let neighbors = net.neighbor_faces_of_face(face);
    if !net.all_faces_have_size(&neighbors, 4) {
        return Ok(None);
    }

    let mut stencil = vec![Point3::new(0.0, 0.0, 0.0); 4 * sides];
    for (wing, &h) in net.face_halfedges(face).iter().enumerate() {
        let ids = [3 + 4 * wing, 1 + 4 * wing, 4 * wing, 2 + 4 * wing];
        let vertices = net
            .command_vertices(h, &[3, 1, 4, 3, 4, 2, 4, 2, 4])
            .ok_or(PolyhedralSplineError::StencilTopology {
                family: PolyhedralSplineFamily::Ngon { sides },
                reason: "n-gon stencil crosses a boundary or non-manifold edge",
            })?;
        for (slot, vertex) in ids.into_iter().zip(vertices) {
            stencil[slot] = net.point(vertex);
        }
    }
    construct_ngon(sides, &stencil).map(Some)
}

pub(super) fn try_t0_face(
    net: &HalfedgeNet<'_>,
    face: usize,
) -> Result<Option<PolyhedralSpline>, PolyhedralSplineError> {
    if net.face_valence(face) != 3 {
        return Ok(None);
    }
    let verts = net.face_vertices(face);
    let n5 = verts
        .iter()
        .filter(|&&v| net.vertex_valence(v) == 5)
        .count();
    let n4 = verts
        .iter()
        .filter(|&&v| net.vertex_valence(v) == 4)
        .count();
    let neighbors = net.neighbor_faces_of_face(face);
    if n5 != 1 || n4 != 2 || neighbors.len() != 7 || !net.all_faces_have_size(&neighbors, 4) {
        return Ok(None);
    }

    let mut tri = verts.to_vec();
    while net.vertex_valence(tri[0]) != 5 {
        tri.rotate_left(1);
    }
    let mut stencil = [Point3::new(0.0, 0.0, 0.0); 14];
    let t8 = tri[2];
    for (slot, vertex) in [4, 7, 8].into_iter().zip(tri) {
        stencil[slot] = net.point(vertex);
    }

    let start = net
        .face_halfedges(face)
        .iter()
        .copied()
        .find(|&h| net.halfedge_from_vertex(h) == t8)
        .ok_or(PolyhedralSplineError::StencilTopology {
            family: PolyhedralSplineFamily::T0,
            reason: "T0 central triangle has no expected 8->4 halfedge",
        })?;

    let mut h = start;
    let (verts, end) = net.command_vertices_with_end(h, &[3, 2, 4]).ok_or(
        PolyhedralSplineError::StencilTopology {
            family: PolyhedralSplineFamily::T0,
            reason: "T0 side stencil crosses a boundary edge",
        },
    )?;
    stencil[5] = net.point(verts[0]);
    h = end;

    let top = repeated(&[3, 1, 1, 4, 1, 4], 2);
    let (verts, end) =
        net.command_vertices_with_end(h, &top)
            .ok_or(PolyhedralSplineError::StencilTopology {
                family: PolyhedralSplineFamily::T0,
                reason: "T0 top stencil crosses a boundary edge",
            })?;
    for (slot, vertex) in [2, 1, 0, 3].into_iter().zip(verts) {
        stencil[slot] = net.point(vertex);
    }
    h = end;

    let bottom = repeated(&[3, 1, 1, 4, 3, 1, 1, 4, 1, 4], 2);
    let (verts, _) = net.command_vertices_with_end(h, &bottom).ok_or(
        PolyhedralSplineError::StencilTopology {
            family: PolyhedralSplineFamily::T0,
            reason: "T0 bottom stencil crosses a boundary edge",
        },
    )?;
    for (slot, vertex) in [6, 10, 11, 12, 13, 9].into_iter().zip(verts) {
        stencil[slot] = net.point(vertex);
    }
    construct_t0(stencil).map(Some)
}

pub(super) fn try_t1_face(
    net: &HalfedgeNet<'_>,
    face: usize,
) -> Result<Option<PolyhedralSpline>, PolyhedralSplineError> {
    if net.face_valence(face) != 5 {
        return Ok(None);
    }
    let verts = net.face_vertices(face);
    let n3 = verts
        .iter()
        .filter(|&&v| net.vertex_valence(v) == 3)
        .count();
    let n4 = verts
        .iter()
        .filter(|&&v| net.vertex_valence(v) == 4)
        .count();
    let neighbors = net.neighbor_faces_of_face(face);
    if n3 != 1 || n4 != 4 || neighbors.len() != 9 || !net.all_faces_have_size(&neighbors, 4) {
        return Ok(None);
    }
    let mut pen = verts.to_vec();
    while net.vertex_valence(pen[0]) != 3 {
        pen.rotate_left(1);
    }
    let mut stencil = [Point3::new(0.0, 0.0, 0.0); 18];
    for (slot, vertex) in [10, 9, 5, 6, 11].into_iter().zip(pen.iter().copied()) {
        stencil[slot] = net.point(vertex);
    }
    let pen_set = pen;
    let mut current = None;
    for &h in net.outgoing_halfedges(pen_set[0]) {
        let candidate = net.halfedge_from_vertex(net.next_halfedge(h));
        if !pen_set.contains(&candidate) {
            stencil[15] = net.point(candidate);
            current = Some(h);
            break;
        }
    }
    let mut h = current.ok_or(PolyhedralSplineError::StencilTopology {
        family: PolyhedralSplineFamily::T1,
        reason: "T1 could not find outgoing 10->15 halfedge",
    })?;
    for ids in [[14, 13, 8], [4, 0, 1], [2, 3, 7], [12, 17, 16]] {
        let (verts, end) = net
            .command_vertices_with_end(h, &[1, 1, 4, 3, 1, 1, 4, 1, 4, 3])
            .ok_or(PolyhedralSplineError::StencilTopology {
                family: PolyhedralSplineFamily::T1,
                reason: "T1 corner stencil crosses a boundary edge",
            })?;
        for (slot, vertex) in ids.into_iter().zip(verts) {
            stencil[slot] = net.point(vertex);
        }
        h = end;
    }
    construct_t1(stencil).map(Some)
}

pub(super) fn try_t2_face(
    net: &HalfedgeNet<'_>,
    face: usize,
) -> Result<Option<PolyhedralSpline>, PolyhedralSplineError> {
    if net.face_valence(face) != 6 {
        return Ok(None);
    }
    let mut vals: Vec<usize> = net
        .face_vertices(face)
        .iter()
        .map(|&v| net.vertex_valence(v))
        .collect();
    if vals.iter().all(|&v| v != 3) {
        return Ok(None);
    }
    while vals.first() != Some(&3) {
        vals.rotate_left(1);
    }
    if vals != [3, 4, 4, 4, 3, 4] && vals != [3, 4, 3, 4, 4, 4] {
        return Ok(None);
    }
    let neighbors = net.neighbor_faces_of_face(face);
    if !net.all_faces_have_size(&neighbors, 4) {
        return Ok(None);
    }
    let mut hex = net.face_vertices(face).to_vec();
    while !(net.vertex_valence(hex[0]) == 3 && net.vertex_valence(hex[2]) == 3) {
        hex.rotate_left(1);
    }
    let mut stencil = [Point3::new(0.0, 0.0, 0.0); 20];
    let v6 = hex[4];
    for (slot, vertex) in [9, 11, 12, 13, 6, 5].into_iter().zip(hex) {
        stencil[slot] = net.point(vertex);
    }
    let start = net
        .face_halfedges(face)
        .iter()
        .copied()
        .find(|&h| net.halfedge_from_vertex(h) == v6)
        .ok_or(PolyhedralSplineError::StencilTopology {
            family: PolyhedralSplineFamily::T2,
            reason: "T2 central hexagon has no expected 6->5 halfedge",
        })?;
    let (_, mut h) = net.command_vertices_with_end(start, &[3, 2, 3]).ok_or(
        PolyhedralSplineError::StencilTopology {
            family: PolyhedralSplineFamily::T2,
            reason: "T2 initial halfedge move crosses a boundary edge",
        },
    )?;
    for ids in [
        vec![1, 0],
        vec![4],
        vec![8],
        vec![10, 15],
        vec![16],
        vec![17],
        vec![18, 19],
        vec![14],
        vec![7, 3],
        vec![2],
    ] {
        let commands = if ids.len() == 2 {
            &[1, 4, 1, 4, 1, 3][..]
        } else {
            &[1, 4, 1, 3][..]
        };
        let (verts, end) = net.command_vertices_with_end(h, commands).ok_or(
            PolyhedralSplineError::StencilTopology {
                family: PolyhedralSplineFamily::T2,
                reason: "T2 ring stencil crosses a boundary edge",
            },
        )?;
        for (slot, vertex) in ids.into_iter().zip(verts) {
            stencil[slot] = net.point(vertex);
        }
        h = end;
    }
    construct_t2(stencil).map(Some)
}

pub(super) fn try_regular_vertex(
    net: &HalfedgeNet<'_>,
    vertex: usize,
) -> Result<Option<PolyhedralSpline>, PolyhedralSplineError> {
    if net.vertex_valence(vertex) != 4 || net.is_boundary_vertex(vertex) {
        return Ok(None);
    }
    let faces = net.faces_around_vertex(vertex);
    if faces.len() != 4 || !net.all_faces_have_size(&faces, 4) {
        return Ok(None);
    }
    let stencil = regular_like_stencil(net, vertex, PolyhedralSplineFamily::RegularQuad)?;
    construct_regular_quad(stencil).map(Some)
}

pub(super) fn try_two_triangles_two_quads_vertex(
    net: &HalfedgeNet<'_>,
    vertex: usize,
) -> Result<Option<PolyhedralSpline>, PolyhedralSplineError> {
    if net.vertex_valence(vertex) != 4 || net.is_boundary_vertex(vertex) {
        return Ok(None);
    }
    let faces = net.faces_around_vertex(vertex);
    if faces.len() != 4 {
        return Ok(None);
    }
    let triangles = faces.iter().filter(|&&f| net.face_valence(f) == 3).count();
    let quads = faces.iter().filter(|&&f| net.face_valence(f) == 4).count();
    if triangles != 2 || quads != 2 {
        return Ok(None);
    }
    let stencil = regular_like_stencil(net, vertex, PolyhedralSplineFamily::TwoTrianglesTwoQuads)?;
    construct_two_triangles_two_quads(stencil).map(Some)
}

pub(super) fn try_extraordinary_vertex(
    net: &HalfedgeNet<'_>,
    vertex: usize,
) -> Result<Option<PolyhedralSpline>, PolyhedralSplineError> {
    let valence = net.vertex_valence(vertex);
    if !(3..=8).contains(&valence) || valence == 4 || net.is_boundary_vertex(vertex) {
        return Ok(None);
    }
    let faces = net.faces_around_vertex(vertex);
    if !net.all_faces_have_size(&faces, 4) {
        return Ok(None);
    }
    let mut stencil = Vec::with_capacity(2 * valence + 1);
    for h in net.ordered_outgoing(vertex) {
        let vertices = net.command_vertices(h, &[1, 4, 1, 4]).ok_or(
            PolyhedralSplineError::StencilTopology {
                family: PolyhedralSplineFamily::ExtraordinaryStar { valence },
                reason: "extraordinary star stencil crosses a boundary edge",
            },
        )?;
        stencil.extend(vertices.into_iter().map(|v| net.point(v)));
    }
    stencil.push(net.point(vertex));
    construct_extraordinary_star(valence, &stencil).map(Some)
}

pub(super) fn try_polar_vertex(
    net: &HalfedgeNet<'_>,
    vertex: usize,
) -> Result<Option<PolyhedralSpline>, PolyhedralSplineError> {
    let sectors = net.vertex_valence(vertex);
    if !(3..=8).contains(&sectors) || net.is_boundary_vertex(vertex) {
        return Ok(None);
    }
    let faces = net.faces_around_vertex(vertex);
    if !net.all_faces_have_size(&faces, 3) {
        return Ok(None);
    }
    let mut stencil = Vec::with_capacity(sectors + 1);
    stencil.push(net.point(vertex));
    for h in net.ordered_outgoing(vertex) {
        let vertices =
            net.command_vertices(h, &[1, 4])
                .ok_or(PolyhedralSplineError::StencilTopology {
                    family: PolyhedralSplineFamily::Polar { sectors },
                    reason: "polar stencil crosses a boundary edge",
                })?;
        stencil.push(net.point(vertices[0]));
    }
    construct_polar(sectors, &stencil).map(Some)
}

fn regular_like_stencil(
    net: &HalfedgeNet<'_>,
    vertex: usize,
    family: PolyhedralSplineFamily,
) -> Result<[Point3; 9], PolyhedralSplineError> {
    let mut stencil = [Point3::new(0.0, 0.0, 0.0); 9];
    stencil[4] = net.point(vertex);
    let mut h = net.ordered_outgoing(vertex)[0];
    for ids in [[1, 0], [3, 6], [7, 8], [5, 2]] {
        let (vertices, end) = net
            .command_vertices_with_end(h, &[1, 4, 1, 4, 1, 3])
            .ok_or(PolyhedralSplineError::StencilTopology {
                family,
                reason: "regular-like stencil crosses a boundary edge",
            })?;
        for (slot, vertex) in ids.into_iter().zip(vertices) {
            stencil[slot] = net.point(vertex);
        }
        h = end;
    }
    Ok(stencil)
}

fn repeated(commands: &[u8], count: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(commands.len() * count);
    for _ in 0..count {
        out.extend_from_slice(commands);
    }
    out
}
