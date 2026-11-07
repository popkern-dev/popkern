//! Tessellating a body of any kind (PK-1107): its faces' triangles
//! tagged by the regions each face lies between, its wire edges as
//! polylines, its acorns as points.
//!
//! A general body's faces do not all bound one solid: a lamina inside a
//! solid or a sheet's face has the void (or one region) on both sides
//! and is seen from both, a face between two solids is inside the body
//! and seen from neither, a wire edge has no face at all. The triangles
//! follow each face's folded normal (from its back to its front) as
//! [`crate::tessellate`] makes them; [`FaceRegions`] says which region
//! lies behind and in front, so a renderer draws a two-sided face from
//! both sides and an analysis reads one region's boundary oriented out
//! of it ([`GeneralMesh::region_triangles`]).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, RegionHandle, RegionKind, VertexHandle,
};

use crate::TessellationOptions;
use crate::error::MeshError;
use crate::mesh::TriMesh;

/// The regions a face lies between.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceRegions {
    /// The face.
    pub face: FaceHandle,
    /// The region behind it (its folded normal points away).
    pub back: RegionHandle,
    /// The region in front of it.
    pub front: RegionHandle,
    /// Seen from both sides: no solid behind or in front, or one region
    /// on both (a sheet's face, a lamina inside a solid).
    pub two_sided: bool,
    /// Between two solid regions: inside the body, seen from neither.
    pub internal: bool,
}

/// A wire edge sampled within the chord height.
#[derive(Debug, Clone, PartialEq)]
pub struct WirePolyline {
    /// The edge.
    pub edge: EdgeHandle,
    /// Points along it, from its first fin's origin to its end.
    pub points: Vec<Point3>,
}

/// A body of any kind, tessellated.
#[derive(Debug, Clone)]
pub struct GeneralMesh {
    /// The faces' triangles ([`crate::tessellate`]'s mesh).
    pub mesh: TriMesh,
    /// Per face of the mesh, the regions it lies between.
    pub faces: Vec<FaceRegions>,
    /// The wire edges.
    pub wires: Vec<WirePolyline>,
    /// The isolated vertices.
    pub acorns: Vec<(VertexHandle, Point3)>,
}

impl GeneralMesh {
    /// The regions triangle `i` lies between, `None` past the mesh.
    pub fn triangle_regions(&self, i: usize) -> Option<FaceRegions> {
        let face = *self.mesh.face_map.get(i)?;
        self.faces.iter().copied().find(|f| f.face == face)
    }

    /// The triangles bounding `region`, each with whether it must be
    /// reversed to face out of the region (it does as made when the
    /// region lies behind its face). A face with the region on both
    /// sides bounds nothing and is left out.
    pub fn region_triangles(&self, region: RegionHandle) -> Vec<(usize, bool)> {
        self.mesh
            .face_map
            .iter()
            .enumerate()
            .filter_map(|(i, &face)| {
                let f = self.faces.iter().find(|r| r.face == face)?;
                match (f.back == region, f.front == region) {
                    (true, false) => Some((i, false)),
                    (false, true) => Some((i, true)),
                    _ => None,
                }
            })
            .collect()
    }
}

/// Tessellate `brep`, whatever its kind: the faces as
/// [`crate::tessellate`] meshes them, each tagged with its regions,
/// the wire edges sampled within `options.chord_height`, the acorns.
///
/// # Errors
///
/// [`MeshError`] when a face does not tessellate.
pub fn tessellate_general(
    brep: &BRep,
    options: &TessellationOptions,
) -> Result<GeneralMesh, MeshError> {
    let mesh = if brep.faces.is_empty() {
        TriMesh::from_face_meshes(Vec::new())
    } else {
        crate::tessellate(brep, options)?
    };
    let solid = |r: RegionHandle| brep.region_kind(r) == Some(RegionKind::Solid);
    let faces = mesh
        .face_meshes
        .iter()
        .filter_map(|fm| {
            let sides = brep.face_sides(fm.face)?;
            let back = brep.shell_region(sides.back)?;
            let front = brep.shell_region(sides.front)?;
            Some(FaceRegions {
                face: fm.face,
                back,
                front,
                two_sided: back == front || (!solid(back) && !solid(front)),
                internal: back != front && solid(back) && solid(front),
            })
        })
        .collect();
    let wires = brep
        .wire_edges()
        .into_iter()
        .map(|edge| WirePolyline {
            edge,
            points: wire_points(brep, edge, options.chord_height),
        })
        .collect();
    let acorns = brep
        .acorns()
        .into_iter()
        .filter_map(|(_, v)| brep.vertices.get(v).map(|d| (v, d.position)))
        .collect();
    Ok(GeneralMesh {
        mesh,
        faces,
        wires,
        acorns,
    })
}

/// Deepest bisection of a wire's parameter range.
const MAX_WIRE_DEPTH: u32 = 12;

/// `edge` sampled by bisection until each chord's midpoint lies within
/// `chord` of the curve (a straight edge: its two ends).
fn wire_points(brep: &BRep, edge: EdgeHandle, chord: f64) -> Vec<Point3> {
    let Some((a, b)) = brep.edge_vertices(edge) else {
        return Vec::new();
    };
    let (pa, pb) = (brep.vertices[a].position, brep.vertices[b].position);
    let Some(c) = brep.edges.get(edge).and_then(|e| e.curve.as_ref()) else {
        return vec![pa, pb];
    };
    let at = |t: f64| c.geometry.point_at(t);
    let mut out = vec![at(c.t_start)];
    bisect(&at, c.t_start, c.t_end, chord.max(1e-12), 0, &mut out);
    out
}

fn bisect(
    at: &dyn Fn(f64) -> Point3,
    t0: f64,
    t1: f64,
    chord: f64,
    depth: u32,
    out: &mut Vec<Point3>,
) {
    let (p0, p1) = (at(t0), at(t1));
    let tm = 0.5 * (t0 + t1);
    let pm = at(tm);
    let mid = p0 + (p1 - p0) * 0.5;
    // A closed curve's two ends coincide: always split the first span.
    if depth < MAX_WIRE_DEPTH && (pm.distance_to(&mid) > chord || depth < 1) {
        bisect(at, t0, tm, chord, depth + 1, out);
        bisect(at, tm, t1, chord, depth + 1, out);
    } else {
        out.push(p1);
    }
}

#[cfg(test)]
#[path = "general_tests.rs"]
mod tests;
