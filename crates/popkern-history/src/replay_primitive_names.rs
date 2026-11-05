//! Role names for a primitive built in its local frame (PK-1083; see
//! `replay_primitive`). Faces first, so edges and vertices can name the
//! faces they bound as ancestors.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, NameId, RoleCategory};
use popkern_core::naming_criteria::{
    GeometricCriteria, capture_edge, capture_face, capture_vertex,
};
use popkern_core::topology::{BRep, FaceHandle};
use popkern_math::PortableMath;

/// Angle about +z in `[0, 2π)`, and height: the order round the axis,
/// then along it.
fn station(direction: Vec3, z: f64) -> (f64, f64) {
    let a = direction
        .y
        .patan2(direction.x)
        .rem_euclid(std::f64::consts::TAU);
    // Quantized so round-off never reorders entities at one station.
    ((a * 1e9).round() / 1e9, (z * 1e9).round() / 1e9)
}

fn mean(points: &[Point3]) -> Point3 {
    let n = points.len().max(1) as f64;
    let s = points.iter().fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| {
        acc + (*p - Point3::new(0.0, 0.0, 0.0))
    });
    Point3::new(0.0, 0.0, 0.0) + s * (1.0 / n)
}

/// A planar face's outward normal.
fn plane_normal(brep: &BRep, face: FaceHandle) -> Option<Vec3> {
    let sof = brep.faces[face].surface.as_ref()?;
    match &sof.geometry {
        SurfaceGeometry::Plane { normal, .. } => {
            Some(if sof.orientation { *normal } else { -*normal })
        }
        _ => None,
    }
}

fn mint(
    ctx: &mut TrackingContext,
    category: RoleCategory,
    index: usize,
    ancestors: Vec<NameId>,
) -> EntityName {
    let name = EntityName {
        operation: ctx.operation,
        role: EntityRole {
            category,
            index: index as u32,
        },
        ancestors,
    };
    ctx.record_generated(name.clone());
    name
}

/// Name every face, edge and vertex of the local primitive. `sector`:
/// the planar faces through the axis are a sector's start (`Cap 2`, at
/// angle 0) and end (`Cap 3`).
pub(super) fn name_by_role(brep: &mut BRep, ctx: &mut TrackingContext, sector: bool) {
    let z = Vec3::new(0.0, 0.0, 1.0);
    let mut caps: Vec<(u32, FaceHandle)> = Vec::new();
    let mut laterals: Vec<((f64, f64), FaceHandle)> = Vec::new();
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    for fh in faces {
        let c = mean(&brep.face_boundary_sample_points(fh, 8));
        let n = plane_normal(brep, fh);
        match n {
            Some(n) if n.dot(&z) < -0.999_999 => caps.push((0, fh)),
            Some(n) if n.dot(&z) > 0.999_999 => caps.push((1, fh)),
            Some(_) if sector => {
                let (a, _) = station(c - Point3::new(0.0, 0.0, 0.0), c.z);
                caps.push((
                    if a < 1e-6 || a > std::f64::consts::TAU - 1e-6 {
                        2
                    } else {
                        3
                    },
                    fh,
                ));
            }
            Some(n) => laterals.push((station(n, c.z), fh)),
            None => laterals.push((station(c - Point3::new(0.0, 0.0, 0.0), c.z), fh)),
        }
    }
    laterals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut face_names = std::collections::HashMap::new();
    let roles = caps
        .iter()
        .map(|(i, fh)| (RoleCategory::Cap, *i as usize, *fh))
        .chain(
            laterals
                .iter()
                .enumerate()
                .map(|(i, (_, fh))| (RoleCategory::Lateral, i, *fh)),
        )
        .collect::<Vec<_>>();
    for (category, index, fh) in roles {
        let name = mint(ctx, category, index, Vec::new());
        let criteria = capture_face(brep, fh).map(GeometricCriteria::Face);
        let id = brep.name_pool.intern(name, criteria);
        brep.faces[fh].name = Some(id);
        face_names.insert(fh, id);
    }
    let face_of = |brep: &BRep, he| {
        brep.half_edges
            .get(he)
            .and_then(|h| brep.loops.get(h.loop_))
            .map(|l| l.face)
    };
    // Edges: by the names of the faces they bound, then round and
    // along the axis.
    let mut edges: Vec<(Vec<NameId>, (f64, f64), _)> = brep
        .edges
        .iter()
        .map(|(eh, e)| {
            let he = e.half_edge;
            let twin = brep.half_edges.get(he).map(|h| h.twin);
            let mut names: Vec<NameId> = [Some(he), twin]
                .into_iter()
                .flatten()
                .filter_map(|h| face_of(brep, h))
                .filter_map(|f| face_names.get(&f).copied())
                .collect();
            names.sort_by_key(|n| n.0);
            names.dedup();
            let mid = e
                .curve
                .as_ref()
                .map(|c| c.geometry.point_at((c.t_start + c.t_end) / 2.0))
                .unwrap_or(Point3::new(0.0, 0.0, 0.0));
            (names, station(mid - Point3::new(0.0, 0.0, 0.0), mid.z), eh)
        })
        .collect();
    edges.sort_by(|a, b| {
        a.0.iter()
            .map(|n| n.0)
            .cmp(b.0.iter().map(|n| n.0))
            .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });
    for (i, (ancestors, _, eh)) in edges.into_iter().enumerate() {
        let name = mint(ctx, RoleCategory::Intersection, i, ancestors);
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let id = brep.name_pool.intern(name, criteria);
        brep.edges[eh].name = Some(id);
    }
    // Vertices: by the faces round them, then by position.
    let mut around: std::collections::HashMap<_, Vec<NameId>> = std::collections::HashMap::new();
    for (he, h) in brep.half_edges.iter() {
        if let Some(id) = face_of(brep, he).and_then(|f| face_names.get(&f).copied()) {
            around.entry(h.origin).or_default().push(id);
        }
    }
    let mut vertices: Vec<(Vec<NameId>, (f64, f64), _)> = brep
        .vertices
        .iter()
        .map(|(vh, v)| {
            let mut names = around.remove(&vh).unwrap_or_default();
            names.sort_by_key(|n| n.0);
            names.dedup();
            let p = v.position;
            (names, station(p - Point3::new(0.0, 0.0, 0.0), p.z), vh)
        })
        .collect();
    vertices.sort_by(|a, b| {
        a.0.iter()
            .map(|n| n.0)
            .cmp(b.0.iter().map(|n| n.0))
            .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });
    for (i, (ancestors, _, vh)) in vertices.into_iter().enumerate() {
        let name = mint(ctx, RoleCategory::Junction, i, ancestors);
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let id = brep.name_pool.intern(name, criteria);
        brep.vertices[vh].name = Some(id);
    }
}
