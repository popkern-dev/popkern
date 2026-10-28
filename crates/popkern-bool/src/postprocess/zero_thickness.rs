//! Zero-thickness self-contact in a served boolean result (PK-964).
//!
//! A manifold B-rep never puts two of its vertices, or two of its
//! edges, in the same place. When a boolean's operands touch along an
//! edge or at a point, every route that can serve the result has had to
//! split the contact into coincident copies to keep the half-edge
//! structure two-manifold: two bodies sharing one point position, two
//! edges sharing one curve. That body is non-manifold in fact — it
//! fails `validate_all`, a slicer opens it — and under decision D1 the
//! boolean declines with the contacts instead of serving it.
//!
//! Detection is direct and route-independent: distinct vertices within
//! the linear tolerance of each other that no edge joins, and distinct
//! edges whose end points and midpoints coincide.
//!
//! Contacts that share no topology at all leave the result as separate
//! shells that merely touch, which neither test sees: those are
//! `shell_contact`'s.

use std::collections::{BTreeSet, HashMap};

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use popkern_core::topology::ContactEntity;

use crate::error::{ContactKind, NonManifoldContact};

fn cell(p: Point3, size: f64) -> (i64, i64, i64) {
    let q = |x: f64| (x / size).floor() as i64;
    (q(p.x), q(p.y), q(p.z))
}

fn edge_ends(brep: &BRep, e: EdgeHandle) -> Option<(VertexHandle, VertexHandle)> {
    let edge = brep.edges.get(e)?;
    let he = brep.half_edges.get(edge.half_edge)?;
    let tw = brep.half_edges.get(he.twin)?;
    Some((he.origin, tw.origin))
}

fn edge_mid(brep: &BRep, e: EdgeHandle) -> Option<Point3> {
    let edge = brep.edges.get(e)?;
    if let Some(c) = &edge.curve {
        return Some(c.geometry.point_at(0.5 * (c.t_start + c.t_end)));
    }
    let (a, b) = edge_ends(brep, e)?;
    let (pa, pb) = (
        brep.vertices.get(a)?.position,
        brep.vertices.get(b)?.position,
    );
    Some(pa + (pb - pa) * 0.5)
}

/// Every zero-thickness self-contact of `brep`, in handle order.
pub(crate) fn zero_thickness_contacts(
    brep: &BRep,
    tolerance: &Tolerance,
) -> Vec<NonManifoldContact> {
    zero_thickness_findings(brep, tolerance)
        .into_iter()
        .map(|(c, _)| c)
        .collect()
}

/// [`zero_thickness_contacts`] with the entities that meet at each.
pub(crate) fn zero_thickness_findings(
    brep: &BRep,
    tolerance: &Tolerance,
) -> Vec<(NonManifoldContact, Vec<ContactEntity>)> {
    let tol = tolerance.linear();
    let size = (tol * 4.0).max(1e-300);
    let mut out = Vec::new();

    // Vertices joined by an edge are allowed to be close (a short edge
    // is a feature, not a contact).
    let mut joined: BTreeSet<(VertexHandle, VertexHandle)> = BTreeSet::new();
    for e in brep.edges.keys() {
        if let Some((a, b)) = edge_ends(brep, e) {
            joined.insert((a.min(b), a.max(b)));
        }
    }
    let mut grid: HashMap<(i64, i64, i64), Vec<VertexHandle>> = HashMap::new();
    let mut vertices: Vec<VertexHandle> = brep.vertices.keys().collect();
    vertices.sort();
    let mut touched = BTreeSet::new();
    for &v in &vertices {
        let p = brep.vertices[v].position;
        let (cx, cy, cz) = cell(p, size);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(list) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &w in list {
                        if brep.vertices[w].position.distance_to(&p) <= tol
                            && !joined.contains(&(v.min(w), v.max(w)))
                            && touched.insert(w)
                        {
                            out.push((
                                NonManifoldContact {
                                    kind: ContactKind::Vertex,
                                    location: brep.vertices[w].position,
                                },
                                vec![ContactEntity::Vertex(w), ContactEntity::Vertex(v)],
                            ));
                        }
                    }
                }
            }
        }
        grid.entry((cx, cy, cz)).or_default().push(v);
    }

    // Coincident edges: same end vertices (or coincident ones), same
    // midpoint.
    let mut by_ends: HashMap<(VertexHandle, VertexHandle), Vec<EdgeHandle>> = HashMap::new();
    let mut edges: Vec<EdgeHandle> = brep.edges.keys().collect();
    edges.sort();
    for &e in &edges {
        if let Some((a, b)) = edge_ends(brep, e) {
            by_ends.entry((a.min(b), a.max(b))).or_default().push(e);
        }
    }
    let mut groups: Vec<&Vec<EdgeHandle>> = by_ends.values().filter(|g| g.len() > 1).collect();
    groups.sort();
    for group in groups {
        let mids: Vec<Option<Point3>> = group.iter().map(|&e| edge_mid(brep, e)).collect();
        for i in 0..group.len() {
            for j in (i + 1)..group.len() {
                if let (Some(a), Some(b)) = (mids[i], mids[j])
                    && a.distance_to(&b) <= tol
                {
                    out.push((
                        NonManifoldContact {
                            kind: ContactKind::Edge,
                            location: a,
                        },
                        vec![ContactEntity::Edge(group[i]), ContactEntity::Edge(group[j])],
                    ));
                }
            }
        }
    }
    out
}

/// Whether `point` lies on a face of `body` among `aabbs`: within the
/// linear tolerance of the face's surface, inside or on its trim. The
/// ray vote of point classification has no on-boundary verdict of its
/// own (a point on a contact line votes either way), so contact points
/// are recognised geometrically.
pub(crate) fn on_body_boundary(
    point: &Point3,
    body: &BRep,
    aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> bool {
    aabbs
        .iter()
        .any(|(face, aabb)| on_face(point, body, *face, aabb, tolerance))
}

pub(crate) fn on_face(
    point: &Point3,
    body: &BRep,
    face: FaceHandle,
    aabb: &Aabb,
    tolerance: &Tolerance,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::validation::FaceUvContainment;
    let linear = tolerance.linear();
    let probe = Aabb {
        min: *point,
        max: *point,
    };
    if !aabb.overlaps(&probe, linear) {
        return false;
    }
    let Some(surface) = body.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let Some((u, v)) =
        popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, point)
    else {
        return false;
    };
    surface.geometry.point_at(u, v).distance_to(point) <= linear
        && body
            .face_uv_trim_containment(face, &surface.geometry, u, v)
            .is_some_and(|c| c != FaceUvContainment::Outside)
}
