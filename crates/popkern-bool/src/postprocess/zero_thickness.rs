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
//! edges whose end points and midpoints coincide, and one vertex
//! shared by separate fans (`pinched_vertex`).
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
    let mut out = Vec::new();

    // Vertices joined by an edge are allowed to be close (a short edge
    // is a feature, not a contact); the set is built at the first close
    // pair, which most results never have.
    let mut joined: Option<BTreeSet<(VertexHandle, VertexHandle)>> = None;
    let joined_pair = |joined: &mut Option<BTreeSet<_>>, pair| {
        joined
            .get_or_insert_with(|| {
                brep.edges
                    .keys()
                    .filter_map(|e| edge_ends(brep, e))
                    .map(|(a, b)| (a.min(b), a.max(b)))
                    .collect()
            })
            .contains(&pair)
    };
    // Close pairs by a sweep over x (each vertex against the later ones
    // within the tolerance in x), then in handle order: the later
    // vertex of each pair names it, as it meets the earlier ones.
    let mut by_x: Vec<(f64, VertexHandle, Point3)> = brep
        .vertices
        .iter()
        .map(|(v, vx)| (vx.position.x, v, vx.position))
        .collect();
    by_x.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut close: Vec<(VertexHandle, VertexHandle)> = Vec::new();
    for (k, &(x, v, p)) in by_x.iter().enumerate() {
        for &(xw, w, q) in &by_x[k + 1..] {
            if xw - x > tol {
                break;
            }
            if q.distance_to(&p) <= tol {
                close.push((v.max(w), v.min(w)));
            }
        }
    }
    close.sort();
    let mut touched = BTreeSet::new();
    for (v, w) in close {
        if !joined_pair(&mut joined, (w, v)) && touched.insert(w) {
            out.push((
                NonManifoldContact {
                    kind: ContactKind::Vertex,
                    location: brep.vertices[w].position,
                },
                vec![ContactEntity::Vertex(w), ContactEntity::Vertex(v)],
            ));
        }
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
    // One vertex shared by separate fans: the contact the result
    // carries without a coincident copy.
    out.extend(super::pinched_vertex::pinched_vertices(brep));
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
