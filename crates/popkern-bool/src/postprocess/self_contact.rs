//! A shell touching itself: two of its faces that share no vertex
//! meeting at a point or along a curve (PK-1155).
//!
//! The shell-contact search compares faces of different shells, and the
//! copy tests read vertices and edges in one place: a result that
//! touches itself inside ONE shell, where nothing of its topology
//! stands at the touch, passed both. Measured: a box united with a body
//! that overlaps it at its foot and, elsewhere, rests a ball on the
//! box's top edge (or on its top face) served one valid solid touching
//! itself at the point, under every policy.
//!
//! In a manifold shell two faces meet only through a vertex or an edge
//! they share, so faces sharing no vertex that come within the
//! tolerance are a contact. Pairs are pruned before the search, which
//! is the shell search's (seeds refined by alternating projection, the
//! signed gap's critical points, each face's edges walked against the
//! other): the boxes must overlap (each holds its face, save the grid
//! box of a face with neither a closed form nor patches), the faces
//! share no vertex, and no bound may prove them apart. One face's
//! box farther than the tolerance from the other's surface clears a
//! pair (`ff_parallel_reach`); a NURBS face is held by the Bézier
//! patches of the knot spans its trim reaches, cleared against the
//! other face by separating axes and halving the patches still near
//! (`self_contact_patches.rs`: the cycloidal disc's walls, two flanks a
//! crest 0.05 wide apart, cost the search seconds a pair); and faces
//! that do not cross meet nearest on a boundary of one unless their
//! normals are parallel there, so two planes, a plane and a cylinder
//! whose normals are never parallel within the tolerance, or two such
//! cylinders, are cleared by their edges alone (`self_contact_clear.rs`:
//! each edge farther than the tolerance from the other face, its plane
//! inside its box or outside the face itself, its box in closed form).
//! Two planar faces bounded by segments are cleared by their exact gap
//! (`self_contact_polygon.rs`: a faceted die's triangle pairs cost the
//! edge test half a millisecond each), and two faces on one surface of
//! reach well over the tolerance by their edges alone, since they meet
//! only where their boundaries do (`self_contact_same_surface.rs`: a
//! torus's pieces, whose boxes lie on each other's surface). Faces
//! meeting near a vertex they share elsewhere are not searched.

use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::error::{ContactKind, NonManifoldContact};

/// Halvings of an edge curve tested against a NURBS face's patches
/// before its walk is left to run.
const CURVE_DEPTH: usize = 8;
use clear::{
    cylinders_clear, edges_clear_of_plane, plane_crosses_cylinder_clear, plane_of, straight,
};

/// Every place where two faces of one shell sharing no vertex touch, in
/// face-handle order, one per location.
pub(crate) fn self_contacts(brep: &BRep, tolerance: &Tolerance) -> Vec<NonManifoldContact> {
    let linear = tolerance.linear_reach();
    // Each face with a box holding it and, on a NURBS surface, the
    // Bézier patches holding it (`self_contact_patches.rs`), whose hull
    // is its box.
    let mut faces: Vec<(FaceHandle, Aabb, Option<Vec<patches::Patch>>)> = brep
        .faces
        .keys()
        .filter_map(|f| {
            let held = match &brep.faces.get(f)?.surface.as_ref()?.geometry {
                SurfaceGeometry::Nurbs(n) => patches::face_patches(brep, f, n),
                _ => None,
            };
            let b = match held.as_deref().and_then(patches::hull) {
                Some(b) => b,
                None => containing_box(brep, f)?,
            };
            Some((f, b, held))
        })
        .collect();
    faces.sort_by_key(|(f, ..)| *f);
    // Each face's boundary (vertices and edges sorted), walked only for
    // a face some pair's boxes reach; vertex loops gathered in one pass.
    let mut lone: Vec<(FaceHandle, VertexHandle)> = brep
        .loops
        .values()
        .filter_map(|l| l.vertex.map(|v| (l.face, v)))
        .collect();
    lone.sort_unstable();
    let mut bounds: Vec<Option<Boundary>> = vec![None; faces.len()];
    let mut polygons: Vec<Option<Option<polygon::Polygon>>> = Vec::new();
    polygons.resize_with(faces.len(), || None);
    let surface = |f: FaceHandle| brep.faces.get(f).and_then(|x| x.surface.as_ref());
    let shell = |f: FaceHandle| brep.faces.get(f).map(|x| x.back);
    let mut out: Vec<NonManifoldContact> = Vec::new();
    for i in 0..faces.len() {
        let (f, box_f, patches_f) = (&faces[i].0, &faces[i].1, &faces[i].2);
        let Some(sf) = surface(*f) else { continue };
        for j in i + 1..faces.len() {
            let (g, box_g, patches_g) = (&faces[j].0, &faces[j].1, &faces[j].2);
            if shell(*f) != shell(*g) || !box_f.overlaps(box_g, linear) {
                continue;
            }
            for k in [i, j] {
                if bounds[k].is_none() {
                    bounds[k] = Some(face_boundary(brep, faces[k].0, &lone));
                }
            }
            let (Some(bf), Some(bg)) = (&bounds[i], &bounds[j]) else {
                continue;
            };
            let Some(sg) = surface(*g) else { continue };
            if shares_any(&bf.vertices, &bg.vertices) {
                // Faces meeting only at a vertex where they face each
                // other (normals opposed) can touch beside it: a ball
                // resting on a box's top 1e-3 from where a rod's rim,
                // the rod fused to the ball, met the ball's own rim
                // under the tolerance served touching itself (PK-1177).
                // The touch counts farther than ten tolerances from the
                // shared vertices; nearer, it is theirs.
                if shares_any(&bf.edges, &bg.edges) {
                    continue;
                }
                let shared: Vec<Point3> = bf
                    .vertices
                    .iter()
                    .filter(|v| bg.vertices.binary_search(v).is_ok())
                    .filter_map(|&v| brep.vertices.get(v).map(|x| x.position))
                    .collect();
                if !shared.iter().any(|&p| facing_at(sf, sg, p)) {
                    continue;
                }
                if let Some(at) = super::shell_contact::face_contact_avoiding(
                    brep,
                    (*f, &sf.geometry),
                    (*g, &sg.geometry),
                    linear,
                    (&shared, 10.0 * linear),
                ) && out
                    .iter()
                    .all(|c| c.location.distance_to(&at) > 10.0 * linear)
                {
                    out.push(NonManifoldContact {
                        kind: ContactKind::Vertex,
                        location: at,
                    });
                }
                continue;
            }
            if crate::parallel::ff_parallel::ff_parallel_reach::faces_apart(
                (&sf.geometry, box_f),
                (&sg.geometry, box_g),
                linear,
            ) {
                continue;
            }
            // Two faces on one surface meet only where their boundaries
            // do (`self_contact_same_surface.rs`).
            if sf.geometry == sg.geometry
                && same_surface::same_surface_clear(
                    brep,
                    &sf.geometry,
                    (&bf.edges, &bg.edges),
                    linear,
                )
            {
                continue;
            }
            // A NURBS face: the patches holding it, each farther than the
            // band from the other face (`self_contact_patches.rs`).
            let off_f = |c: Point3, r: f64| clear::ball_off_face(brep, *f, c, r);
            let off_g = |c: Point3, r: f64| clear::ball_off_face(brep, *g, c, r);
            let planar = |s: &SurfaceGeometry| plane_of(s).is_some();
            let (off_f, off_g): (Option<&patches::OffFace>, Option<&patches::OffFace>) = (
                planar(&sf.geometry).then_some(&off_f as &patches::OffFace),
                planar(&sg.geometry).then_some(&off_g as &patches::OffFace),
            );
            let apart = match (patches_f, patches_g) {
                (Some(a), Some(b)) => {
                    patches::patches_apart(a, &patches::Other::Patches(b), linear)
                }
                (Some(a), None) => patches::patches_apart(
                    a,
                    &patches::Other::Face(&sg.geometry, box_g, off_g),
                    linear,
                ),
                (None, Some(b)) => patches::patches_apart(
                    b,
                    &patches::Other::Face(&sf.geometry, box_f, off_f),
                    linear,
                ),
                (None, None) => false,
            };
            if apart {
                continue;
            }
            // A planar face and a cylindrical one whose normals are never
            // parallel (the axis not in the plane's directions): faces of
            // a valid shell do not cross, so their nearest points lie on
            // a boundary of one, else on both interiors with the segment
            // between normal to both. Cleared when no edge of the
            // cylinder's face comes within the band of the plane inside
            // the planar face's box, nor an edge of the planar face
            // within the band of the cylinder's face's box.
            let (side_f, side_g) = (
                (*f, &sf.geometry, box_f, bf.edges.as_slice()),
                (*g, &sg.geometry, box_g, bg.edges.as_slice()),
            );
            if plane_crosses_cylinder_clear(brep, side_f, side_g, linear)
                || plane_crosses_cylinder_clear(brep, side_g, side_f, linear)
                || cylinders_clear(brep, side_f, side_g, linear)
            {
                continue;
            }
            // Planar faces touch only where a boundary of one meets the
            // other: no edge of either near the other's plane within its
            // box clears the pair in closed form.
            // Two faces bounded by segments: their exact gap
            // (`self_contact_polygon.rs`).
            if let (Some(pf), Some(pg)) = (plane_of(&sf.geometry), plane_of(&sg.geometry)) {
                for k in [i, j] {
                    if polygons[k].is_none() {
                        polygons[k] = Some(polygon::polygon_of(brep, faces[k].0));
                    }
                }
                let cleared = match (&polygons[i], &polygons[j]) {
                    (Some(Some(a)), Some(Some(b))) => polygon::polygons_apart(a, b, linear),
                    _ => {
                        edges_clear_of_plane(brep, &bf.edges, pg, (*g, box_g), linear)
                            && edges_clear_of_plane(brep, &bg.edges, pf, (*f, box_f), linear)
                    }
                };
                if cleared {
                    continue;
                }
            }
            // Each NURBS face read on its own piece of its surface.
            let found = {
                let piece = |h: FaceHandle, sh: &SurfaceGeometry| match sh {
                    SurfaceGeometry::Nurbs(n) => {
                        patches::face_piece(brep, h, n).map(SurfaceGeometry::Nurbs)
                    }
                    _ => None,
                };
                let (pf, pg) = (piece(*f, &sf.geometry), piece(*g, &sg.geometry));
                // An edge walked against the other face is skipped where
                // its control hull, halved while near, is apart from that
                // face's patches or its analytic surface and box: a slot
                // wall walked against a thread flank beside it cost 1.5 s
                // a pair (PK-1172).
                let clear = |walked: FaceHandle, c: &CurveOnEdge| {
                    let (other, surface, bx) = if walked == *f {
                        (patches_g.as_deref(), &sg.geometry, box_g)
                    } else {
                        (patches_f.as_deref(), &sf.geometry, box_f)
                    };
                    let curve = c.geometry.to_nurbs_in_range((c.t_start, c.t_end));
                    match other {
                        Some(theirs) => patches::curve_apart(&curve, theirs, linear, CURVE_DEPTH),
                        None => crate::parallel::ff_parallel::ff_parallel_reach::curve_clear_of_surface_face(
                            &curve, surface, bx, linear,
                        ),
                    }
                };
                super::shell_contact::face_contact_on(
                    brep,
                    (*f, pf.as_ref().unwrap_or(&sf.geometry)),
                    (*g, pg.as_ref().unwrap_or(&sg.geometry)),
                    linear,
                    Some(&clear),
                )
            };
            let Some(at) = found else {
                continue;
            };
            if out
                .iter()
                .all(|c| c.location.distance_to(&at) > 10.0 * linear)
            {
                out.push(NonManifoldContact {
                    kind: ContactKind::Vertex,
                    location: at,
                });
            }
        }
    }
    out
}

/// A box holding `face`, not trim-scale where that costs: a sphere's or
/// torus's whole surface (exact for the whole surface) (the trim-scale box searches the surface's
/// critical points inside the trim, most of this test's cost on a
/// filleted body), the boolean's face box otherwise. A looser box only
/// adds pairs for the closed-form bounds to clear.
fn containing_box(brep: &BRep, face: FaceHandle) -> Option<Aabb> {
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let around = |c: Point3, r: f64| {
        let d = Vec3::new(r, r, r);
        Aabb {
            min: c - d,
            max: c + d,
        }
    };
    match surface {
        SurfaceGeometry::Sphere { center, radius, .. } => Some(around(*center, radius.abs())),
        // The whole torus: along a coordinate axis it reaches its tube
        // circle's radius R times the sine of that axis's angle to its
        // own, plus the tube's r.
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let a = axis.normalized()?;
            let (big, small) = (major_radius.abs(), minor_radius.abs());
            let reach = |c: f64| big * (1.0 - c * c).max(0.0).sqrt() + small;
            let d = Vec3::new(reach(a.x), reach(a.y), reach(a.z));
            Some(Aabb {
                min: *center - d,
                max: *center + d,
            })
        }
        // Their faces lie within their boundaries' extent (a plane's
        // trivially; a cylinder's or cone's extreme along any direction
        // lies on a ruling the boundary crosses): what the face box takes
        // too, without its loop samples.
        SurfaceGeometry::Plane { .. }
        | SurfaceGeometry::Cylinder { .. }
        | SurfaceGeometry::Cone { .. } => boundary_extent(brep, face),
        _ => crate::intersect::face_aabb(brep, face),
    }
}

/// The box of `face`'s boundary curves: each edge's ends and its
/// conservative extent.
fn boundary_extent(brep: &BRep, face: FaceHandle) -> Option<Aabb> {
    let f = brep.faces.get(face)?;
    let mut points = Vec::new();
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let start = brep.loops.get(l)?.half_edge;
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            points.push(brep.vertices.get(h.origin)?.position);
            // A straight edge lies between its vertices.
            if let Some(c) = brep.edges.get(h.edge).and_then(|e| e.curve.as_ref())
                && !straight(c)
            {
                crate::face_box::conservative_curve_extent(c, &mut points);
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    Aabb::from_points(&points)
}

/// A face's vertices and edges, each sorted.
#[derive(Clone, Default)]
struct Boundary {
    vertices: Vec<VertexHandle>,
    edges: Vec<EdgeHandle>,
}

/// `face`'s boundary: its loops' vertices and edges and the vertices of
/// its vertex loops (`lone`, sorted by face).
fn face_boundary(brep: &BRep, face: FaceHandle, lone: &[(FaceHandle, VertexHandle)]) -> Boundary {
    let mut out = Boundary::default();
    let Some(f) = brep.faces.get(face) else {
        return out;
    };
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(l).map(|x| x.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            out.vertices.push(h.origin);
            out.edges.push(h.edge);
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    let from = lone.partition_point(|(g, _)| *g < face);
    out.vertices.extend(
        lone[from..]
            .iter()
            .take_while(|(g, _)| *g == face)
            .map(|(_, v)| *v),
    );
    out.vertices.sort_unstable();
    out.vertices.dedup();
    out.edges.sort_unstable();
    out.edges.dedup();
    out
}

/// Whether the faces on `a` and `b` face each other at `p`: their
/// outward normals there opposed within a hundredth of a radian's
/// cosine (a tangent pinch, the faces' interiors on either side of one
/// tangent plane).
fn facing_at(
    a: &popkern_core::geometry::surface::SurfaceOnFace,
    b: &popkern_core::geometry::surface::SurfaceOnFace,
    p: Point3,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    let outward = |s: &popkern_core::geometry::surface::SurfaceOnFace| {
        let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(&s.geometry, &p)?;
        let n = s.geometry.normal_at(u, v);
        Some(if s.orientation { n } else { n * -1.0 })
    };
    match (outward(a), outward(b)) {
        (Some(x), Some(y)) => x.dot(&y) < -0.99995,
        _ => false,
    }
}

/// Whether two sorted lists share an element.
fn shares_any<T: Ord>(a: &[T], b: &[T]) -> bool {
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => return true,
        }
    }
    false
}

#[path = "self_contact_clear.rs"]
mod clear;
#[path = "self_contact_patches.rs"]
mod patches;
#[path = "self_contact_polygon.rs"]
mod polygon;
#[path = "self_contact_same_surface.rs"]
mod same_surface;

#[cfg(test)]
#[path = "self_contact_tests.rs"]
mod tests;
