//! The sheet primitives (PK-1083): one planar face with free edges
//! (the canonical sheet: each boundary edge's other half-edge null and
//! the edge in `boundary_edges`, as `remove_faces` leaves them), so a
//! sheet thickens, offsets and cuts like one stripped from a solid.
//!
//! Built through the exact profile path (arcs stay arcs, splines stay
//! splines, pcurves authored by the sketch compiler): the profile's
//! face comes as a front face and a surfaceless back face sharing its
//! edges; the back face is retired to the boundary. A fill takes a
//! closed planar chain of wire edges, in any order and direction, and
//! spells it as a profile on the chain's own plane. Names: the face
//! `Cap 0`, its edges `Intersection i` and vertices `Junction i` in
//! loop order, each with the face as ancestor.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
use popkern_core::naming_criteria::{
    GeometricCriteria, capture_edge, capture_face, capture_vertex,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::plane::SketchPlane;
use popkern_sketch::profile::{AnalyticProfile, PathSegment, ProfileShape};

use super::{failed, invalid};
use crate::error::ReplayError;
use crate::working_bodies::WorkingBodies;

/// The profile's face on `plane` as a one-face sheet.
pub(super) fn profile_sheet(
    shape: &str,
    profile: &AnalyticProfile<f64>,
    plane: SketchPlane,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let facing = plane.normal;
    let mut brep = profile
        .to_brep(plane, tolerance, None)
        .map_err(|e| failed(shape, "profile", e))?
        .brep;
    // The back face points against the plane's normal (or has no
    // surface, the core constructors' form).
    let back: Vec<_> = brep
        .faces
        .iter()
        .filter(|(_, f)| match f.surface.as_ref() {
            None => true,
            Some(s) => match s.geometry {
                popkern_core::geometry::surface::SurfaceGeometry::Plane { normal, .. } => {
                    let n = if s.orientation { normal } else { -normal };
                    n.dot(&facing) < 0.0
                }
                _ => false,
            },
        })
        .map(|(h, _)| h)
        .collect();
    if back.len() != 1 || brep.faces.len() != 2 {
        return Err(failed(
            shape,
            "sheet",
            format!(
                "the profile's face came as {} faces, {} facing back; a front and a back were expected",
                brep.faces.len(),
                back.len()
            ),
        ));
    }
    brep.retire_faces_to_boundary(&back)
        .map_err(|e| failed(shape, "sheet", e))?;
    Ok(brep)
}

/// The face a closed planar chain of wire edges bounds.
pub(super) fn fill(
    shape: &str,
    edges: &[EntityName],
    bodies: &WorkingBodies,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let join = tolerance.linear() * 10.0;
    let mut ids: Vec<_> = bodies.keys().copied().collect();
    ids.sort_by_key(|b| b.0);
    let mut curves = Vec::with_capacity(edges.len());
    for (i, name) in edges.iter().enumerate() {
        let found = ids.iter().filter_map(|id| bodies.get(id)).find_map(|b| {
            crate::resolve::resolve_edge(name, b)
                .ok()
                .and_then(|e| b.edges[e].curve.clone())
        });
        let curve = found.ok_or_else(|| invalid(shape, "edges", i as f64, "edge not found"))?;
        curves.push(curve);
    }
    let chain = crate::curve_chain::chain(curves, join).map_err(|e| match e {
        crate::curve_chain::ChainError::Empty => {
            invalid(shape, "edges", 0.0, "needs at least one edge")
        }
        crate::curve_chain::ChainError::Broken { chained } => invalid(
            shape,
            "edges",
            chained as f64,
            "the edges do not chain end to end",
        ),
    })?;
    let (head, _) = chain[0].ends();
    let (_, tail) = chain[chain.len() - 1].ends();
    if head.distance_to(&tail) > join {
        return Err(invalid(shape, "edges", 0.0, "the edges do not close"));
    }
    // The chain's plane: Newell over samples, every sample on it.
    let samples: Vec<Point3> = chain
        .iter()
        .flat_map(|p| (0..16).map(move |k| p.at(k as f64 / 16.0)))
        .collect();
    let mut normal = Vec3::new(0.0, 0.0, 0.0);
    for (i, p) in samples.iter().enumerate() {
        let q = samples[(i + 1) % samples.len()];
        normal = normal
            + Vec3::new(
                (p.y - q.y) * (p.z + q.z),
                (p.z - q.z) * (p.x + q.x),
                (p.x - q.x) * (p.y + q.y),
            );
    }
    let normal = normal
        .normalized()
        .ok_or_else(|| invalid(shape, "edges", 0.0, "the chain encloses no area"))?;
    let origin = samples[0];
    if let Some(off) = samples
        .iter()
        .map(|p| (*p - origin).dot(&normal).abs())
        .find(|d| *d > tolerance.linear())
    {
        return Err(invalid(
            shape,
            "edges",
            off,
            "the edges do not lie in one plane (a fill needs a planar wire)",
        ));
    }
    let plane = SketchPlane::from_normal(origin, normal).map_err(|e| failed(shape, "plane", e))?;
    let uv = |p: Point3| {
        let d = p - plane.origin;
        [d.dot(&plane.u_axis), d.dot(&plane.v_axis)]
    };
    let profile = match chain.as_slice() {
        [only]
            if matches!(only.curve.geometry, CurveGeometry::Circle { .. })
                || (matches!(only.curve.geometry, CurveGeometry::Arc { .. })
                    && head.distance_to(&tail) <= join) =>
        {
            let (CurveGeometry::Circle { center, radius, .. }
            | CurveGeometry::Arc { center, radius, .. }) = only.curve.geometry
            else {
                return Err(invalid(shape, "edges", 0.0, "not a circle"));
            };
            AnalyticProfile {
                shape: ProfileShape::Circle { radius },
                center: uv(center),
            }
        }
        _ => {
            let segments = chain.iter().map(|p| segment(p, &uv)).collect();
            AnalyticProfile::new(ProfileShape::Path {
                start: uv(head),
                segments,
            })
        }
    };
    profile_sheet(shape, &profile, plane, tolerance)
}

/// One chain piece as a profile segment on the plane.
fn segment(
    piece: &crate::curve_chain::ChainPiece,
    uv: &dyn Fn(Point3) -> [f64; 2],
) -> PathSegment<f64> {
    let (_, to) = piece.ends();
    match piece.curve.geometry {
        CurveGeometry::Line { .. } => PathSegment::Line { to: uv(to) },
        CurveGeometry::Arc { .. } | CurveGeometry::Circle { .. } => PathSegment::Arc {
            via: uv(piece.at(0.5)),
            to: uv(to),
        },
        _ => {
            let nurbs = piece.nurbs();
            let n = nurbs.control_points.len();
            PathSegment::Spline {
                control: nurbs.control_points[1..n - 1]
                    .iter()
                    .map(|p| uv(*p))
                    .collect(),
                to: uv(to),
                degree: nurbs.degree,
                weights: nurbs.weights.clone(),
                knots: nurbs.knots.clone(),
            }
        }
    }
}

/// Name the sheet's face, then its edges and vertices in loop order.
pub(super) fn name_sheet(brep: &mut BRep, ctx: &mut TrackingContext) {
    let mint = |ctx: &mut TrackingContext, category, index: usize, ancestors| {
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
    };
    let Some(face) = brep.faces.keys().next() else {
        return;
    };
    let name = mint(ctx, RoleCategory::Cap, 0, Vec::new());
    let criteria = capture_face(brep, face).map(GeometricCriteria::Face);
    let face_id = brep.name_pool.intern(name, criteria);
    brep.faces[face].name = Some(face_id);
    let mut loop_edges = Vec::new();
    let mut loop_vertices = Vec::new();
    let loops: Vec<_> = std::iter::once(brep.faces[face].outer_loop)
        .chain(brep.faces[face].inner_loops.iter().copied())
        .collect();
    for l in loops {
        let start = brep.loops[l].half_edge;
        let mut he = start;
        loop {
            let h = &brep.half_edges[he];
            if !loop_edges.contains(&h.edge) {
                loop_edges.push(h.edge);
            }
            if !loop_vertices.contains(&h.origin) {
                loop_vertices.push(h.origin);
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    for (i, eh) in loop_edges.into_iter().enumerate() {
        let name = mint(ctx, RoleCategory::Intersection, i, vec![face_id]);
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let id = brep.name_pool.intern(name, criteria);
        brep.edges[eh].name = Some(id);
    }
    for (i, vh) in loop_vertices.into_iter().enumerate() {
        let name = mint(ctx, RoleCategory::Junction, i, vec![face_id]);
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let id = brep.name_pool.intern(name, criteria);
        brep.vertices[vh].name = Some(id);
    }
}
