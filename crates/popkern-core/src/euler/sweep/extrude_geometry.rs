//! Geometry authoring for the extrude sweep (PK-316).
//!
//! `BRep::extrude` builds topology with Euler operators; the faces it
//! creates through `mef` INHERIT the base face's surface and its new
//! edges carry no curves at all. Every consumer therefore received a
//! body whose lateral faces all claimed the base plane (measured:
//! validate() red with a 1.41 loop gap on the history sketch→extrude
//! path; the cube/prism constructors survived only because they stomp
//! everything with `attach_geometry` afterwards). This pass authors
//! the geometry the sweep actually created, incrementally — only the
//! new edges, the side faces, and the cap are touched, so extruding a
//! face of a larger model never disturbs unrelated geometry.
//!
//! Authoring rules:
//! - vertical edges: lines from each base vertex to its extruded copy;
//! - top edges: the corresponding base edge's curve translated by the
//!   displacement (a chord line when the base edge carries no curve);
//! - side faces: a plane from the loop winding when the base edge is
//!   straight; a true `Cylinder` when the base edge is a circular arc
//!   whose plane normal is parallel to the displacement (PK-377 — an
//!   extruded arc IS a cylinder, and leaving it as a NURBS degrades
//!   every downstream consumer to approximate paths); a ruled NURBS
//!   (base curve × displacement) otherwise — geometrically exact for
//!   any base curve;
//! - cap face: the base surface translated by the displacement;
//! - orientation flags: derived from the outer-loop winding via a
//!   Newell-normal test at authoring time (the flag is DERIVED state —
//!   PK-244);
//! - pcurves: recomputed per touched face.

use std::collections::{HashMap, HashSet};

use crate::constructors::{line_between, plane_from_points};
use crate::error::CoreError;
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Transform3, Vec3};
use crate::topology::brep::{transform_curve_on_edge, transform_surface};
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::extrude::ExtrudeResult;

/// Author curves and surfaces for the entities an extrude created.
pub(super) fn author_extrude_geometry(
    brep: &mut BRep,
    result: &ExtrudeResult,
    base_face: FaceHandle,
    displacement: Vec3,
) -> Result<(), CoreError> {
    // Capture the base surface FIRST: the input face handle ends the
    // extrude as one of the SIDE faces (see the naming logic in
    // `extrude.rs`), so its surface is replaced by step 2 below.
    let base_surface = brep
        .faces
        .get(base_face)
        .and_then(|f| f.surface.as_ref())
        .map(|sof| sof.geometry.clone());
    let new_vertex_set: HashSet<VertexHandle> = result.new_vertices.iter().copied().collect();
    let new_edge_set: HashSet<EdgeHandle> = result.new_edges.iter().copied().collect();
    let translation = Transform3::translation(displacement.x, displacement.y, displacement.z);

    // Top edges keyed by the base edge they were minted above (recorded
    // by the Euler pass — see `ExtrudeResult::top_edge_of_base`).
    let base_of_top: HashMap<EdgeHandle, EdgeHandle> = result
        .top_edge_of_base
        .iter()
        .map(|&(base, top)| (top, base))
        .collect();

    // 1. New edge curves.
    for &eh in &result.new_edges {
        let he = brep.edges[eh].half_edge;
        let Some((va, vb)) = edge_vertices(brep, he) else {
            continue;
        };
        let pa = brep.vertices[va].position;
        let pb = brep.vertices[vb].position;
        let a_new = new_vertex_set.contains(&va);
        let b_new = new_vertex_set.contains(&vb);
        let curve = if a_new && b_new {
            // Top edge: translate its own base edge's curve.
            let base_curve = base_of_top
                .get(&eh)
                .and_then(|&be| brep.edges.get(be))
                .and_then(|ed| ed.curve.clone());
            match base_curve {
                Some(mut c) => {
                    transform_curve_on_edge(&mut c, &translation, 1.0);
                    align_curve_to_endpoints(&mut c, pa, pb);
                    c
                }
                None => line_between(pa, pb),
            }
        } else {
            // Vertical edge: a line from the base vertex upward.
            line_between(pa, pb)
        };
        brep.edges[eh].curve = Some(curve);
    }

    // 2. Side face surfaces.
    for &side in &result.side_faces {
        let Some(base_he) = base_half_edge_of_side(brep, side, &new_edge_set) else {
            continue;
        };
        let base_edge = brep.half_edges[base_he].edge;
        let base_curve = brep.edges[base_edge].curve.clone();
        let straight = match &base_curve {
            None => true,
            Some(c) => matches!(c.geometry, CurveGeometry::Line { .. }),
        };
        let angular_tol = brep.tolerance.angular();
        let surface = if straight {
            side_plane_from_loop(brep, side)?
        } else {
            match base_curve
                .as_ref()
                .and_then(|c| analytic_side_surface(c, displacement, angular_tol))
                .or_else(|| ruled_side_surface(base_curve.as_ref(), displacement))
            {
                Some(s) => s,
                None => side_plane_from_loop(brep, side)?,
            }
        };
        let orientation = winding_orientation(brep, side, &surface);
        brep.faces[side].surface = Some(SurfaceOnFace {
            geometry: surface,
            orientation,
        });
    }

    // 3. Cap face surface: the base surface translated.
    let cap_surface = match base_surface.as_ref() {
        Some(geometry) => transform_surface(geometry, &translation, 1.0),
        None => side_plane_from_loop(brep, result.cap_face)?,
    };
    let cap_orientation = winding_orientation(brep, result.cap_face, &cap_surface);
    brep.faces[result.cap_face].surface = Some(SurfaceOnFace {
        geometry: cap_surface,
        orientation: cap_orientation,
    });

    // 4. Pcurves for the touched faces.
    let mut touched = result.side_faces.clone();
    touched.push(result.cap_face);
    for f in touched {
        brep.recompute_face_pcurves(f)?;
    }
    Ok(())
}

/// Both endpoint vertices of an edge via its stored half-edge.
fn edge_vertices(
    brep: &BRep,
    he: crate::topology::HalfEdgeHandle,
) -> Option<(VertexHandle, VertexHandle)> {
    let hed = brep.half_edges.get(he)?;
    let twin = brep.half_edges.get(hed.twin)?;
    Some((hed.origin, twin.origin))
}

/// Reverse a translated curve when its parameter direction opposes
/// the stored half-edge's vertex order (the translated base curve may
/// run either way relative to the top edge's canonical direction).
fn align_curve_to_endpoints(curve: &mut CurveOnEdge, start: Point3, end: Point3) {
    let c_start = curve.geometry.point_at(curve.t_start);
    if c_start.distance_to(&start) <= c_start.distance_to(&end) {
        return;
    }
    std::mem::swap(&mut curve.t_start, &mut curve.t_end);
}

/// The side face's base half-edge: the one whose edge existed before
/// the extrude.
fn base_half_edge_of_side(
    brep: &BRep,
    side: FaceHandle,
    new_edges: &HashSet<EdgeHandle>,
) -> Option<crate::topology::HalfEdgeHandle> {
    let start = brep.loops.get(brep.faces.get(side)?.outer_loop)?.half_edge;
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let hed = brep.half_edges.get(cur)?;
        if !new_edges.contains(&hed.edge) {
            return Some(cur);
        }
        cur = hed.next;
        if cur == start {
            break;
        }
    }
    None
}

/// A plane from the first three loop vertices, oriented by the loop
/// winding (the `attach_geometry` convention).
fn side_plane_from_loop(brep: &BRep, face: FaceHandle) -> Result<SurfaceGeometry, CoreError> {
    let start = brep
        .loops
        .get(
            brep.faces
                .get(face)
                .ok_or_else(|| CoreError::InvalidHandle("side face".into()))?
                .outer_loop,
        )
        .ok_or_else(|| CoreError::InvalidHandle("side loop".into()))?
        .half_edge;
    let he0 = &brep.half_edges[start];
    let he1 = &brep.half_edges[he0.next];
    let he2 = &brep.half_edges[he1.next];
    plane_from_points(
        brep.vertices[he0.origin].position,
        brep.vertices[he1.origin].position,
        brep.vertices[he2.origin].position,
    )
}

/// Analytic promotion of the side surface (PK-377): a circular-arc base
/// edge extruded along its own plane normal sweeps a cylinder of the
/// arc's centre, axis, and radius. The axis keeps the ARC's normal
/// sign (not the displacement's) so the chart's angular direction
/// agrees with the arc's own parameter sense. `angular_tol` is the
/// body's explicit angular tolerance: the sine of the angle between
/// the arc normal and the displacement must not exceed it. Every
/// other curved edge, and any arc extruded obliquely (an elliptic
/// cylinder), returns `None` and keeps the exact ruled NURBS.
fn analytic_side_surface(
    base: &CurveOnEdge,
    displacement: Vec3,
    angular_tol: f64,
) -> Option<SurfaceGeometry> {
    let CurveGeometry::Arc {
        center,
        normal,
        radius,
        ..
    } = &base.geometry
    else {
        return None;
    };
    let axis = normal.normalized()?;
    let dir = displacement.normalized()?;
    if axis.cross(&dir).length() > angular_tol || !radius.is_finite() || *radius <= 0.0 {
        return None;
    }
    Some(SurfaceGeometry::Cylinder {
        origin: *center,
        axis,
        radius: *radius,
    })
}

/// Ruled surface `base_curve × displacement` as a degree (p,1) NURBS —
/// geometrically exact for any base curve (the ruling is a rigid
/// translation of every control point).
fn ruled_side_surface(base: Option<&CurveOnEdge>, displacement: Vec3) -> Option<SurfaceGeometry> {
    let base = base?;
    let nurbs = base.geometry.to_nurbs_in_range((base.t_start, base.t_end));
    let n = nurbs.control_points.len();
    let mut control_points = Vec::with_capacity(2 * n);
    let mut weights = Vec::with_capacity(2 * n);
    // Storage is u-major (index = i_u * num_v + i_v) with v as the
    // ruling direction.
    for (cp, w) in nurbs.control_points.iter().zip(nurbs.weights.iter()) {
        control_points.push(*cp);
        control_points.push(*cp + displacement);
        weights.push(*w);
        weights.push(*w);
    }
    Some(SurfaceGeometry::Nurbs(
        crate::geometry::nurbs_surface::NurbsSurface {
            degree_u: nurbs.degree,
            degree_v: 1,
            knots_u: nurbs.knots.clone(),
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            num_u: n,
            num_v: 2,
            control_points,
            weights,
        },
    ))
}

/// Orientation flag from the outer-loop winding: Newell normal of the
/// loop's SAMPLED outline against the surface normal at the outline's
/// centroid. Every half-edge's curve is sampled in traversal
/// direction — the vertex polygon alone is degenerate for a curved
/// wall (a half-cylinder's four corners span the chord plane, whose
/// normal is perpendicular to the cylinder normal at those corners:
/// the flag was a coin toss on every arc-sided extrusion, PK-377) —
/// and the centroid lands mid-wall, where the surface normal is
/// unambiguous.
fn winding_orientation(brep: &BRep, face: FaceHandle, surface: &SurfaceGeometry) -> bool {
    const SAMPLES: usize = 16;
    let Some(fd) = brep.faces.get(face) else {
        return true;
    };
    let Some(ld) = brep.loops.get(fd.outer_loop) else {
        return true;
    };
    let start = ld.half_edge;
    let mut pts: Vec<Point3> = Vec::new();
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let Some(hed) = brep.half_edges.get(cur) else {
            break;
        };
        let origin = brep.vertices[hed.origin].position;
        match brep.edges.get(hed.edge).and_then(|e| e.curve.as_ref()) {
            Some(c) => {
                // Sample strictly inside the curve's own range: a
                // malformed edge range must degrade to its vertex,
                // never panic.
                let (r0, r1) = c.geometry.parameter_range();
                let (lo, hi) = (r0.min(r1), r0.max(r1));
                let clamp = |t: f64| {
                    if lo.is_finite() && hi.is_finite() {
                        t.clamp(lo, hi)
                    } else {
                        t
                    }
                };
                let (ta, tb) = (clamp(c.t_start), clamp(c.t_end));
                let forward = c.geometry.point_at(ta).distance_to(&origin)
                    <= c.geometry.point_at(tb).distance_to(&origin);
                for i in 0..SAMPLES {
                    let f = i as f64 / SAMPLES as f64;
                    let t = if forward {
                        ta + (tb - ta) * f
                    } else {
                        tb + (ta - tb) * f
                    };
                    pts.push(c.geometry.point_at(t));
                }
            }
            None => pts.push(origin),
        }
        cur = hed.next;
        if cur == start {
            break;
        }
    }
    if pts.len() < 3 {
        return true;
    }
    let mut newell = Vec3::new(0.0, 0.0, 0.0);
    let mut centroid = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..pts.len() {
        let a = pts[i];
        let b = pts[(i + 1) % pts.len()];
        newell.x += (a.y - b.y) * (a.z + b.z);
        newell.y += (a.z - b.z) * (a.x + b.x);
        newell.z += (a.x - b.x) * (a.y + b.y);
        centroid = centroid + Vec3::new(a.x, a.y, a.z);
    }
    let n = pts.len() as f64;
    let centroid = Point3::new(centroid.x / n, centroid.y / n, centroid.z / n);
    let Some((u, v)) = crate::geometry::uv_projection::project_point_to_uv(surface, &centroid)
    else {
        return true;
    };
    surface.normal_at(u, v).dot(&newell) >= 0.0
}
