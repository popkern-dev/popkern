//! Topology queries and small geometry helpers for the fillet overflow
//! splice (`strip_overflow_surgery`, PK-1189): faces by surface, loops,
//! neighbours, edges through points, and the edge splits and new edges
//! the splice makes.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::strip_overflow::OverflowPlan;

/// Whether `face` has an edge shared with `other` passing within `reach`
/// of `p`.
pub(super) fn edge_through(
    b: &BRep,
    face: FaceHandle,
    other: FaceHandle,
    p: &Point3,
    reach: f64,
) -> bool {
    face_half_edges(b, face).into_iter().any(|h| {
        mate_face(b, h) == Some(other)
            && b.edges
                .get(b.half_edges[h].edge)
                .and_then(|e| e.curve.as_ref())
                .is_some_and(|c| {
                    let tt = popkern_core::geometry::curve::closest_point_in_window(
                        &c.geometry,
                        p,
                        c.t_start,
                        c.t_end,
                    );
                    c.geometry.point_at(tt).distance_to(p) <= reach
                })
    })
}

pub(super) fn surface_of(b: &BRep, face: FaceHandle) -> Option<&SurfaceGeometry> {
    b.faces.get(face)?.surface.as_ref().map(|s| &s.geometry)
}

pub(super) fn faces_with_surface(b: &BRep, surface: &SurfaceGeometry) -> Vec<FaceHandle> {
    b.faces
        .iter()
        .filter(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| same_surface(&s.geometry, surface, SAME_SURFACE_REACH))
        })
        .map(|(h, _)| h)
        .collect()
}

/// How far apart two descriptions of one surface may sit (the apply
/// phase re-anchors a strip's cylinder along its axis and may flip it).
pub(super) const SAME_SURFACE_REACH: f64 = 1e-9;

/// Whether two surfaces are one: a cylinder or torus with the same axis
/// line either way and the same radii, a plane with the same normal line
/// and offset, anything else by value.
pub(super) fn same_surface(a: &SurfaceGeometry, c: &SurfaceGeometry, reach: f64) -> bool {
    let reach = reach.max(SAME_SURFACE_REACH);
    let parallel = |u: &popkern_core::geometry::Vec3, v: &popkern_core::geometry::Vec3| {
        u.cross(v).length() <= reach
    };
    let on_line = |p: &Point3, o: &Point3, axis: &popkern_core::geometry::Vec3| {
        let d = *p - *o;
        (d - *axis * d.dot(axis)).length() <= reach
    };
    match (a, c) {
        (
            SurfaceGeometry::Cylinder {
                origin: o1,
                axis: a1,
                radius: r1,
            },
            SurfaceGeometry::Cylinder {
                origin: o2,
                axis: a2,
                radius: r2,
            },
        ) => (r1 - r2).abs() <= reach && parallel(a1, a2) && on_line(o2, o1, a1),
        (
            SurfaceGeometry::Torus {
                center: c1,
                axis: a1,
                major_radius: m1,
                minor_radius: n1,
                ..
            },
            SurfaceGeometry::Torus {
                center: c2,
                axis: a2,
                major_radius: m2,
                minor_radius: n2,
                ..
            },
        ) => {
            c1.distance_to(c2) <= reach
                && parallel(a1, a2)
                && (m1 - m2).abs() <= reach
                && (n1 - n2).abs() <= reach
        }
        (
            SurfaceGeometry::Plane {
                origin: o1,
                normal: n1,
                ..
            },
            SurfaceGeometry::Plane {
                origin: o2,
                normal: n2,
                ..
            },
        ) => parallel(n1, n2) && (*o2 - *o1).dot(n1).abs() <= reach,
        _ => a == c,
    }
}

/// The shared support face in the applied body: the planned handle when
/// it survives, else (a splice before this one merged it away) the face
/// on the support surface with an edge through P1.
pub(super) fn support_face(
    b: &BRep,
    plan: &OverflowPlan,
    reach: f64,
) -> Result<FaceHandle, String> {
    if b.faces.contains_key(plan.face_f) {
        return Ok(plan.face_f);
    }
    faces_with_surface(b, &plan.support_surface)
        .into_iter()
        .find(|&face| {
            face_edges(b, face).into_iter().any(|e| {
                b.edges
                    .get(e)
                    .and_then(|e| e.curve.as_ref())
                    .is_some_and(|c| {
                        let tt = popkern_core::geometry::curve::closest_point_in_window(
                            &c.geometry,
                            &plan.p1,
                            c.t_start,
                            c.t_end,
                        );
                        c.geometry.point_at(tt).distance_to(&plan.p1) <= reach
                    })
            })
        })
        .ok_or_else(|| "the shared support face did not survive the apply phase".into())
}

pub(super) fn loop_half_edges(
    b: &BRep,
    l: popkern_core::topology::LoopHandle,
) -> Vec<HalfEdgeHandle> {
    let mut out = Vec::new();
    let Some(lp) = b.loops.get(l) else {
        return out;
    };
    let start = lp.half_edge;
    let mut h = start;
    loop {
        out.push(h);
        let Some(he) = b.half_edges.get(h) else {
            break;
        };
        h = he.next;
        if h == start || out.len() > 100_000 {
            break;
        }
    }
    out
}

pub(super) fn face_half_edges(b: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let Some(f) = b.faces.get(face) else {
        return Vec::new();
    };
    std::iter::once(f.outer_loop)
        .chain(f.inner_loops.iter().copied())
        .flat_map(|l| loop_half_edges(b, l))
        .collect()
}

pub(super) fn face_edges(b: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let mut edges: Vec<EdgeHandle> = face_half_edges(b, face)
        .into_iter()
        .map(|h| b.half_edges[h].edge)
        .collect();
    edges.sort();
    edges.dedup();
    edges
}

pub(super) fn mate_face(b: &BRep, h: HalfEdgeHandle) -> Option<FaceHandle> {
    let m = b.half_edges.get(h)?.mate()?;
    Some(b.loops.get(b.half_edges.get(m)?.loop_)?.face)
}

pub(super) fn face_neighbours(b: &BRep, face: FaceHandle) -> Vec<FaceHandle> {
    let mut out: Vec<FaceHandle> = face_half_edges(b, face)
        .into_iter()
        .filter_map(|h| mate_face(b, h))
        .collect();
    out.sort();
    out.dedup();
    out
}

pub(super) fn half_edge_from(
    b: &BRep,
    face: FaceHandle,
    v: VertexHandle,
) -> Option<HalfEdgeHandle> {
    face_half_edges(b, face)
        .into_iter()
        .find(|&h| b.half_edges[h].origin == v)
}

pub(super) fn edge_between(b: &BRep, a: VertexHandle, c: VertexHandle) -> Option<EdgeHandle> {
    b.half_edges.iter().find_map(|(h, he)| {
        let next = b.half_edges.get(he.next)?;
        (he.origin == a && next.origin == c).then(|| b.half_edges[h].edge)
    })
}

/// The face on the side of `edge` that `keep` accepts.
pub(super) fn side_face(
    b: &BRep,
    edge: EdgeHandle,
    keep: impl Fn(FaceHandle) -> bool,
) -> Option<FaceHandle> {
    let h = b.edges.get(edge)?.half_edge;
    let here = b.loops.get(b.half_edges.get(h)?.loop_)?.face;
    let there = mate_face(b, h);
    [Some(here), there]
        .into_iter()
        .flatten()
        .find(|&face| keep(face))
}

pub(super) fn reversed(curve: &NurbsCurve) -> NurbsCurve {
    let (a, z) = curve.parameter_range();
    let mut out = curve.clone();
    out.control_points.reverse();
    out.weights.reverse();
    out.knots = curve.knots.iter().rev().map(|k| a + z - k).collect();
    out
}

pub(super) fn valence(b: &BRep, v: VertexHandle) -> usize {
    let mut edges: Vec<EdgeHandle> = b
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == v)
        .map(|(_, he)| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    edges.len()
}

/// Whether a NURBS curve is a straight segment: every control point on
/// the chord between its ends.
pub(super) fn is_straight(curve: &NurbsCurve, tolerance: &Tolerance) -> bool {
    let (Some(a), Some(z)) = (curve.control_points.first(), curve.control_points.last()) else {
        return false;
    };
    let Some(dir) = (*z - *a).normalized() else {
        return false;
    };
    curve.control_points.iter().all(|p| {
        let d = *p - *a;
        (d - dir * d.dot(&dir)).length() <= tolerance.linear()
    })
}

/// The vertex at `point` on an edge between `faces` and `others`: an end
/// vertex within the reach, else a new one splitting the edge there.
pub(super) fn split_at(
    b: &mut BRep,
    faces: &[FaceHandle],
    others: &[FaceHandle],
    point: Point3,
    reach: f64,
    tolerance: &Tolerance,
) -> Result<VertexHandle, String> {
    let mut best: Option<(f64, EdgeHandle)> = None;
    for h in faces.iter().flat_map(|&face| face_half_edges(b, face)) {
        if !mate_face(b, h).is_some_and(|m| others.contains(&m)) {
            continue;
        }
        let e = b.half_edges[h].edge;
        let Some(c) = b.edges.get(e).and_then(|e| e.curve.as_ref()) else {
            continue;
        };
        let tt = popkern_core::geometry::curve::closest_point_in_window(
            &c.geometry,
            &point,
            c.t_start,
            c.t_end,
        );
        let d = c.geometry.point_at(tt).distance_to(&point);
        if best.is_none_or(|(bd, _)| d < bd) {
            best = Some((d, e));
        }
    }
    let (d, e) = best.ok_or("no edge between the faces")?;
    if d > reach {
        return Err(format!(
            "the crossing lies {d:.3e} off the edge between the faces"
        ));
    }
    let h = b.edges[e].half_edge;
    for v in [
        b.half_edges[h].origin,
        b.half_edges[b.half_edges[h].next].origin,
    ] {
        if b.vertices[v].position.distance_to(&point) <= reach {
            return Ok(v);
        }
    }
    let split = b
        .split_edge_at_point(e, point, tolerance)
        .map_err(|err| format!("split: {err}"))?;
    let landed = b.vertices[split.vertex].position.distance_to(&point);
    if landed > reach {
        return Err(format!("the split landed {landed:.3e} from its point"));
    }
    Ok(split.vertex)
}

/// A new edge from `a` to `c` with `curve`, splitting the face carrying
/// `surface` that holds both.
pub(super) fn make_edge_in(
    b: &mut BRep,
    surface: &SurfaceGeometry,
    a: VertexHandle,
    c: VertexHandle,
    curve: &NurbsCurve,
) -> Result<EdgeHandle, String> {
    let faces = faces_with_surface(b, surface);
    make_edge_among(b, &faces, a, c, curve)
}

pub(super) fn make_edge_among(
    b: &mut BRep,
    faces: &[FaceHandle],
    a: VertexHandle,
    c: VertexHandle,
    curve: &NurbsCurve,
) -> Result<EdgeHandle, String> {
    for &face in faces {
        let (Some(ha), Some(hc)) = (half_edge_from(b, face, a), half_edge_from(b, face, c)) else {
            continue;
        };
        if b.half_edges[ha].loop_ != b.half_edges[hc].loop_ {
            continue;
        }
        let made = b.mef(ha, hc).map_err(|e| format!("mef: {e}"))?;
        // The stored curve runs along the edge's primary half-edge.
        let primary = b.edges[made.edge].half_edge;
        let from = b.vertices[b.half_edges[primary].origin].position;
        let (c0, c1) = curve.parameter_range();
        let forward =
            curve.point_at(c0).distance_to(&from) <= curve.point_at(c1).distance_to(&from);
        let geometry = if forward {
            curve.clone()
        } else {
            reversed(curve)
        };
        let (g0, g1) = geometry.parameter_range();
        if let Some(edge) = b.edges.get_mut(made.edge) {
            edge.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(geometry),
                t_start: g0,
                t_end: g1,
            });
        }
        return Ok(made.edge);
    }
    Err("no face holds both ends of the new edge in one loop".into())
}

/// Contract the zero-length bridge whose half-edge is `h`, keeping the
/// rail's vertex.
pub(super) fn contract_onto_rail(
    b: &mut BRep,
    h: HalfEdgeHandle,
    rail_vertex: VertexHandle,
) -> Result<(), String> {
    let edge = b.half_edges.get(h).ok_or("bridge gone")?.edge;
    let primary = b.edges[edge].half_edge;
    let candidates = [
        primary,
        b.half_edges[primary].mate().ok_or("bridge has no mate")?,
    ];
    let from_rail = candidates
        .into_iter()
        .find(|&c| b.half_edges[c].origin == rail_vertex)
        .ok_or("bridge does not start at the rail")?;
    b.kve(from_rail).map_err(|e| format!("kve: {e}"))?;
    Ok(())
}
