//! The corner regions `replay_delete_faces` closes a deleted face's gap
//! with: along a straight spine, about a circular one, and at a corner
//! ball. Each returns the region, the corner its neighbours meet in and
//! the middle of the face's section: the face's outward normal there
//! points toward the corner when the region is added.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_sketch::profile::PathSegment;

use super::{edge_points, extrude_section, failed, frame_of, revolve_section};
use crate::error::ReplayError;

/// The region, the corner and the face's section middle.
type Region = (BRep, Point3, Point3);

/// The face's boundary edges with the face across each (`None` for a
/// seam, both sides on the face).
pub(super) fn boundary(
    brep: &BRep,
    face: FaceHandle,
) -> Vec<(popkern_core::topology::EdgeHandle, Option<FaceHandle>)> {
    let Some(f) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            let across = h
                .mate()
                .and_then(|m| brep.half_edges.get(m))
                .and_then(|m| brep.loops.get(m.loop_))
                .map(|l| l.face)
                .filter(|other| *other != face);
            out.push((h.edge, across));
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    out
}

/// The direction of the face's longest straight edge.
pub(super) fn longest_straight_edge(brep: &BRep, face: FaceHandle) -> Option<Vec3> {
    boundary(brep, face)
        .into_iter()
        .filter_map(|(e, _)| edge_points(brep, e))
        .filter(|[a, m, b]| {
            let mid = *a + (*b - *a) * 0.5;
            mid.distance_to(m) <= 1e-9 * (1.0 + a.distance_to(b))
        })
        .max_by(|x, y| x[0].distance_to(&x[2]).total_cmp(&y[0].distance_to(&y[2])))
        .and_then(|[a, _, b]| (b - a).normalized())
}

/// The point where the lines through `p` along `u` and through `q`
/// along `v` meet, in the plane.
fn meet(p: [f64; 2], u: [f64; 2], q: [f64; 2], v: [f64; 2]) -> Option<[f64; 2]> {
    let det = u[0] * v[1] - u[1] * v[0];
    if det.abs() < 1e-12 {
        return None;
    }
    let t = ((q[0] - p[0]) * v[1] - (q[1] - p[1]) * v[0]) / det;
    Some([p[0] + t * u[0], p[1] + t * u[1]])
}

/// The section triangle `p_a` → corner → `p_b`, and the middle of the
/// face's section (on its arc, or its chord), in section coordinates.
fn section(
    p_a: [f64; 2],
    corner: [f64; 2],
    p_b: [f64; 2],
    arc: Option<([f64; 2], f64)>,
) -> Option<(Vec<PathSegment<f64>>, [f64; 2])> {
    let mid = [0.5 * (p_a[0] + p_b[0]), 0.5 * (p_a[1] + p_b[1])];
    let face_mid = match arc {
        Some((c, r)) => {
            let (dx, dy) = (mid[0] - c[0], mid[1] - c[1]);
            let l = (dx * dx + dy * dy).sqrt();
            if l < 1e-12 {
                return None; // a half circle: its side is not decided
            }
            [c[0] + r * dx / l, c[1] + r * dy / l]
        }
        None => mid,
    };
    // The whole triangle to the corner: the part of it across the face
    // from the corner is body where the region is added (a convex blend
    // bulges toward its corner) and air where it is taken away, so the
    // boolean reads the same solid, and the face's own curved surface,
    // which the body holds exactly, never meets itself in it.
    Some((
        vec![
            PathSegment::Line { to: corner },
            PathSegment::Line { to: p_b },
            PathSegment::Line { to: p_a },
        ],
        face_mid,
    ))
}

/// A round (cylinder about `arc`'s axis) or a bevel along the straight
/// spine `d`: its two edges along the spine, each meeting a plane that
/// holds the spine direction; its other edges in planes square to it.
pub(super) fn straight(
    brep: &BRep,
    face: FaceHandle,
    d: Vec3,
    arc: Option<(Point3, f64)>,
    tolerance: Tolerance,
) -> Result<Region, ReplayError> {
    let lin = tolerance.linear();
    let anchor = arc.map(|(o, _)| o).unwrap_or_else(|| {
        brep.face_outer_loop_positions(face)
            .first()
            .copied()
            .unwrap_or(Point3::new(0.0, 0.0, 0.0))
    });
    let (d, e1, frame) = frame_of(anchor, d).ok_or_else(|| failed("degenerate spine"))?;
    let e2 = d.cross(&e1);
    let to2 = |p: Point3| [(p - anchor).dot(&e1), (p - anchor).dot(&e2)];
    let along = |p: Point3| (p - anchor).dot(&d);
    let mut long: Vec<([f64; 2], FaceHandle)> = Vec::new();
    let (mut z0, mut z1) = (f64::INFINITY, f64::NEG_INFINITY);
    for (edge, across) in boundary(brep, face) {
        let pts = edge_points(brep, edge).ok_or_else(|| failed("an edge has no curve"))?;
        for p in &pts {
            z0 = z0.min(along(*p));
            z1 = z1.max(along(*p));
        }
        let chord = pts[2] - pts[0];
        let lies_along = chord.cross(&d).length() <= lin * (1.0 + chord.length())
            && chord.length() > lin
            && (pts[1] - pts[0]).cross(&d).length() <= lin * (1.0 + chord.length());
        let squared = pts.iter().all(|p| (along(*p) - along(pts[0])).abs() <= lin);
        if lies_along {
            let across = across.ok_or_else(|| failed("a seam along the spine"))?;
            long.push((to2(pts[0]), across));
        } else if !squared {
            return Err(failed("an end of the blend is not square to its spine"));
        }
    }
    let [(p_a, fa), (p_b, fb)] = long[..] else {
        return Err(failed(format!(
            "a straight blend has two edges along its spine, not {}",
            long.len()
        )));
    };
    let line_dir = |f: FaceHandle| -> Result<[f64; 2], ReplayError> {
        match brep
            .faces
            .get(f)
            .and_then(|x| x.surface.as_ref())
            .map(|s| &s.geometry)
        {
            Some(SurfaceGeometry::Plane { normal, .. }) if normal.dot(&d).abs() <= 1e-9 => {
                let t = normal.cross(&d);
                Ok([t.dot(&e1), t.dot(&e2)])
            }
            _ => Err(failed(
                "a neighbour along the spine is not a plane holding it",
            )),
        }
    };
    let corner = meet(p_a, line_dir(fa)?, p_b, line_dir(fb)?)
        .ok_or_else(|| failed("the neighbours are parallel"))?;
    let arc2 = arc.map(|(o, r)| (to2(o), r));
    let (segments, mid) =
        section(p_a, corner, p_b, arc2).ok_or_else(|| failed("a half-round blend"))?;
    let region = extrude_section(p_a, segments, (z0, z1 - z0), &frame, tolerance)?;
    let at = |q: [f64; 2]| anchor + e1 * q[0] + e2 * q[1] + d * (0.5 * (z0 + z1));
    Ok((region, at(corner), at(mid)))
}

/// A round (torus, its tube `tube` = (major, minor) about `centre`) or a
/// bevel (cone) about the axis through `centre` along `axis`: its edges
/// circles about the axis, each meeting a plane square to it or a
/// cylinder or cone about it.
pub(super) fn circular(
    brep: &BRep,
    face: FaceHandle,
    centre: Point3,
    axis: Vec3,
    tube: Option<(f64, f64)>,
    tolerance: Tolerance,
) -> Result<Region, ReplayError> {
    let lin = tolerance.linear();
    let (a, x, frame) = frame_of(centre, axis).ok_or_else(|| failed("degenerate axis"))?;
    let rh = |p: Point3| {
        let h = (p - centre).dot(&a);
        [((p - centre) - a * h).length(), h]
    };
    let mut rims: Vec<([f64; 2], FaceHandle)> = Vec::new();
    for (edge, across) in boundary(brep, face) {
        let Some(across) = across else {
            continue; // the seam
        };
        let pts = edge_points(brep, edge).ok_or_else(|| failed("an edge has no curve"))?;
        let q = pts.map(rh);
        if q.iter()
            .any(|s| (s[0] - q[0][0]).abs() > lin || (s[1] - q[0][1]).abs() > lin)
        {
            return Err(failed(
                "an edge of the blend is not a circle about its axis",
            ));
        }
        rims.push((q[0], across));
    }
    let [(p_a, fa), (p_b, fb)] = rims[..] else {
        return Err(failed(format!(
            "a blend about an axis has two rims, not {}",
            rims.len()
        )));
    };
    let line_dir = |f: FaceHandle, p: [f64; 2]| -> Result<[f64; 2], ReplayError> {
        let parallel = |v: &Vec3| v.normalized().is_some_and(|v| v.cross(&a).length() <= 1e-9);
        let on_axis = |o: &Point3| ((*o - centre) - a * (*o - centre).dot(&a)).length() <= lin;
        match brep
            .faces
            .get(f)
            .and_then(|x| x.surface.as_ref())
            .map(|s| &s.geometry)
        {
            Some(SurfaceGeometry::Plane { normal, .. }) if parallel(normal) => Ok([1.0, 0.0]),
            Some(SurfaceGeometry::Cylinder { origin, axis, .. })
                if parallel(axis) && on_axis(origin) =>
            {
                Ok([0.0, 1.0])
            }
            Some(SurfaceGeometry::Cone { apex, axis, .. }) if parallel(axis) && on_axis(apex) => {
                let h = (*apex - centre).dot(&a);
                Ok([p[0], p[1] - h])
            }
            _ => Err(failed(
                "a neighbour is not a plane square to the axis, nor a cylinder or cone about it",
            )),
        }
    };
    let corner = meet(p_a, line_dir(fa, p_a)?, p_b, line_dir(fb, p_b)?)
        .ok_or_else(|| failed("the neighbours are parallel"))?;
    let arc2 = tube.map(|(major, minor)| ([major, 0.0], minor));
    let (segments, mid) =
        section(p_a, corner, p_b, arc2).ok_or_else(|| failed("a half-round blend"))?;
    if [p_a, corner, p_b].iter().any(|p| p[0] < -lin) {
        return Err(failed("the corner falls across the axis"));
    }
    let region = revolve_section(p_a, segments, &frame, tolerance)?;
    let at = |q: [f64; 2]| centre + x * q[0] + a * q[1];
    Ok((region, at(corner), at(mid)))
}

/// A corner ball: the patch of the sphere about `centre` whose three
/// corners lie along orthogonal directions from it; the region is the
/// cube those directions span out to `radius`, less the ball.
pub(super) fn corner_ball(
    brep: &BRep,
    face: FaceHandle,
    centre: Point3,
    radius: f64,
    tolerance: Tolerance,
) -> Result<Region, ReplayError> {
    let corners = brep.face_outer_loop_positions(face);
    if corners.len() != 3 {
        return Err(failed("a corner ball has three corners"));
    }
    let e: Vec<Vec3> = corners
        .iter()
        .filter_map(|p| (*p - centre).normalized())
        .collect();
    if e.len() != 3 || (0..3).any(|i| e[i].dot(&e[(i + 1) % 3]).abs() > 1e-9) {
        return Err(failed(
            "a corner ball's corners lie along orthogonal directions",
        ));
    }
    let (z, x) = if e[0].cross(&e[1]).dot(&e[2]) > 0.0 {
        (e[2], e[0])
    } else {
        (e[2], e[1])
    };
    let frame =
        Transform3::from_axis_placement(centre, z, x).map_err(|er| failed(er.to_string()))?;
    let square = vec![
        PathSegment::Line { to: [radius, 0.0] },
        PathSegment::Line {
            to: [radius, radius],
        },
        PathSegment::Line { to: [0.0, radius] },
        PathSegment::Line { to: [0.0, 0.0] },
    ];
    let cube = extrude_section([0.0, 0.0], square, (0.0, radius), &frame, tolerance)?;
    let ball = BRep::sphere(radius, tolerance)
        .and_then(|b| b.transform(&Transform3::translation(centre.x, centre.y, centre.z)))
        .map_err(|er| failed(er.to_string()))?;
    let region = popkern_bool::subtraction(cube, ball, tolerance, None)
        .map_err(|er| failed(format!("the corner ball's region: {er}")))?;
    let out = e[0] + e[1] + e[2];
    let mid = centre + out.normalized().unwrap_or(out) * radius;
    Ok((region, centre + out * radius, mid))
}

/// That every neighbour of `face` stands square to it along their shared
/// edges, so it runs on along its own surface as the face moves: at
/// each edge's middle, the neighbour's normal lies in the face's tangent
/// plane.
pub(super) fn require_square_neighbours(
    brep: &BRep,
    face: FaceHandle,
    tolerance: Tolerance,
) -> Result<(), ReplayError> {
    use popkern_core::geometry::surface::SurfaceEval;
    let normal = |f: FaceHandle, p: &Point3| -> Option<Vec3> {
        let s = &brep.faces.get(f)?.surface.as_ref()?.geometry;
        let (u, v) = s.closest_point(p)?;
        s.normal_at(u, v).normalized()
    };
    for (edge, across) in boundary(brep, face) {
        let Some(across) = across else {
            continue;
        };
        let [_, mid, _] = edge_points(brep, edge).ok_or_else(|| failed("an edge has no curve"))?;
        let (Some(n_face), Some(n_across)) = (normal(face, &mid), normal(across, &mid)) else {
            return Err(failed("a normal at an edge cannot be read"));
        };
        if n_face.dot(&n_across).abs() > tolerance.angular() {
            return Err(failed(
                "a neighbour leans onto the face; only faces whose neighbours stand square to them move",
            ));
        }
    }
    Ok(())
}

/// The region `face` sweeps moving `distance` along its outward normal:
/// a flat face's prism (the face swept along its normal), a cylinder
/// face's ring between its radius and the moved one over its axial
/// span (a hole's wall, outward normal toward the axis, widens as it
/// moves in).
pub(super) fn face_sweep(
    brep: &BRep,
    face: FaceHandle,
    distance: f64,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    use popkern_core::geometry::surface::SurfaceEval;
    let sof = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .ok_or_else(|| failed("the face has no surface"))?;
    let sign = if sof.orientation { 1.0 } else { -1.0 };
    let points = brep.face_boundary_sample_points(face, 8);
    let first = *points
        .first()
        .ok_or_else(|| failed("the face has no boundary"))?;
    match &sof.geometry {
        SurfaceGeometry::Plane { normal, .. } => {
            let dir = *normal * (sign * distance.signum());
            let (_, _, frame) = frame_of(first, dir).ok_or_else(|| failed("degenerate normal"))?;
            let to_local = frame.inverse().map_err(|e| failed(e.to_string()))?;
            let local = brep
                .transform(&to_local)
                .map_err(|e| failed(e.to_string()))?;
            let n = points.len() as f64;
            let c = points
                .iter()
                .map(|p| to_local.apply_point(p))
                .fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| {
                    acc + (p - Point3::new(0.0, 0.0, 0.0))
                })
                * (1.0 / n);
            let start = Point3::new(c.x, c.y, 0.0);
            let path = popkern_core::geometry::NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![start, start + Vec3::new(0.0, 0.0, distance.abs())],
                weights: vec![1.0, 1.0],
            };
            let swept = popkern_sweep::sweep(
                &popkern_sweep::SweepProfile::Face {
                    brep: Box::new(local),
                    face,
                },
                &path,
                &popkern_sweep::SweepOptions::default(),
                tolerance,
                None,
            )
            .map_err(|e| failed(format!("the face's sweep: {e}")))?;
            swept
                .brep
                .transform(&frame)
                .map_err(|e| failed(e.to_string()))
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let a = axis.normalized().ok_or_else(|| failed("degenerate axis"))?;
            let (h0, h1) = points
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                    let h = (*p - *origin).dot(&a);
                    (lo.min(h), hi.max(h))
                });
            // The outward normal at the face: radial, away from the axis
            // for a boss and toward it for a hole.
            let (u, v) = sof
                .geometry
                .closest_point(&first)
                .ok_or_else(|| failed("no chart point"))?;
            let n = sof.geometry.normal_at(u, v) * sign;
            let radial = (first - *origin) - a * (first - *origin).dot(&a);
            let outward_from_axis = n.dot(&radial) > 0.0;
            let moved = if outward_from_axis {
                radius + distance
            } else {
                radius - distance
            };
            if moved <= tolerance.linear() {
                return Err(failed("the moved cylinder closes on its axis"));
            }
            let (inner, outer) = (radius.min(moved), radius.max(moved));
            let base = *origin + a * h0;
            let (_, _, frame) = frame_of(base, a).ok_or_else(|| failed("degenerate axis"))?;
            let solid = |r: f64| {
                BRep::cylinder(r, h1 - h0, tolerance)
                    .and_then(|c| c.transform(&frame))
                    .map_err(|e| failed(e.to_string()))
            };
            popkern_bool::subtraction(solid(outer)?, solid(inner)?, tolerance, None)
                .map_err(|e| failed(format!("the wall's ring: {e}")))
        }
        _ => Err(failed("only flat faces and cylinders move")),
    }
}
