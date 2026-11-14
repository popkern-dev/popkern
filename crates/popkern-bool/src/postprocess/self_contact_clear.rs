//! Closed-form clears of a face pair for the self-contact search
//! (PK-1155): two faces of one valid shell sharing no vertex do not
//! cross, so their nearest points lie on a boundary of one, or on both
//! interiors with the segment between normal to both. Planes and
//! cylinders whose normals cannot be parallel within the band are
//! cleared by their boundaries alone.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Whether `c` is a segment: a line, or a NURBS of degree one between
/// two points (its weights move the parameter, not the points).
pub(super) fn straight(c: &CurveOnEdge) -> bool {
    match &c.geometry {
        CurveGeometry::Line { .. } => true,
        CurveGeometry::Nurbs(n) => n.degree == 1 && n.control_points.len() == 2,
        _ => false,
    }
}

/// A face of a pair as the closed-form clears read it: its handle,
/// surface, a box holding it and its boundary's edges.
pub(super) type Side<'a> = (FaceHandle, &'a SurfaceGeometry, &'a Aabb, &'a [EdgeHandle]);

/// A plane's point and unit normal.
pub(super) fn plane_of(s: &SurfaceGeometry) -> Option<(Point3, Vec3)> {
    match s {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((*origin, normal.normalized()?)),
        _ => None,
    }
}

/// Whether no part of any of `edges` within `tol` of the plane lies in
/// the box `near` (padded by `tol`): each edge's box (its curve's
/// conservative extent) lies wholly beyond the band on one side, or a
/// segment's stretch inside the band lies outside the box, or outside
/// the planar face `face` itself (its middle outside the face's region
/// and farther than the band and half the stretch from its boundary:
/// then no point of it crosses the boundary, and its distance to the
/// face is its distance to the boundary).
pub(super) fn edges_clear_of_plane(
    brep: &BRep,
    edges: &[EdgeHandle],
    (origin, n): (Point3, Vec3),
    (face, near): (FaceHandle, &Aabb),
    tol: f64,
) -> bool {
    edges.iter().all(|&e| {
        brep.edges
            .get(e)
            .and_then(|x| x.curve.as_ref())
            .is_some_and(|c| {
                curve_clear_of_plane(c, (origin, n), near, tol)
                    || stretch_off_face(brep, face, c, (origin, n), tol)
                    || curve_off_face(brep, face, c, (origin, n), tol)
            })
    })
}

/// Whether the stretch of the straight edge `c` within `tol` of the
/// plane of planar `face` lies off the face by more than `tol` (see
/// [`edges_clear_of_plane`]).
fn stretch_off_face(
    brep: &BRep,
    face: FaceHandle,
    c: &CurveOnEdge,
    (origin, n): (Point3, Vec3),
    tol: f64,
) -> bool {
    if !straight(c) {
        return false;
    }
    let (p, q) = (c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end));
    let Some((a, b)) = band_stretch(p, q, (origin, n), tol) else {
        return true;
    };
    ball_off_face(brep, face, a + (b - a) * 0.5, tol + 0.5 * a.distance_to(&b))
}

/// Halvings of a curved edge whose pieces lie in the plane's band.
const OFF_FACE_DEPTH: usize = 8;

/// Whether the curved edge `c` lies farther than `tol` from the planar
/// `face` (on the plane `(origin, n)`): each piece of its NURBS form, of
/// positive weights, halved while undecided, lies beyond the band on
/// one side (its control points do), or within the ball about its
/// control box's middle of its half-diagonal, off the face by `tol`
/// ([`ball_off_face`]). A torus's meridian lying in the plane beside
/// the planar face is in the band everywhere, off the face throughout.
fn curve_off_face(
    brep: &BRep,
    face: FaceHandle,
    c: &CurveOnEdge,
    (origin, n): (Point3, Vec3),
    tol: f64,
) -> bool {
    use popkern_core::geometry::nurbs_curve::NurbsCurve;
    fn piece_off(
        brep: &BRep,
        face: FaceHandle,
        curve: &NurbsCurve,
        (origin, n): (Point3, Vec3),
        tol: f64,
        depth: usize,
    ) -> bool {
        let Some(hull) = Aabb::from_points(&curve.control_points) else {
            return false;
        };
        let side = |p: &Point3| (*p - origin).dot(&n);
        if curve.control_points.iter().all(|p| side(p) > tol)
            || curve.control_points.iter().all(|p| side(p) < -tol)
        {
            return true;
        }
        let middle = hull.min + (hull.max - hull.min) * 0.5;
        if ball_off_face(
            brep,
            face,
            middle,
            tol + 0.5 * hull.min.distance_to(&hull.max),
        ) {
            return true;
        }
        if depth == 0 {
            return false;
        }
        let p = curve.degree;
        let (Some(&t0), Some(&t1)) = (
            curve.knots.get(p),
            curve
                .knots
                .len()
                .checked_sub(p + 1)
                .and_then(|i| curve.knots.get(i)),
        ) else {
            return false;
        };
        curve.split_at(0.5 * (t0 + t1)).is_ok_and(|(a, b)| {
            piece_off(brep, face, &a, (origin, n), tol, depth - 1)
                && piece_off(brep, face, &b, (origin, n), tol, depth - 1)
        })
    }
    if straight(c) {
        return false;
    }
    let sub = c.geometry.to_nurbs_in_range((c.t_start, c.t_end));
    sub.weights.iter().all(|&w| w > 0.0)
        && piece_off(brep, face, &sub, (origin, n), tol, OFF_FACE_DEPTH)
}

/// Whether every point within `reach - pad` of `centre` lies farther
/// than `pad` from planar `face`, `reach` the ball's radius and the pad:
/// the centre's foot on the plane lies outside the face's region and
/// farther than `reach` from its boundary. The boundary lies in the
/// plane, so each point of the ball projects within its radius of the
/// foot, outside the region, and lies at least as far from the boundary
/// as its projection, by more than the pad.
pub(super) fn ball_off_face(brep: &BRep, face: FaceHandle, centre: Point3, reach: f64) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let Some((u, v)) =
        popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, &centre)
    else {
        return false;
    };
    let foot = surface.geometry.point_at(u, v);
    crate::classify_uv_loops::point_in_face_region_exact(brep, face, u, v) == Some(false)
        && !brep.point_near_face_boundary(face, &foot, reach)
}

/// The stretch of the segment `p`..`q` within `tol` of the plane, `None`
/// when it has none.
fn band_stretch(
    p: Point3,
    q: Point3,
    (origin, n): (Point3, Vec3),
    tol: f64,
) -> Option<(Point3, Point3)> {
    let (dp, dq) = ((p - origin).dot(&n), (q - origin).dot(&n));
    if (dq - dp).abs() < 1e-300 {
        return (dp.abs() <= tol).then_some((p, q));
    }
    let at = |d: f64| ((d - dp) / (dq - dp)).clamp(0.0, 1.0);
    let (s0, s1) = (at(-tol).min(at(tol)), at(-tol).max(at(tol)));
    if s1 <= s0 && !(dp.abs() <= tol || dq.abs() <= tol) {
        return None;
    }
    Some((p + (q - p) * s0, p + (q - p) * s1))
}

/// Whether a planar face (`plane`, its box and edges) and a cylindrical
/// face (`other`) provably lie farther than `tol` apart: the cylinder's
/// axis makes an angle with the plane (then no normals of the two are
/// parallel, and the faces' nearest points lie on a boundary), no edge
/// of the cylinder's face comes within `tol` of the plane inside the
/// planar face's box, and no edge of the planar face within `tol` of
/// the cylinder face's box (each curve halved over its control hull).
pub(super) fn plane_crosses_cylinder_clear(
    brep: &BRep,
    (plane_face, plane, plane_box, plane_edges): Side,
    (_, other, other_box, other_edges): Side,
    tol: f64,
) -> bool {
    let (
        Some(pp),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        },
    ) = (plane_of(plane), other)
    else {
        return false;
    };
    let Some(a) = axis.normalized() else {
        return false;
    };
    // The axis in the plane's directions: the normals are parallel on the
    // rulings nearest and farthest from the plane, |d - r| and d + r from
    // it (d the axis's distance, read beside the face's box; a tilt under
    // the threshold moves it by at most the tilt across the box).
    let tilt = a.dot(&pp.1).abs();
    let normals_apart = tilt > AXIS_IN_PLANE || {
        let middle = other_box.min + (other_box.max - other_box.min) * 0.5;
        let on_axis = *origin + a * (middle - *origin).dot(&a);
        let d = (on_axis - pp.0).dot(&pp.1).abs();
        let slack = tilt * other_box.min.distance_to(&other_box.max);
        (d - radius.abs()).abs() > tol + slack && d + radius.abs() > tol + slack
    };
    normals_apart
        && edges_clear_of_plane(brep, other_edges, pp, (plane_face, plane_box), tol)
        && plane_edges.iter().all(|&e| {
            brep.edges
                .get(e)
                .and_then(|x| x.curve.as_ref())
                .is_some_and(|c| {
                    let sub = c.geometry.to_nurbs_in_range((c.t_start, c.t_end));
                    crate::parallel::ff_parallel::ff_parallel_reach::curve_clear_of_surface_face(
                        &sub, other, other_box, tol,
                    )
                })
        })
}

/// Whether two cylindrical faces provably lie farther than `tol`
/// apart. Their normals are parallel only along `m`, the common
/// perpendicular of their axes (for parallel axes, the direction from
/// one to the other), on the rulings `r` from each axis along `m`:
/// interior points with parallel normals lie `|D ± r1 ± r2|` apart, `D`
/// the axes' separation along `m` (near-parallel axes read as
/// parallel, with the tilt across the faces' boxes as slack). When
/// none is within `tol` the faces' nearest points lie on a boundary,
/// and no edge of either may come within `tol` of the other's face.
pub(super) fn cylinders_clear(
    brep: &BRep,
    (_, s1, box1, edges1): Side,
    (_, s2, box2, edges2): Side,
    tol: f64,
) -> bool {
    let (
        SurfaceGeometry::Cylinder {
            origin: o1,
            axis: x1,
            radius: r1,
        },
        SurfaceGeometry::Cylinder {
            origin: o2,
            axis: x2,
            radius: r2,
        },
    ) = (s1, s2)
    else {
        return false;
    };
    let (Some(a1), Some(a2)) = (x1.normalized(), x2.normalized()) else {
        return false;
    };
    let tilt = a1.cross(&a2).length();
    let mid = |b: &Aabb| b.min + (b.max - b.min) * 0.5;
    let (c1, c2) = (mid(box1), mid(box2));
    let (p1, p2) = (
        *o1 + a1 * (c1 - *o1).dot(&a1),
        *o2 + a2 * (c2 - *o2).dot(&a2),
    );
    let (m, slack) = if tilt > AXIS_IN_PLANE {
        (a1.cross(&a2) * (1.0 / tilt), 0.0)
    } else {
        let off = (p2 - p1) - a1 * (p2 - p1).dot(&a1);
        let diag =
            box1.min.distance_to(&box1.max) + box2.min.distance_to(&box2.max) + c1.distance_to(&c2);
        // Coaxial: any direction across the axis (the separation is 0).
        let across = || {
            let e = if a1.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            (e - a1 * e.dot(&a1)).normalized()
        };
        let Some(m) = off.normalized().or_else(across) else {
            return false;
        };
        (m, tilt * diag)
    };
    let d = (p2 - p1).dot(&m);
    let (r1, r2) = (r1.abs(), r2.abs());
    let normals_apart = [d + r1 + r2, d + r1 - r2, d - r1 + r2, d - r1 - r2]
        .iter()
        .all(|x| x.abs() > tol + slack);
    let edges_clear = |edges: &[EdgeHandle], surface: &SurfaceGeometry, bx: &Aabb| {
        edges.iter().all(|&e| {
            brep.edges
                .get(e)
                .and_then(|x| x.curve.as_ref())
                .is_some_and(|c| {
                    let sub = c.geometry.to_nurbs_in_range((c.t_start, c.t_end));
                    crate::parallel::ff_parallel::ff_parallel_reach::curve_clear_of_surface_face(
                        &sub, surface, bx, tol,
                    )
                })
        })
    };
    normals_apart && edges_clear(edges1, s2, box2) && edges_clear(edges2, s1, box1)
}

/// A cylinder's axis this close to the plane's directions (the cosine
/// of its angle to the normal) may be tangent to it: not cleared.
const AXIS_IN_PLANE: f64 = 1e-6;

/// Whether no point of `c` within `tol` of the plane lies in the box
/// `near` padded by `tol` (see [`edges_clear_of_plane`]).
pub(super) fn curve_clear_of_plane(
    c: &CurveOnEdge,
    (origin, n): (Point3, Vec3),
    near: &Aabb,
    tol: f64,
) -> bool {
    let side = |p: &Point3| (*p - origin).dot(&n);
    {
        let (p, q) = (c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end));
        if straight(c) {
            // The stretch of the segment within the band.
            return band_stretch(p, q, (origin, n), tol).is_none_or(|(a, b)| {
                Aabb::from_points(&[a, b]).is_none_or(|stretch| !stretch.overlaps(near, tol))
            });
        }
        // The edge lies in its control hull (a NURBS of positive
        // weights over the edge's span) or, for a conic, between its
        // centre's height plus and minus its axes' amplitude.
        let heights: Option<(f64, f64)> = match &c.geometry {
            CurveGeometry::Nurbs(_) => {
                let sub = c.geometry.to_nurbs_in_range((c.t_start, c.t_end));
                (sub.weights.iter().all(|&w| w > 0.0) && !sub.control_points.is_empty()).then(
                    || {
                        sub.control_points
                            .iter()
                            .map(&side)
                            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), d| {
                                (lo.min(d), hi.max(d))
                            })
                    },
                )
            }
            CurveGeometry::Circle { center, .. }
            | CurveGeometry::Arc { center, .. }
            | CurveGeometry::Ellipse { center, .. } => {
                use popkern_core::geometry::PortableMath;
                let u = c.geometry.point_at(0.0) - *center;
                let v = c.geometry.point_at(std::f64::consts::FRAC_PI_2) - *center;
                let amp = u.dot(&n).phypot(v.dot(&n));
                let h = side(center);
                Some((h - amp, h + amp))
            }
            _ => None,
        };
        if let Some((lo, hi)) = heights
            && (lo > tol || hi < -tol)
        {
            return true;
        }
        let mut points = vec![p, q];
        crate::face_box::conservative_curve_extent(c, &mut points);
        let Some(b) = Aabb::from_points(&points) else {
            return false;
        };
        let corners = [
            Point3::new(b.min.x, b.min.y, b.min.z),
            Point3::new(b.max.x, b.min.y, b.min.z),
            Point3::new(b.min.x, b.max.y, b.min.z),
            Point3::new(b.max.x, b.max.y, b.min.z),
            Point3::new(b.min.x, b.min.y, b.max.z),
            Point3::new(b.max.x, b.min.y, b.max.z),
            Point3::new(b.min.x, b.max.y, b.max.z),
            Point3::new(b.max.x, b.max.y, b.max.z),
        ];
        let d = corners.map(|c| side(&c));
        d.iter().all(|&x| x > tol) || d.iter().all(|&x| x < -tol) || !b.overlaps(near, tol)
    }
}
