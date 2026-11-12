use super::clear::{curve_clear_of_plane, straight};
use super::*;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::transform::Transform3;

fn plane() -> (Point3, Vec3) {
    (
        Point3::new(0.1, -0.2, 0.05),
        Vec3::new(0.3, -0.5, 0.8).normalized().expect("normal"),
    )
}

fn near() -> Aabb {
    Aabb {
        min: Point3::new(-1.0, -0.8, -0.6),
        max: Point3::new(0.9, 1.1, 0.7),
    }
}

fn segment(p: Point3, q: Point3) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: p,
            direction: q - p,
        },
        t_start: 0.0,
        t_end: 1.0,
    }
}

fn nurbs(points: Vec<Point3>, weights: Vec<f64>) -> CurveOnEdge {
    let degree = points.len() - 1;
    let mut knots = vec![0.0; degree + 1];
    knots.extend(vec![1.0; degree + 1]);
    CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree,
            knots,
            control_points: points,
            weights,
        }),
        t_start: 0.0,
        t_end: 1.0,
    }
}

fn arc(c: Point3, normal: Vec3, radius: f64, a: f64, b: f64) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center: c,
            normal,
            radius,
            start_angle: a,
            end_angle: b,
        },
        t_start: a,
        t_end: b,
    }
}

/// Every curve the plane test clears keeps each of 2000 points either
/// farther than the band from the plane or outside the padded box:
/// segments, straight NURBS of unequal weights, arcs and a rational
/// quadratic over a sweep of placements, each form both cleared and
/// kept somewhere.
#[test]
fn a_cleared_curve_never_comes_within_the_band_inside_the_box() {
    let tol = 0.05;
    let (origin, n) = plane();
    let b = near();
    let inside = |p: &Point3| {
        p.x >= b.min.x - tol
            && p.x <= b.max.x + tol
            && p.y >= b.min.y - tol
            && p.y <= b.max.y + tol
            && p.z >= b.min.z - tol
            && p.z <= b.max.z + tol
    };
    let forms: [&dyn Fn(Point3) -> CurveOnEdge; 4] = [
        &|c| segment(c, c + Vec3::new(0.7, 0.4, -0.3)),
        &|c| nurbs(vec![c, c + Vec3::new(-0.2, 0.9, 0.4)], vec![1.0, 3.0]),
        &|c| arc(c, Vec3::new(0.0, 0.6, 0.8), 0.6, 0.3, 2.6),
        &|c| {
            nurbs(
                vec![
                    c,
                    c + Vec3::new(0.5, 0.6, 0.1),
                    c + Vec3::new(1.0, 0.0, -0.2),
                ],
                vec![1.0, 0.5, 1.0],
            )
        },
    ];
    for (k, form) in forms.iter().enumerate() {
        let (mut cleared, mut kept) = (0, 0);
        for i in 0..9 {
            for j in 0..9 {
                for l in 0..5 {
                    let c = Point3::new(
                        -2.0 + 0.45 * f64::from(i),
                        -2.0 + 0.45 * f64::from(j),
                        -1.2 + 0.5 * f64::from(l),
                    );
                    let curve = form(c);
                    if !curve_clear_of_plane(&curve, (origin, n), &b, tol) {
                        kept += 1;
                        continue;
                    }
                    cleared += 1;
                    for s in 0..=2000 {
                        let t =
                            curve.t_start + (curve.t_end - curve.t_start) * f64::from(s) / 2000.0;
                        let p = curve.geometry.point_at(t);
                        let d = (p - origin).dot(&n).abs();
                        assert!(
                            d > tol || !inside(&p),
                            "form {k} at {c:?}: {p:?} is {d} from the plane"
                        );
                    }
                }
            }
        }
        assert!(
            cleared > 20 && kept > 10,
            "form {k}: cleared {cleared}, kept {kept}"
        );
    }
}

/// A degree-one NURBS between two points is read as its segment: its
/// weights move the parameter, not the points.
#[test]
fn a_straight_nurbs_is_a_segment() {
    let p = Point3::new(0.0, 0.0, 0.0);
    let q = Point3::new(1.0, 2.0, 0.5);
    assert!(straight(&nurbs(vec![p, q], vec![1.0, 4.0])));
    assert!(straight(&segment(p, q)));
    assert!(!straight(&nurbs(
        vec![p, Point3::new(0.5, 1.5, 0.0), q],
        vec![1.0; 3]
    )));
}

/// The whole-surface boxes of a turned ball and torus hold every point
/// of their surfaces.
#[test]
fn whole_surface_boxes_hold_their_surfaces() {
    let tol = popkern_core::tolerance::Tolerance::default();
    let turn = Transform3::rotation_around_axis(
        Point3::new(1.0, -2.0, 0.5),
        Vec3::new(0.4, 1.0, -0.7),
        0.9,
    )
    .expect("axis");
    for body in [
        BRep::torus(5.0, 1.5, tol).expect("torus"),
        BRep::sphere(3.0, tol).expect("ball"),
    ] {
        let body = body.transform(&turn).expect("turn");
        for f in body.faces.keys() {
            let surface = &body.faces[f].surface.as_ref().expect("surface").geometry;
            if !matches!(
                surface,
                SurfaceGeometry::Torus { .. } | SurfaceGeometry::Sphere { .. }
            ) {
                continue;
            }
            let b = containing_box(&body, f).expect("box");
            let ((u0, u1), (v0, v1)) = surface.parameter_range();
            for i in 0..=60 {
                for j in 0..=60 {
                    let u = u0 + (u1 - u0) * f64::from(i) / 60.0;
                    let v = v0 + (v1 - v0) * f64::from(j) / 60.0;
                    let p = surface.point_at(u, v);
                    let slack = 1e-9;
                    assert!(
                        p.x >= b.min.x - slack
                            && p.x <= b.max.x + slack
                            && p.y >= b.min.y - slack
                            && p.y <= b.max.y + slack
                            && p.z >= b.min.z - slack
                            && p.z <= b.max.z + slack,
                        "{p:?} outside {b:?}"
                    );
                }
            }
        }
    }
}

/// A body of `parts` side by side (`BRep::append_body`): the clears read
/// faces and edges alone, not which shell holds them.
fn side_by_side(parts: &[BRep]) -> BRep {
    let mut out = parts[0].clone();
    for p in &parts[1..] {
        out.append_body(p);
    }
    out
}

fn prism_z(profile: &[(f64, f64)], z0: f64, height: f64) -> BRep {
    let tol = popkern_core::tolerance::Tolerance::default();
    let points: Vec<Point3> = profile
        .iter()
        .map(|&(x, y)| Point3::new(x, y, z0))
        .collect();
    BRep::prism(&points, Vec3::new(0.0, 0.0, 1.0), height, tol).expect("prism")
}

fn cylinder_at(radius: f64, height: f64, turn: Option<Transform3>, at: Vec3) -> BRep {
    let tol = popkern_core::tolerance::Tolerance::default();
    let mut c = BRep::cylinder(radius, height, tol).expect("cylinder");
    if let Some(t) = turn {
        c = c.transform(&t).expect("turn");
    }
    c.transform(&Transform3::translation(at.x, at.y, at.z))
        .expect("place")
}

/// The face whose surface `pick` accepts.
fn face_where(body: &BRep, pick: impl Fn(&SurfaceGeometry) -> bool) -> FaceHandle {
    body.faces
        .keys()
        .find(|&f| pick(&body.faces[f].surface.as_ref().expect("surface").geometry))
        .expect("face")
}

fn top_plane(z: f64) -> impl Fn(&SurfaceGeometry) -> bool {
    move |s| matches!(s, SurfaceGeometry::Plane { origin, normal, .. } if (origin.z - z).abs() < 1e-9 && normal.z.abs() > 0.99)
}

fn lateral(s: &SurfaceGeometry) -> bool {
    matches!(s, SurfaceGeometry::Cylinder { .. })
}

/// Whether the cylinder pair clear (`cylinders`) or the plane and
/// cylinder clear proves faces `f` and `g` apart.
fn clears(body: &BRep, f: FaceHandle, g: FaceHandle, cylinders: bool) -> bool {
    let tol = popkern_core::tolerance::Tolerance::default().linear_reach();
    let read = |h: FaceHandle| {
        let b = containing_box(body, h).expect("box");
        let edges = face_boundary(body, h, &[]).edges;
        (h, b, edges)
    };
    let ((hf, bf, ef), (hg, bg, eg)) = (read(f), read(g));
    let s = |h: FaceHandle| &body.faces[h].surface.as_ref().expect("surface").geometry;
    let (side_f, side_g) = (
        (hf, s(hf), &bf, ef.as_slice()),
        (hg, s(hg), &bg, eg.as_slice()),
    );
    if cylinders {
        clear::cylinders_clear(body, side_f, side_g, tol)
    } else {
        clear::plane_crosses_cylinder_clear(body, side_f, side_g, tol)
    }
}

/// A plane face and a cylinder face are cleared only when apart: a
/// cylinder across the plane far outside the face (no normals parallel)
/// and one lying 0.01 above it (the rulings' distances 0.01 and 2.01)
/// are; one through the face, one lying on it, and one along it 1e-7
/// above are not.
#[test]
fn a_plane_and_a_cylinder_are_cleared_only_when_apart() {
    let slab = prism_z(
        &[(-3.0, -3.0), (3.0, -3.0), (3.0, 3.0), (-3.0, 3.0)],
        0.0,
        2.0,
    );
    let along_x = Some(Transform3::rotation_y(std::f64::consts::FRAC_PI_2));
    let cases = [
        (
            "across, far out",
            cylinder_at(1.0, 4.0, None, Vec3::new(10.0, 0.0, -1.0)),
            true,
        ),
        (
            "through the face",
            cylinder_at(1.0, 2.0, None, Vec3::new(0.0, 0.0, 1.0)),
            false,
        ),
        (
            "lying on it",
            cylinder_at(1.0, 4.0, along_x, Vec3::new(-2.0, 0.0, 3.0)),
            false,
        ),
        (
            "lying 0.01 above",
            cylinder_at(1.0, 4.0, along_x, Vec3::new(-2.0, 0.0, 3.01)),
            true,
        ),
        (
            "lying 1e-7 above",
            cylinder_at(1.0, 4.0, along_x, Vec3::new(-2.0, 0.0, 3.0 + 1e-7)),
            false,
        ),
    ];
    for (label, cyl, apart) in cases {
        let body = side_by_side(&[slab.clone(), cyl]);
        let (top, wall) = (
            face_where(&body, top_plane(2.0)),
            face_where(&body, lateral),
        );
        assert_eq!(clears(&body, top, wall, false), apart, "{label}");
        assert!(
            !clears(&body, wall, top, false),
            "{label}: the cylinder is not a plane"
        );
    }
}

/// A cylinder's seam crossing the plane of a triangular face inside the
/// face's box but outside the face is cleared by the face itself; the
/// same cylinder through the face is not.
#[test]
fn a_seam_crossing_the_plane_outside_the_face_but_inside_its_box_is_cleared() {
    let wedge = prism_z(&[(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)], 0.0, 2.0);
    for (at, apart) in [((3.0, 3.0), true), ((1.0, 1.0), false)] {
        let cyl = cylinder_at(0.5, 2.0, None, Vec3::new(at.0, at.1, 1.0));
        let body = side_by_side(&[wedge.clone(), cyl]);
        let (top, wall) = (
            face_where(&body, top_plane(2.0)),
            face_where(&body, lateral),
        );
        assert_eq!(clears(&body, top, wall, false), apart, "cylinder at {at:?}");
    }
}

/// Two cylinders are cleared only when apart: parallel a unit apart and
/// crossed half a unit apart are; parallel or crossed and touching are
/// not.
#[test]
fn two_cylinders_are_cleared_only_when_apart() {
    let along_x = Some(Transform3::rotation_y(std::f64::consts::FRAC_PI_2));
    let upright = cylinder_at(1.0, 4.0, None, Vec3::new(0.0, 0.0, 0.0));
    let cases = [
        (
            "parallel, a unit apart",
            cylinder_at(1.0, 4.0, None, Vec3::new(3.0, 0.0, 0.0)),
            true,
        ),
        (
            "parallel, touching",
            cylinder_at(1.0, 4.0, None, Vec3::new(2.0, 0.0, 0.0)),
            false,
        ),
        (
            "crossed, half apart",
            cylinder_at(1.0, 4.0, along_x, Vec3::new(-2.0, 2.5, 2.0)),
            true,
        ),
        (
            "crossed, touching",
            cylinder_at(1.0, 4.0, along_x, Vec3::new(-2.0, 2.0, 2.0)),
            false,
        ),
    ];
    for (label, other, apart) in cases {
        let body = side_by_side(&[upright.clone(), other]);
        let walls: Vec<FaceHandle> = body
            .faces
            .keys()
            .filter(|&f| lateral(&body.faces[f].surface.as_ref().expect("surface").geometry))
            .collect();
        assert_eq!(walls.len(), 2, "{label}");
        assert_eq!(clears(&body, walls[0], walls[1], true), apart, "{label}");
    }
}

/// The distance from `p` to the triangle `(0, 0)`, `(4, 0)`, `(0, 4)` at
/// height 2.
fn triangle_distance(p: &Point3) -> f64 {
    let (x, y, h) = (p.x, p.y, p.z - 2.0);
    if x >= 0.0 && y >= 0.0 && x + y <= 4.0 {
        return h.abs();
    }
    let corners = [(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)];
    (0..3)
        .map(|k| {
            let ((ax, ay), (bx, by)) = (corners[k], corners[(k + 1) % 3]);
            let (dx, dy) = (bx - ax, by - ay);
            let t = (((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy)).clamp(0.0, 1.0);
            (x - ax - t * dx).hypot(y - ay - t * dy).hypot(h)
        })
        .fold(f64::MAX, f64::min)
}

/// A ball called off a planar face lies farther than the pad from it at
/// every point of a lattice through it, balls above and below the plane
/// included (the boundary is measured from the centre's foot); balls
/// off it and on it are both met.
#[test]
fn a_ball_called_off_a_planar_face_is_off_it() {
    let wedge = prism_z(&[(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)], 0.0, 2.0);
    let top = face_where(&wedge, top_plane(2.0));
    let pad = 0.05;
    let (mut off, mut on) = (0, 0);
    for i in 0..12 {
        for j in 0..12 {
            for k in 0..3 {
                let c = Point3::new(
                    -1.5 + 0.6 * f64::from(i),
                    -1.5 + 0.6 * f64::from(j),
                    2.0 + 0.4 * (f64::from(k) - 1.0),
                );
                let r = 0.2 + 0.15 * f64::from((i + j + k) % 4);
                if !clear::ball_off_face(&wedge, top, c, pad + r) {
                    on += 1;
                    continue;
                }
                off += 1;
                for a in -3..=3 {
                    for b in -3..=3 {
                        for d in -3..=3 {
                            let v = Vec3::new(f64::from(a), f64::from(b), f64::from(d)) * (r / 3.0);
                            if v.length() > r {
                                continue;
                            }
                            let p = c + v;
                            assert!(triangle_distance(&p) > pad, "ball at {c:?} r {r}: {p:?}");
                        }
                    }
                }
            }
        }
    }
    assert!(off > 50 && on > 50, "off {off}, on {on}");
}
