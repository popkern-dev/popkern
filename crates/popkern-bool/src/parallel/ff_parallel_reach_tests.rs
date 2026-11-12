use super::*;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{Point3, Vec3};

/// The nine-point rational quadratic circle of radius `r` about the z
/// axis at height `z`, centred at `(cx, cy)`.
fn circle(cx: f64, cy: f64, z: f64, r: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let corners = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
    ];
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: corners
            .iter()
            .map(|(x, y)| Point3::new(cx + r * x, cy + r * y, z))
            .collect(),
        weights: (0..9).map(|i| if i % 2 == 1 { w } else { 1.0 }).collect(),
    }
}

fn segment(p: Point3, q: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p, q],
        weights: vec![1.0, 1.0],
    }
}

fn plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.4),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn boss_bottom() -> Aabb {
    Aabb {
        min: Point3::new(-2.0, -2.0, 1.4),
        max: Point3::new(2.0, 2.0, 1.4),
    }
}

/// The distance from `p` to box `b`.
fn point_gap(b: &Aabb, p: &Point3) -> f64 {
    let g = |lo: f64, hi: f64, v: f64| (lo - v).max(v - hi).max(0.0);
    let (dx, dy, dz) = (
        g(b.min.x, b.max.x, p.x),
        g(b.min.y, b.max.y, p.y),
        g(b.min.z, b.max.z, p.z),
    );
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// The disc's wall circle round the boss's bottom face is clear of it,
/// though its control box holds the face: the halving proves it.
#[test]
fn a_circle_round_a_face_is_clear_of_it() {
    let pad = 1e-4;
    assert!(curve_clear_of_face(
        &circle(0.0, 0.0, 1.4, 30.0),
        &plane(),
        &boss_bottom(),
        pad
    ));
    // A circle through the face, or passing within the pad, is not.
    assert!(!curve_clear_of_face(
        &circle(0.0, 0.0, 1.4, 2.0),
        &plane(),
        &boss_bottom(),
        pad
    ));
    assert!(!curve_clear_of_face(
        &circle(0.0, 0.0, 1.4 + 0.5 * pad, 30.0),
        &plane(),
        &Aabb {
            min: Point3::new(-2.0, 29.0, 1.4),
            max: Point3::new(2.0, 31.0, 1.4),
        },
        pad,
    ));
}

/// A segment crossing the face, or ending within the pad, is not clear;
/// one beyond it is.
#[test]
fn segments_are_clear_only_beyond_the_pad() {
    let pad = 1e-4;
    let b = boss_bottom();
    let s = plane();
    let crossing = segment(Point3::new(-5.0, 0.0, 1.4), Point3::new(5.0, 0.0, 1.4));
    assert!(!curve_clear_of_face(&crossing, &s, &b, pad));
    let near = segment(
        Point3::new(2.0 + 0.5 * pad, -5.0, 1.4),
        Point3::new(2.0 + 0.5 * pad, 5.0, 1.4),
    );
    assert!(!curve_clear_of_face(&near, &s, &b, pad));
    let far = segment(
        Point3::new(2.0 + 2.0 * pad, -5.0, 1.4),
        Point3::new(2.0 + 2.0 * pad, 5.0, 1.4),
    );
    assert!(curve_clear_of_face(&far, &s, &b, pad));
}

/// A face on a surface whose box may fall short of it, and a curve of a
/// weight not positive, are never cleared.
#[test]
fn only_contained_boxes_and_positive_weights_clear() {
    let pad = 1e-4;
    let far = circle(0.0, 0.0, 1.4, 30.0);
    let nurbs = SurfaceGeometry::Nurbs(popkern_core::geometry::nurbs_surface::NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(-2.0, -2.0, 1.4),
            Point3::new(-2.0, 2.0, 1.4),
            Point3::new(2.0, -2.0, 1.4),
            Point3::new(2.0, 2.0, 1.4),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    assert!(!curve_clear_of_face(&far, &nurbs, &boss_bottom(), pad));
    let mut negative = far.clone();
    negative.weights[1] = -0.5;
    assert!(!curve_clear_of_face(
        &negative,
        &plane(),
        &boss_bottom(),
        pad
    ));
}

/// Whatever is cleared lies farther than the pad from the box at every
/// sample: circles and segments about boxes over a sweep of placements.
#[test]
fn a_cleared_curve_never_comes_within_the_pad() {
    let pad = 0.05;
    let mut cleared = 0;
    let mut kept = 0;
    for i in 0..12 {
        for j in 0..12 {
            let (cx, cy) = (-3.0 + 0.5 * f64::from(i), -3.0 + 0.5 * f64::from(j));
            let b = Aabb {
                min: Point3::new(-1.0, -0.5, 0.0),
                max: Point3::new(1.0, 0.5, 0.2),
            };
            let curves = [
                circle(cx, cy, 0.1, 1.5),
                circle(cx, cy, 0.3, 0.7),
                segment(Point3::new(cx, cy, 0.1), Point3::new(cy, -cx, 0.25)),
            ];
            for c in &curves {
                if !curve_clear_of_face(c, &plane(), &b, pad) {
                    kept += 1;
                    continue;
                }
                cleared += 1;
                let (t0, t1) = c.parameter_range();
                for k in 0..=2000 {
                    let p = c.point_at(t0 + (t1 - t0) * f64::from(k) / 2000.0);
                    assert!(point_gap(&b, &p) > pad, "{c:?} at {p:?} within {pad}");
                }
            }
        }
    }
    assert!(cleared > 50 && kept > 50, "cleared {cleared}, kept {kept}");
}

/// The distance from `p` to the whole of a plane, sphere, cylinder, cone
/// or torus.
fn surface_distance(s: &SurfaceGeometry, p: &Point3) -> f64 {
    match s {
        SurfaceGeometry::Plane { origin, normal, .. } => (*p - *origin)
            .dot(&normal.normalized().expect("normal"))
            .abs(),
        SurfaceGeometry::Sphere { center, radius, .. } => (p.distance_to(center) - radius).abs(),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let a = axis.normalized().expect("axis");
            let d = *p - *origin;
            ((d - a * d.dot(&a)).length() - radius).abs()
        }
        // The nearer of the two generator rays in the half-plane through
        // the axis, each foot clamped to its ray.
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let a = axis.normalized().expect("axis");
            let v = *p - *apex;
            let (h, rho) = (v.dot(&a), (v - a * v.dot(&a)).length());
            [1.0, -1.0]
                .iter()
                .map(|side: &f64| {
                    let (dh, dr) = (side * half_angle.cos(), half_angle.sin());
                    let t = (h * dh + rho * dr).max(0.0);
                    (h - t * dh).hypot(rho - t * dr)
                })
                .fold(f64::MAX, f64::min)
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let a = axis.normalized().expect("axis");
            let v = *p - *center;
            let (h, rho) = (v.dot(&a), (v - a * v.dot(&a)).length());
            ((rho - major_radius).hypot(h) - minor_radius).abs()
        }
        _ => unreachable!(),
    }
}

/// Whatever box is cleared of a surface lies farther than the pad from
/// it at every point of a dense lattice; boxes inside, outside and
/// across a turned plane, sphere, cylinder, cone (both nappes) and
/// torus are each cleared and kept.
#[test]
fn a_box_cleared_of_a_surface_never_comes_within_the_pad() {
    let pad = 0.05;
    let axis = Vec3::new(0.3, -0.5, 0.8).normalized().expect("axis");
    let o = Point3::new(0.4, -1.0, 2.0);
    let surfaces = [
        SurfaceGeometry::Plane {
            origin: o,
            normal: axis,
            u_axis: Vec3::new(0.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 0.0),
        },
        SurfaceGeometry::Sphere {
            center: o,
            radius: 2.0,
            axis,
        },
        SurfaceGeometry::Cylinder {
            origin: o,
            axis,
            radius: 2.0,
        },
        SurfaceGeometry::Cone {
            apex: o,
            axis,
            half_angle: 0.6,
        },
        SurfaceGeometry::Torus {
            center: o,
            axis,
            major_radius: 2.5,
            minor_radius: 0.8,
        },
    ];
    for s in &surfaces {
        let (mut cleared, mut kept) = (0, 0);
        for i in 0..9 {
            for j in 0..9 {
                for k in 0..9 {
                    let lo = Point3::new(
                        -4.0 + f64::from(i),
                        -5.0 + f64::from(j),
                        -2.0 + f64::from(k),
                    );
                    let size = 0.3 + 0.2 * f64::from((i + j + k) % 4);
                    let b = Aabb {
                        min: lo,
                        max: Point3::new(lo.x + size, lo.y + 0.5 * size, lo.z + size),
                    };
                    if !box_clear_of_surface(&b, s, pad) {
                        kept += 1;
                        continue;
                    }
                    cleared += 1;
                    let n = 8;
                    for a in 0..=n {
                        for c in 0..=n {
                            for d in 0..=n {
                                let f = |lo: f64, hi: f64, t: i32| {
                                    lo + (hi - lo) * f64::from(t) / f64::from(n)
                                };
                                let p = Point3::new(
                                    f(b.min.x, b.max.x, a),
                                    f(b.min.y, b.max.y, c),
                                    f(b.min.z, b.max.z, d),
                                );
                                assert!(surface_distance(s, &p) > pad, "{s:?} {b:?} at {p:?}");
                            }
                        }
                    }
                }
            }
        }
        assert!(
            cleared > 50 && kept > 20,
            "{s:?}: cleared {cleared}, kept {kept}"
        );
    }
}

/// The boss's bottom face against the disc's wall: apart, by either
/// face's box; a face whose box may fall short of it decides nothing.
#[test]
fn the_boss_bottom_and_the_disc_wall_are_apart() {
    let wall = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 30.0,
    };
    let wall_box = Aabb {
        min: Point3::new(-30.0, -30.0, 0.0),
        max: Point3::new(30.0, 30.0, 1.5),
    };
    assert!(faces_apart(
        (&plane(), &boss_bottom()),
        (&wall, &wall_box),
        1e-4
    ));
    assert!(faces_apart(
        (&wall, &wall_box),
        (&plane(), &boss_bottom()),
        1e-4
    ));
    // The disc's top plane at z = 1.5 meets the wall: not apart.
    let top = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.5),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert!(!faces_apart((&top, &wall_box), (&wall, &wall_box), 1e-4));
    let nurbs = SurfaceGeometry::Nurbs(popkern_core::geometry::nurbs_surface::NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(-2.0, -2.0, 1.4),
            Point3::new(-2.0, 2.0, 1.4),
            Point3::new(2.0, -2.0, 1.4),
            Point3::new(2.0, 2.0, 1.4),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    // The NURBS face's box is not trusted, and the wall's box reaches
    // the NURBS surface's extension: not decided.
    assert!(!faces_apart(
        (&nurbs, &boss_bottom()),
        (&wall, &wall_box),
        1e-4
    ));
}
