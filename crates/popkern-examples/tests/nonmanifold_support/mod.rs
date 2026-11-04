//! The non-manifold corpus cases (PK-967), shared by the decline and
//! split-policy corpus tests.

#![allow(dead_code)]

use popkern_bool::{BoolError, BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Expect {
    /// Served valid with these shells and this volume.
    Valid {
        shells: usize,
        volume: f64,
    },
    NonManifold,
    Empty,
}

#[derive(Debug, Clone, Copy)]
pub enum Shape {
    /// Axis-aligned cube: min corner and edge length.
    Cube((f64, f64, f64), f64),
    /// Cylinder along +x: radius, length, and the centre of its start
    /// cap.
    RodX(f64, f64, (f64, f64, f64)),
    /// Unit cube turned 45° about x, its lowest edge along x at the
    /// given (y, z).
    EdgeDownX(f64, f64),
    /// Sphere: radius and centre.
    Ball(f64, (f64, f64, f64)),
}

pub fn build(shape: Shape, place: &Transform3, tol: Tolerance) -> BRep {
    let local = match shape {
        Shape::Cube(lo, s) => BRep::unit_cube(tol)
            .expect("cube")
            .transform(&Transform3::scale_uniform(s))
            .expect("scale")
            .transform(&Transform3::translation(lo.0, lo.1, lo.2))
            .expect("translate"),
        Shape::RodX(r, len, c) => BRep::cylinder(r, len, tol)
            .expect("cylinder")
            .transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2))
            .expect("rotate")
            .transform(&Transform3::translation(c.0, c.1, c.2))
            .expect("translate"),
        Shape::EdgeDownX(y, z) => BRep::unit_cube(tol)
            .expect("cube")
            .transform(&Transform3::translation(0.0, -0.5, -0.5))
            .expect("centre")
            .transform(&Transform3::rotation_x(std::f64::consts::FRAC_PI_4))
            .expect("rotate")
            .transform(&Transform3::translation(
                0.0,
                y,
                z + std::f64::consts::FRAC_1_SQRT_2,
            ))
            .expect("translate"),
        Shape::Ball(r, c) => BRep::sphere(r, tol)
            .expect("sphere")
            .transform(&Transform3::translation(c.0, c.1, c.2))
            .expect("translate"),
    };
    local.transform(place).expect("place")
}

pub struct Case {
    pub name: &'static str,
    pub a: Shape,
    pub b: Shape,
    pub op: BoolOperation,
    pub expect: Expect,
}

pub fn cases() -> Vec<Case> {
    use BoolOperation::*;
    use Expect::*;
    use Shape::*;
    let unit = Cube((0.0, 0.0, 0.0), 1.0);
    let big = Cube((0.0, 0.0, 0.0), 2.0);
    vec![
        // Edge contact.
        Case {
            name: "edge_touch_union",
            a: unit,
            b: Cube((1.0, 1.0, 0.0), 1.0),
            op: Union,
            expect: NonManifold,
        },
        Case {
            name: "edge_touch_intersection",
            a: unit,
            b: Cube((1.0, 1.0, 0.0), 1.0),
            op: Intersection,
            expect: Empty,
        },
        Case {
            name: "edge_touch_subtraction",
            a: unit,
            b: Cube((1.0, 1.0, 0.0), 1.0),
            op: Subtraction,
            expect: Valid {
                shells: 1,
                volume: 1.0,
            },
        },
        // Vertex contact.
        Case {
            name: "vertex_touch_union",
            a: unit,
            b: Cube((1.0, 1.0, 1.0), 1.0),
            op: Union,
            expect: NonManifold,
        },
        Case {
            name: "vertex_touch_intersection",
            a: unit,
            b: Cube((1.0, 1.0, 1.0), 1.0),
            op: Intersection,
            expect: Empty,
        },
        Case {
            name: "vertex_touch_subtraction",
            a: unit,
            b: Cube((1.0, 1.0, 1.0), 1.0),
            op: Subtraction,
            expect: Valid {
                shells: 1,
                volume: 1.0,
            },
        },
        // Face contact: full, partial.
        Case {
            name: "face_flush_union",
            a: unit,
            b: Cube((1.0, 0.0, 0.0), 1.0),
            op: Union,
            expect: Valid {
                shells: 1,
                volume: 2.0,
            },
        },
        Case {
            name: "face_partial_union",
            a: unit,
            b: Cube((1.0, 0.5, 0.5), 1.0),
            op: Union,
            expect: Valid {
                shells: 1,
                volume: 2.0,
            },
        },
        Case {
            name: "face_partial_intersection",
            a: unit,
            b: Cube((1.0, 0.5, 0.5), 1.0),
            op: Intersection,
            expect: Empty,
        },
        // Contained tool flush on 1, 2, 3 faces (PK-965), and inside.
        Case {
            name: "flush3_notch",
            a: big,
            b: Cube((1.0, 0.0, 0.0), 1.0),
            op: Subtraction,
            expect: Valid {
                shells: 1,
                volume: 7.0,
            },
        },
        Case {
            name: "flush2_notch",
            a: big,
            b: Cube((1.0, 0.0, 0.5), 1.0),
            op: Subtraction,
            expect: Valid {
                shells: 1,
                volume: 7.0,
            },
        },
        Case {
            name: "flush1_pocket",
            a: big,
            b: Cube((0.5, 0.5, 1.0), 1.0),
            op: Subtraction,
            expect: Valid {
                shells: 1,
                volume: 7.0,
            },
        },
        Case {
            name: "flush0_cavity",
            a: big,
            b: Cube((0.5, 0.5, 0.5), 1.0),
            op: Subtraction,
            expect: Valid {
                shells: 2,
                volume: 7.0,
            },
        },
        Case {
            name: "flush3_union",
            a: big,
            b: Cube((1.0, 0.0, 0.0), 1.0),
            op: Union,
            expect: Valid {
                shells: 1,
                volume: 8.0,
            },
        },
        Case {
            name: "flush3_intersection",
            a: big,
            b: Cube((1.0, 0.0, 0.0), 1.0),
            op: Intersection,
            expect: Valid {
                shells: 1,
                volume: 1.0,
            },
        },
        // Overlap leaving a manifold result (control).
        Case {
            name: "overlap_union",
            a: unit,
            b: Cube((0.5, 0.5, 0.5), 1.0),
            op: Union,
            expect: Valid {
                shells: 1,
                volume: 1.875,
            },
        },
        Case {
            name: "overlap_subtraction",
            a: unit,
            b: Cube((0.5, 0.5, 0.5), 1.0),
            op: Subtraction,
            expect: Valid {
                shells: 1,
                volume: 0.875,
            },
        },
        // Line contact: a rod resting on a box top.
        Case {
            name: "line_contact_intersection",
            a: unit,
            b: RodX(0.5, 1.0, (0.0, 0.5, 1.5)),
            op: Intersection,
            expect: Empty,
        },
        Case {
            name: "line_contact_union",
            a: unit,
            b: RodX(0.5, 1.0, (0.0, 0.5, 1.5)),
            op: Union,
            expect: NonManifold,
        },
        // Contacts sharing no topology: an edge lying on a face, a ball
        // resting on a face.
        Case {
            name: "edge_on_face_union",
            a: unit,
            b: EdgeDownX(0.5, 1.0),
            op: Union,
            expect: NonManifold,
        },
        Case {
            name: "edge_on_face_intersection",
            a: unit,
            b: EdgeDownX(0.5, 1.0),
            op: Intersection,
            expect: Empty,
        },
        Case {
            name: "ball_on_face_union",
            a: unit,
            b: Ball(0.4, (0.5, 0.5, 1.4)),
            op: Union,
            expect: NonManifold,
        },
        Case {
            name: "ball_on_face_subtraction",
            a: unit,
            b: Ball(0.4, (0.5, 0.5, 1.4)),
            op: Subtraction,
            expect: Valid {
                shells: 1,
                volume: 1.0,
            },
        },
    ]
}

/// Two-step case: two overshooting cuts from a block leave two blocks
/// joined only along an edge — a non-manifold edge.
pub fn two_cut_nm_edge(place: &Transform3, tol: Tolerance) -> Result<BRep, BoolError> {
    let s1 = boolean(
        BoolOperation::Subtraction,
        build(Shape::Cube((0.0, 0.0, 0.0), 2.0), place, tol),
        build(Shape::Cube((1.0, -0.5, -0.5), 1.5), place, tol),
        tol,
        None,
    )?;
    boolean(
        BoolOperation::Subtraction,
        s1,
        build(Shape::Cube((-0.5, 1.0, -0.5), 1.5), place, tol),
        tol,
        None,
    )
}

pub fn placements() -> Vec<(&'static str, Transform3, f64)> {
    let rot = Transform3::rotation_around_axis(
        Point3::new(0.3, -0.2, 0.1),
        Vec3::new(1.0, 2.0, 3.0),
        0.37,
    )
    .expect("rotation");
    vec![
        ("identity", Transform3::identity(), 1.0),
        ("rotated", rot, 1.0),
        ("mm", Transform3::scale_uniform(1.0e-2), 1.0e-6),
        ("far", Transform3::translation(1.0e4, -2.0e4, 3.0e4), 1.0),
    ]
}
