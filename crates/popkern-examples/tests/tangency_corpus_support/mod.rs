//! The tangency corpus cases (PK-1124): operands that touch, are
//! tangent, or meet in curves that cross, each with its closed-form
//! answer per non-manifold policy and the self-turns that move its
//! operands' seams and poles without changing the answer.

#![allow(dead_code)]

use std::f64::consts::{FRAC_PI_2, PI};

use popkern_bool::BoolOperation;
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

mod more;
mod one_shell;

/// What a boolean should give.
#[derive(Debug, Clone, PartialEq)]
pub enum Expect {
    /// A valid served body with this volume and area.
    Body { volume: f64, area: f64 },
    /// A decline with `NonManifoldResult` naming exactly these points.
    NonManifold(Vec<Point3>),
    /// An empty result.
    Empty,
}

/// Whether a case's row must pass now (a gate) or is the target of a
/// later stage (recorded, and must not pass silently: a target that
/// starts passing is promoted).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Gate,
    /// The stage whose exit gate this row is.
    Target(&'static str),
}

/// A turn of one operand about a symmetry of it: the answer is
/// unchanged, but its seam and poles move.
#[derive(Debug, Clone, Copy)]
pub struct SelfTurn {
    pub operand: usize,
    pub origin: (f64, f64, f64),
    pub axis: (f64, f64, f64),
    pub angle: f64,
}

#[derive(Clone)]
pub struct Case {
    pub name: &'static str,
    pub op: BoolOperation,
    pub build: fn(Tolerance) -> (BRep, BRep),
    /// Under `NonManifoldPolicy::Decline`.
    pub decline: Expect,
    /// Under `NonManifoldPolicy::General`, when it differs.
    pub general: Option<Expect>,
    pub turns: Vec<SelfTurn>,
    /// The `decline` row's status.
    pub status: Status,
    /// The `general` row's status.
    pub general_status: Status,
}

fn moved(b: BRep, t: &Transform3) -> BRep {
    b.transform(t).expect("transform")
}

/// The box from `lo` with edge lengths `size`, as a frustum of two
/// equal rectangles (a scaled unit cube would need a non-uniform scale,
/// which `BRep::transform` refuses).
fn box_at(lo: (f64, f64, f64), size: (f64, f64, f64), tol: Tolerance) -> BRep {
    let (x0, y0, z0) = lo;
    let (x1, y1, z1) = (x0 + size.0, y0 + size.1, z0 + size.2);
    let rect = |z: f64| {
        vec![
            Point3::new(x0, y0, z),
            Point3::new(x1, y0, z),
            Point3::new(x1, y1, z),
            Point3::new(x0, y1, z),
        ]
    };
    BRep::polygon_frustum(&rect(z0), &rect(z1), tol).expect("box")
}

fn ball(r: f64, c: (f64, f64, f64), tol: Tolerance) -> BRep {
    moved(
        BRep::sphere(r, tol).expect("sphere"),
        &Transform3::translation(c.0, c.1, c.2),
    )
}

/// A half-space as a large block: the side of the plane through `p`
/// with outward normal `n` that the normal points away from, reaching
/// `reach` past `p` in every direction.
fn half_space(p: (f64, f64, f64), n: (f64, f64, f64), reach: f64, tol: Tolerance) -> BRep {
    let block = box_at(
        (-reach, -reach, -2.0 * reach),
        (2.0 * reach, 2.0 * reach, 2.0 * reach),
        tol,
    );
    // Turn +z onto n, then move the top face's centre to p.
    let n = Vec3::new(n.0, n.1, n.2).normalized().expect("normal");
    let z = Vec3::new(0.0, 0.0, 1.0);
    let axis = z.cross(&n);
    let turned = if axis.length() < 1e-12 {
        if n.z > 0.0 {
            block
        } else {
            moved(block, &Transform3::rotation_x(PI))
        }
    } else {
        let angle = z.dot(&n).clamp(-1.0, 1.0).acos();
        moved(
            block,
            &Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), axis, angle)
                .expect("axis"),
        )
    };
    moved(turned, &Transform3::translation(p.0, p.1, p.2))
}

fn turn(operand: usize, origin: (f64, f64, f64), axis: (f64, f64, f64), angle: f64) -> SelfTurn {
    SelfTurn {
        operand,
        origin,
        axis,
        angle,
    }
}

/// Turns of a ball about its centre: off every contact.
fn ball_turns(operand: usize, c: (f64, f64, f64)) -> Vec<SelfTurn> {
    vec![
        turn(operand, c, (1.0, 0.7, 0.3), 0.645_771_823_237_9),
        turn(operand, c, (0.0, 1.0, 0.0), FRAC_PI_2),
        turn(operand, c, (1.0, 0.0, 0.0), FRAC_PI_2),
    ]
}

fn sphere_volume(r: f64) -> f64 {
    4.0 / 3.0 * PI * r * r * r
}

pub fn cases() -> Vec<Case> {
    use BoolOperation::*;
    use Expect::*;
    let face_centres = vec![
        Point3::new(5.0, 5.0, 0.0),
        Point3::new(5.0, 5.0, 10.0),
        Point3::new(5.0, 0.0, 5.0),
        Point3::new(5.0, 10.0, 5.0),
        Point3::new(0.0, 5.0, 5.0),
        Point3::new(10.0, 5.0, 5.0),
    ];
    let inscribed = |tol| {
        (
            box_at((0.0, 0.0, 0.0), (10.0, 10.0, 10.0), tol),
            ball(5.0, (5.0, 5.0, 5.0), tol),
        )
    };
    vec![
        // The cube [0, 10]³ and its inscribed ball, touching each face
        // at its centre from inside.
        Case {
            name: "inscribed-ball-common",
            op: Intersection,
            build: inscribed,
            decline: Body {
                volume: sphere_volume(5.0),
                area: 100.0 * PI,
            },
            general: None,
            turns: ball_turns(1, (5.0, 5.0, 5.0)),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        Case {
            name: "inscribed-ball-union",
            op: Union,
            build: inscribed,
            decline: Body {
                volume: 1000.0,
                area: 600.0,
            },
            general: None,
            turns: ball_turns(1, (5.0, 5.0, 5.0)),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // The cavity touches the outside at the six face centres: a
        // wall of no thickness there.
        Case {
            name: "inscribed-ball-cut",
            op: Subtraction,
            build: inscribed,
            decline: NonManifold(face_centres),
            general: Some(Body {
                volume: 1000.0 - sphere_volume(5.0),
                area: 600.0 + 100.0 * PI,
            }),
            turns: ball_turns(1, (5.0, 5.0, 5.0)),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // A ball resting on the cube's top face from outside.
        Case {
            name: "ball-resting-common",
            op: Intersection,
            build: |tol| {
                (
                    box_at((0.0, 0.0, 0.0), (10.0, 10.0, 10.0), tol),
                    ball(5.0, (5.0, 5.0, 15.0), tol),
                )
            },
            decline: Empty,
            general: None,
            turns: ball_turns(1, (5.0, 5.0, 15.0)),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        Case {
            name: "ball-resting-union",
            op: Union,
            build: |tol| {
                (
                    box_at((0.0, 0.0, 0.0), (10.0, 10.0, 10.0), tol),
                    ball(5.0, (5.0, 5.0, 15.0), tol),
                )
            },
            decline: NonManifold(vec![Point3::new(5.0, 5.0, 10.0)]),
            general: Some(Body {
                volume: 1000.0 + sphere_volume(5.0),
                area: 600.0 + 100.0 * PI,
            }),
            turns: ball_turns(1, (5.0, 5.0, 15.0)),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // The cube [−5, 5]³ and the ball of radius 5√2 touching its
        // twelve edges at their midpoints: six disks and eight sphere
        // patches, manifold.
        Case {
            name: "midsphere-common",
            op: Intersection,
            build: |tol| {
                (
                    box_at((-5.0, -5.0, -5.0), (10.0, 10.0, 10.0), tol),
                    ball(5.0 * 2f64.sqrt(), (0.0, 0.0, 0.0), tol),
                )
            },
            decline: Body {
                volume: 965.068_858_214_997_2,
                area: 547.466_718_050_061_6,
            },
            general: None,
            turns: ball_turns(1, (0.0, 0.0, 0.0)),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // A ball of radius 10 and a cylinder of radius 5 on x = 5,
        // y = 0, |z| ≤ 12, internally tangent at (10, 0, 0): Viviani's
        // figure eight crosses itself there.
        Case {
            name: "viviani-common",
            op: Intersection,
            build: viviani,
            decline: Body {
                volume: 1_205.506_213_504_306_5,
                area: 628.318_530_717_958_7,
            },
            general: None,
            turns: viviani_turns(),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        Case {
            name: "viviani-union",
            op: Union,
            build: viviani,
            decline: Body {
                volume: 4_868.239_583_435_96,
                area: 1_539.380_400_258_998_6,
            },
            general: None,
            turns: viviani_turns(),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        Case {
            name: "viviani-cut",
            op: Subtraction,
            build: viviani,
            decline: NonManifold(vec![Point3::new(10.0, 0.0, 0.0)]),
            general: Some(Body {
                volume: 2_983.283_991_282_084_6,
                area: 1_428.318_530_717_958_5,
            }),
            turns: viviani_turns(),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // A torus (R 5, r 3) cut by the bitangent plane 3x − 4z = 0,
        // keeping z ≤ 3x/4: Villarceau's two circles cross at the two
        // tangent points; at (−16/5, 0, −12/5) the result is pinched.
        Case {
            name: "villarceau-half-torus",
            op: Intersection,
            build: |tol| {
                let torus = moved(
                    BRep::torus(5.0, 3.0, tol).expect("torus"),
                    &Transform3::rotation_z(FRAC_PI_2),
                );
                (
                    torus,
                    half_space((0.0, 0.0, 0.0), (-3.0, 0.0, 4.0), 20.0, tol),
                )
            },
            decline: NonManifold(vec![Point3::new(-3.2, 0.0, -2.4)]),
            general: Some(Body {
                volume: 444.132_198_049_021_16,
                area: 408.438_242_912_009_2,
            }),
            turns: vec![turn(0, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.4)],
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // A cone of radius 12 tangent to a ball of radius 15 along its
        // base rim: the ice-cream cone.
        Case {
            name: "ice-cream-cone",
            op: Union,
            build: |tol| {
                let cone = moved(
                    moved(
                        BRep::cone(12.0, 16.0, tol).expect("cone"),
                        &Transform3::rotation_x(PI),
                    ),
                    &Transform3::translation(0.0, 0.0, -9.0),
                );
                (ball(15.0, (0.0, 0.0, 0.0), tol), cone)
            },
            decline: Body {
                volume: 15_079.644_737_231_007,
                area: 3_015.928_947_446_201,
            },
            general: None,
            turns: vec![
                turn(0, (0.0, 0.0, 0.0), (1.0, 0.7, 0.3), 0.645_771_823_237_9),
                turn(1, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.7),
            ],
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // A spherical dimple whose rim touches the plate's top edge at
        // (20, 10, 5): the top face touches itself there, the solid is
        // one wedge (a pinch use, F3).
        Case {
            name: "dimple-touching-edge",
            op: Subtraction,
            build: |tol| {
                (
                    box_at((0.0, 0.0, 0.0), (20.0, 20.0, 5.0), tol),
                    ball(4.0, (16.0, 10.0, 5.0), tol),
                )
            },
            decline: Body {
                volume: 2000.0 - 2.0 / 3.0 * PI * 64.0,
                // The plate's area less the dimple's opening, plus the
                // hemisphere.
                area: 2.0 * 400.0 + 4.0 * 100.0 - 16.0 * PI + 32.0 * PI,
            },
            general: None,
            turns: ball_turns(1, (16.0, 10.0, 5.0)),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // A ball halved through its poles (a singular case that serves).
        Case {
            name: "ball-halved-through-its-poles",
            op: Subtraction,
            build: |tol| {
                (
                    ball(5.0, (0.0, 0.0, 0.0), tol),
                    half_space((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 20.0, tol),
                )
            },
            decline: Body {
                volume: 261.799_387_799_149_43,
                area: 235.619_449_019_234_48,
            },
            general: None,
            turns: vec![turn(0, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.5)],
            status: Status::Gate,
            general_status: Status::Gate,
        },
    ]
    .into_iter()
    .chain(more::more_cases())
    .chain(one_shell::one_shell_cases())
    .collect()
}

/// Viviani's operands, the ball's poles on the y axis and its seam
/// through (−10, 0, 0), the cylinder's seam on its side at x = 0.
fn viviani(tol: Tolerance) -> (BRep, BRep) {
    // The sphere's poles are on z and its seam at +x: turn z onto y
    // (poles at ±y), then a half turn about y (seam through −x).
    let ball = moved(
        moved(
            BRep::sphere(10.0, tol).expect("sphere"),
            &Transform3::rotation_x(-FRAC_PI_2),
        ),
        &Transform3::rotation_y(PI),
    );
    // The cylinder runs z ∈ [0, 24] with its seam at +x: a half turn
    // about its axis puts the seam at −x, then move it to x = 5.
    let rod = moved(
        moved(
            BRep::cylinder(5.0, 24.0, tol).expect("cylinder"),
            &Transform3::rotation_z(PI),
        ),
        &Transform3::translation(5.0, 0.0, -12.0),
    );
    (ball, rod)
}

fn viviani_turns() -> Vec<SelfTurn> {
    vec![
        turn(0, (0.0, 0.0, 0.0), (1.0, 0.7, 0.3), 0.645_771_823_237_9),
        turn(0, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), FRAC_PI_2),
        turn(1, (5.0, 0.0, 0.0), (0.0, 0.0, 1.0), PI),
    ]
}

/// The scene placements every case runs under: identity, two rigid
/// turns, a far translation and two scales.
pub fn placements() -> Vec<(&'static str, Transform3, f64)> {
    let t = |p: Transform3, s: f64| (p, s);
    let turned = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.7, 0.3),
        0.645_771_823_237_9,
    )
    .expect("axis");
    vec![
        ("identity", t(Transform3::identity(), 1.0)),
        ("turned", t(turned, 1.0)),
        ("quarter-z", t(Transform3::rotation_z(FRAC_PI_2), 1.0)),
        (
            "far",
            t(Transform3::translation(1234.5, -678.25, 910.125), 1.0),
        ),
        ("milli", t(Transform3::scale_uniform(1e-3), 1e-3)),
        ("kilo", t(Transform3::scale_uniform(1e3), 1e3)),
    ]
    .into_iter()
    .map(|(n, (p, s))| (n, p, s))
    .collect()
}

/// Apply a self-turn to operand `which` of `(a, b)`.
pub fn apply_turn(pair: (BRep, BRep), t: &SelfTurn) -> (BRep, BRep) {
    let r = Transform3::rotation_around_axis(
        Point3::new(t.origin.0, t.origin.1, t.origin.2),
        Vec3::new(t.axis.0, t.axis.1, t.axis.2),
        t.angle,
    )
    .expect("axis");
    match t.operand {
        0 => (moved(pair.0, &r), pair.1),
        _ => (pair.0, moved(pair.1, &r)),
    }
}
