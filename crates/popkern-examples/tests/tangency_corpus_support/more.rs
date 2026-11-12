//! More tangency corpus cases (the file-size rule): the cylinder
//! inscribed through a cube, the dimple family's cone member, a ball
//! resting on a cube's edge, and Viviani's operands with the ball's
//! poles on the curve.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use popkern_bool::BoolOperation;
use popkern_core::geometry::Transform3;
use popkern_core::topology::BRep;

use popkern_core::geometry::Point3;

use super::{Case, Expect, Status, ball, ball_turns, box_at, moved, sphere_volume, turn, viviani};

pub(super) fn more_cases() -> Vec<Case> {
    use BoolOperation::*;
    use Expect::*;
    vec![
        // The cube [-5, 5]³ united with a cylinder of radius 5 on the x
        // axis from x = -10 to 10, tangent to four of its faces along
        // lines: the cube and the cylinder's two ends. Its seam lies on
        // the contact line at z = -5 as built; the turns move it off
        // and onto another.
        Case {
            name: "cylinder-inscribed-through-a-cube",
            op: Union,
            build: |tol| {
                let rod = moved(
                    moved(
                        BRep::cylinder(5.0, 20.0, tol).expect("cylinder"),
                        &Transform3::rotation_y(FRAC_PI_2),
                    ),
                    &Transform3::translation(-10.0, 0.0, 0.0),
                );
                (box_at((-5.0, -5.0, -5.0), (10.0, 10.0, 10.0), tol), rod)
            },
            decline: Body {
                volume: 1000.0 + 250.0 * PI,
                area: 914.159_265_358_979_3,
            },
            general: None,
            turns: vec![
                turn(1, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), FRAC_PI_4),
                turn(1, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), FRAC_PI_2),
                turn(1, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.3),
            ],
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // A conical dimple (rim radius 4 on the plate's top face, apex
        // 4 below) whose rim touches the plate's top edge at (20, 10,
        // 5): the top face touches itself there, a pinch use (F3). Its
        // seam runs through the pinch as built; the turns move it off.
        Case {
            name: "cone-dimple-touching-edge",
            op: Subtraction,
            build: |tol| {
                let cone = moved(
                    moved(
                        BRep::cone(4.0, 4.0, tol).expect("cone"),
                        &Transform3::rotation_x(PI),
                    ),
                    &Transform3::translation(16.0, 10.0, 5.0),
                );
                (box_at((0.0, 0.0, 0.0), (20.0, 20.0, 5.0), tol), cone)
            },
            decline: Body {
                volume: 2000.0 - 64.0 * PI / 3.0,
                // The plate less the dimple's opening, plus the cone's
                // side (slant 4√2).
                area: 1200.0 - 16.0 * PI + 16.0 * 2f64.sqrt() * PI,
            },
            general: None,
            turns: vec![
                turn(1, (16.0, 10.0, 5.0), (0.0, 0.0, 1.0), 0.7),
                turn(1, (16.0, 10.0, 5.0), (0.0, 0.0, 1.0), PI),
            ],
            status: Status::Gate,
            general_status: Status::Gate,
        },
        // A ball of radius 3 resting on the cube's top edge along y at
        // (10, 5, 10): it touches the top and side faces only at their
        // common edge, where the surfaces pass through each other
        // beyond it (an edge grazing a face, R5.8). Its union served
        // the two shells touching at a point.
        Case {
            name: "ball-resting-on-an-edge-union",
            op: Union,
            build: ball_on_edge,
            decline: NonManifold(vec![Point3::new(10.0, 5.0, 10.0)]),
            general: Some(Body {
                volume: 1000.0 + sphere_volume(3.0),
                area: 600.0 + 36.0 * PI,
            }),
            turns: ball_turns(1, ball_on_edge_centre()),
            status: Status::Gate,
            general_status: Status::Gate,
        },
        Case {
            name: "ball-resting-on-an-edge-common",
            op: Intersection,
            build: ball_on_edge,
            decline: Empty,
            general: None,
            turns: ball_turns(1, ball_on_edge_centre()),
            status: Status::Gate,
            general_status: Status::Gate,
        },
    ]
    .into_iter()
    .chain(viviani_poles_on_the_curve())
    .collect()
}

/// Viviani's operands (`viviani`) with the ball as built, its poles at
/// (0, 0, ±10) on the curve and its seam meridian through the double
/// point: each lobe of the figure eight passes through a pole, where
/// the sphere chart's longitude turns half a turn. The scan declined
/// such a branch to the marching path, which ran past two minutes and
/// at kilometre scale served an empty cut (PK-1169). The turns keep
/// the poles on the curve with the seam turned off the double point and
/// to the far side, and turn the rod's seam (the line through both
/// poles) off them.
fn viviani_poles_on_the_curve() -> Vec<Case> {
    use BoolOperation::*;
    use Expect::*;
    let build: fn(popkern_core::tolerance::Tolerance) -> (BRep, BRep) =
        |tol| (BRep::sphere(10.0, tol).expect("sphere"), viviani(tol).1);
    let turns = || {
        vec![
            turn(0, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), FRAC_PI_2),
            turn(0, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), PI),
            turn(1, (5.0, 0.0, 0.0), (0.0, 0.0, 1.0), FRAC_PI_2),
        ]
    };
    let row = |name, op, decline, general| Case {
        name,
        op,
        build,
        decline,
        general,
        turns: turns(),
        status: Status::Gate,
        general_status: Status::Gate,
    };
    vec![
        row(
            "viviani-common-poles-on-the-curve",
            Intersection,
            Body {
                volume: 1_205.506_213_504_306_5,
                area: 628.318_530_717_958_7,
            },
            None,
        ),
        row(
            "viviani-union-poles-on-the-curve",
            Union,
            Body {
                volume: 4_868.239_583_435_96,
                area: 1_539.380_400_258_998_6,
            },
            None,
        ),
        row(
            "viviani-cut-poles-on-the-curve",
            Subtraction,
            NonManifold(vec![Point3::new(10.0, 0.0, 0.0)]),
            Some(Body {
                volume: 2_983.283_991_282_084_6,
                area: 1_428.318_530_717_958_5,
            }),
        ),
    ]
}

/// The centre of a ball of radius 3 touching the edge y-wise at x = z
/// = 10 from outside, off both faces at 45°.
fn ball_on_edge_centre() -> (f64, f64, f64) {
    let d = 3.0 / std::f64::consts::SQRT_2;
    (10.0 + d, 5.0, 10.0 + d)
}

fn ball_on_edge(tol: popkern_core::tolerance::Tolerance) -> (BRep, BRep) {
    (
        box_at((0.0, 0.0, 0.0), (10.0, 10.0, 10.0), tol),
        ball(3.0, ball_on_edge_centre(), tol),
    )
}
