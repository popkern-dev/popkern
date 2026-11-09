//! More tangency corpus cases (the file-size rule): the cylinder
//! inscribed through a cube, and the dimple family's cone member.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use popkern_bool::BoolOperation;
use popkern_core::geometry::Transform3;
use popkern_core::topology::BRep;

use super::{Case, Expect, Status, box_at, moved, turn};

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
            status: Status::Target("PK-1144"),
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
            status: Status::Target("PK-1143"),
            general_status: Status::Gate,
        },
    ]
}
