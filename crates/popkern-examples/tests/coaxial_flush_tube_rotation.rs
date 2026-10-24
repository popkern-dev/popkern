//! PK-932: a cylinder minus a coaxial cylinder of the same height (a
//! tube with flush caps) under rigid rotations of both operands.
//!
//! The inner wall's rims are the cap x wall records themselves. Under a
//! general rotation the wall's chart seam moved to u 2.80 while the
//! bottom rim record ran u 0..2π: the clip's boundary-riding test
//! compared the record against the rim on one period branch only, read
//! the seam as an exit, kept half the rim, and the serve declined on a
//! curveless edge. The test now compares on every period branch.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::body_volume;

#[test]
fn a_flush_tube_serves_its_volume_under_rotation() {
    let tol = Tolerance::default();
    let want = std::f64::consts::PI * (1.0 - 0.25);
    let rotations = [
        ("unrotated", Transform3::identity()),
        ("z 0.7", Transform3::rotation_z(0.7)),
        ("x 1.1", Transform3::rotation_x(1.1)),
        ("y -0.4", Transform3::rotation_y(-0.4)),
        (
            "mixed",
            Transform3::rotation_x(0.3)
                .compose(&Transform3::rotation_y(1.9))
                .compose(&Transform3::rotation_z(-2.2)),
        ),
        (
            "mixed 2",
            Transform3::rotation_x(1.3).compose(&Transform3::rotation_z(0.4)),
        ),
        (
            "mixed 3",
            Transform3::rotation_y(2.5).compose(&Transform3::rotation_x(-0.8)),
        ),
    ];
    for (name, r) in rotations {
        let outer = BRep::cylinder(1.0, 1.0, tol)
            .and_then(|c| c.transform(&r))
            .expect("outer");
        let inner = BRep::cylinder(0.5, 1.0, tol)
            .and_then(|c| c.transform(&r))
            .expect("inner");
        let body = boolean(BoolOperation::Subtraction, outer, inner, tol, None)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let errors = body.validate_all();
        assert!(errors.is_empty(), "{name}: {errors:?}");
        assert_eq!(body.faces.len(), 4, "{name}: two walls, two annuli");
        let v = body_volume(&body, &tol).expect("volume");
        assert!(
            (v - want).abs() < 1e-9 * want,
            "{name}: volume {v} vs {want}"
        );
    }
}
