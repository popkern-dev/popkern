//! A boolean keeps the imprecision ceiling its operands were built to.
//!
//! A uniform scale-up grows the model's imprecision budget with the body
//! (`BRep::transform`): what construction left between an edge and its
//! faces grows a thousandfold at a thousand. A boolean of such bodies
//! stamped with the bare default ceiling (100 linear) declined a piece of
//! an operand's 9.1e-4 rim edge measured at 9.5e-4 (generated corpus
//! scene 1440); the result now carries its operands' ceiling.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn unit_box(at: Point3) -> BRep {
    let profile = [
        at,
        Point3::new(at.x + 1.0, at.y, at.z),
        Point3::new(at.x + 1.0, at.y + 1.0, at.z),
        Point3::new(at.x, at.y + 1.0, at.z),
    ];
    BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .expect("box")
}

#[test]
fn a_boolean_of_scaled_bodies_keeps_their_imprecision_ceiling() {
    let scale = Transform3::scale_uniform(1000.0);
    let a = unit_box(Point3::new(0.0, 0.0, 0.0))
        .transform(&scale)
        .expect("scaled a");
    let b = unit_box(Point3::new(0.5, 0.3, 0.2))
        .transform(&scale)
        .expect("scaled b");
    let ceiling = a.max_imprecision.max(b.max_imprecision);
    // The scale grew the operands' budget past the boolean's default.
    let default = Tolerance::default().linear() * 100.0;
    assert!(ceiling > 10.0 * default, "operand ceiling {ceiling:e}");
    for op in [
        BoolOperation::Union,
        BoolOperation::Intersection,
        BoolOperation::Subtraction,
    ] {
        let result = boolean(op, a.clone(), b.clone(), Tolerance::default(), None).expect("served");
        assert!(
            result.max_imprecision >= ceiling,
            "{op:?}: result ceiling {:e} under the operands' {ceiling:e}",
            result.max_imprecision
        );
    }
}
