//! Material-sense resolution pins.

use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn unit_cube_material_sense_is_determinate_and_flips_with_body_orientation() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let seed = brep.faces.keys().next().unwrap();
    let sense = brep
        .component_material_sense(seed)
        .expect("cube sense determinate");

    let mut reversed = brep.clone();
    reversed.reverse_body_orientation();
    let seed_r = reversed.faces.keys().next().unwrap();
    let sense_r = reversed
        .component_material_sense(seed_r)
        .expect("reversed cube sense determinate");
    assert_ne!(
        sense, sense_r,
        "reversing the body orientation must flip the material sense"
    );
}
